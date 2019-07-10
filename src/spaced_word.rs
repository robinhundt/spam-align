use std::collections::{HashMap, HashSet};
use std::convert::TryFrom;
use std::error::Error;
use std::str::FromStr;

use itertools::Itertools;
use rand::prelude::*;
use serde::de::Unexpected::Seq;

use crate::data_loaders::{Alignment, Sequence, SequencePosition};
use crate::score::score_prot_pairwise;

type BoxError = Box<dyn Error>;

#[derive(PartialEq, Eq, Debug, Copy, Clone, Hash)]
pub enum Position {
    Match,
    DontCare,
}

impl Position {
    pub fn is_match(self) -> bool {
        self == Position::Match
    }

    pub fn is_dont_care(self) -> bool {
        self == Position::DontCare
    }
}

impl TryFrom<&u8> for Position {
    type Error = Box<dyn Error>;

    fn try_from(value: &u8) -> Result<Self, Self::Error> {
        let converted = match value {
            0 => Position::DontCare,
            1 => Position::Match,
            _ => Err("Invalid token in provided pattern, only 0 and 1 allowed")?,
        };
        Ok(converted)
    }
}

#[derive(Debug, PartialEq, Eq, Hash)]
pub struct Pattern {
    weight: usize,
    positions: Vec<Position>,
}

pub type Word = Vec<u8>;

#[allow(clippy::len_without_is_empty)]
impl Pattern {
    pub fn new(binary_pattern: &[u8]) -> Result<Self, BoxError> {
        let positions = binary_pattern
            .iter()
            .map(Position::try_from)
            .collect::<Result<Vec<Position>, _>>()?;
        let weight = positions
            .iter()
            .filter(|&&pos| pos == Position::Match)
            .count();
        Ok(Pattern { weight, positions })
    }

    pub fn weight(&self) -> usize {
        self.weight
    }

    pub fn positions(&self) -> &[Position] {
        &self.positions
    }

    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// returns the spaced word and the word of the don't care positions
    pub fn match_slice(&self, seq_slice: &[u8]) -> Result<(Word, Word), BoxError> {
        if seq_slice.len() < self.len() {
            Err("Slice must be longer than pattern")?;
        }

        Ok(self.positions().iter().zip(seq_slice.iter()).fold(
            (
                Vec::with_capacity(self.weight()),
                Vec::with_capacity(self.len() - self.weight()),
            ),
            |(mut match_vec, mut dont_care_vec), (pattern_pos, seq_pos)| {
                match pattern_pos {
                    Position::Match => match_vec.push(*seq_pos),
                    Position::DontCare => dont_care_vec.push(*seq_pos),
                };
                (match_vec, dont_care_vec)
            },
        ))
    }
}

impl FromStr for Pattern {
    type Err = Box<dyn Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Pattern::new(
            &s.chars()
                .map(|ch| ch.to_digit(10).expect("Pattern must contain only 0 and 1") as u8)
                .collect::<Vec<_>>(),
        )
    }
}

pub fn generate_random_patterns(
    count: usize,
    min_len: usize,
    max_len: usize,
    min_weight: usize,
    max_weight: usize,
) -> Vec<Pattern> {
    if min_len < max_weight {
        panic!("min_len: {} is less than max_weight: {}", min_len, max_len);
    }
    let mut rng = StdRng::seed_from_u64(13546515);
    (0..count)
        .map(|_| {
            let len = rng.gen_range(min_len, max_len + 1);
            let weight = rng.gen_range(min_weight, max_weight + 1);
            let mut match_positions: HashSet<usize> = HashSet::with_capacity(weight);
            while match_positions.len() < weight {
                let match_pos = rng.gen_range(0, len);
                match_positions.insert(match_pos);
            }
            Pattern {
                weight,
                positions: (0..len)
                    .map(|pos| match match_positions.contains(&pos) {
                        true => Position::Match,
                        false => Position::DontCare,
                    })
                    .collect(),
            }
        })
        .collect()
}

#[derive(Clone, Debug)]
pub struct SpacedWord {
    pub seq: usize,
    pub pos_in_seq: usize,
    pub dont_care_word: Word,
}

#[derive(Clone, Debug)]
pub struct SpacedWordMatch {
    pub fst_word: SpacedWord,
    pub snd_word: SpacedWord,
    pub score: i32,
}

#[derive(Debug, Clone)]
pub enum MatchConsistency {
    Consistent,
    Inconsistent,
    PartiallyConsistent,
}

impl SpacedWordMatch {
    pub fn is_consistent(&self, spaced_word_len: usize, alignment: &Alignment) -> MatchConsistency {
        let fst_range = self.fst_word.pos_in_seq..self.fst_word.pos_in_seq + spaced_word_len;
        let snd_range = self.snd_word.pos_in_seq..self.snd_word.pos_in_seq + spaced_word_len;

        let positions = fst_range
            .zip(snd_range)
            .map(|(pos1, pos2)| {
                let pos1 = SequencePosition {
                    seq: self.fst_word.seq,
                    pos: pos1,
                };
                let pos2 = SequencePosition {
                    seq: self.snd_word.seq,
                    pos: pos2,
                };
                (pos1, pos2)
            })
            .collect_vec();

        let all_positions_aligned = || {
            positions
                .iter()
                .all(|(pos1, pos2)| alignment.pos_aligned(pos1, pos2))
        };

        let partial_positions_aligned = || {
            positions
                .iter()
                .any(|(pos1, pos2)| alignment.pos_aligned(pos1, pos2))
        };

        if all_positions_aligned() {
            MatchConsistency::Consistent
        } else if partial_positions_aligned() {
            MatchConsistency::PartiallyConsistent
        } else {
            MatchConsistency::Inconsistent
        }
    }
}

pub fn find_word_matches(
    pattern: &Pattern,
    sequences: &[Sequence],
) -> impl Iterator<Item = SpacedWordMatch> {
    let match_buckets = find_word_match_buckets(pattern, sequences);
    match_buckets
        .into_iter()
        .map(|(_, spaced_words)| best_match(spaced_words))
        .filter(Option::is_some)
        .map(Option::unwrap)
}

fn best_match(spaced_words: Vec<SpacedWord>) -> Option<SpacedWordMatch> {
    let spaced_word_match = spaced_words
        .into_iter()
        .tuple_combinations::<(SpacedWord, SpacedWord)>()
        .filter(|(word1, word2)| word1.seq != word2.seq)
        .max_by_key(|(word1, word2)| {
            score_prot_pairwise(&word1.dont_care_word, &word2.dont_care_word)
        });

    spaced_word_match.map(|(fst_word, snd_word)| {
        let score = score_prot_pairwise(&fst_word.dont_care_word, &snd_word.dont_care_word);
        SpacedWordMatch {
            fst_word,
            snd_word,
            score,
        }
    })
}

pub fn find_word_match_buckets(
    pattern: &Pattern,
    sequences: &[Sequence],
) -> HashMap<Word, Vec<SpacedWord>> {
    let mut spaced_words_in_seqs: HashMap<Word, Vec<SpacedWord>> = HashMap::new();

    let spaced_word_iter = sequences
        .into_iter()
        .enumerate()
        .flat_map(|(seq_idx, seq)| {
            word_matches_in_single_sequence(pattern, &seq.data).map(
                move |(match_word, dont_care_word, pos)| {
                    (
                        match_word,
                        SpacedWord {
                            seq: seq_idx,
                            pos_in_seq: pos,
                            dont_care_word,
                        },
                    )
                },
            )
        });
    for (match_slice, seq_pos) in spaced_word_iter {
        if let Some(word_matches) = spaced_words_in_seqs.get_mut(&match_slice) {
            word_matches.push(seq_pos);
        } else {
            spaced_words_in_seqs.insert(match_slice, vec![seq_pos]);
        }
    }

    spaced_words_in_seqs
}

fn word_matches_in_single_sequence<'b, 'a: 'b>(
    pattern: &'a Pattern,
    sequence: &'b [u8],
) -> impl Iterator<Item = (Word, Word, usize)> + 'b {
    sequence
        .windows(pattern.len())
        .enumerate()
        .map(move |(start_pos, slice)| {
            let (match_word, dont_care_word) = pattern
                .match_slice(slice)
                .expect("Bug in word_matches_in_single_sequence");
            (match_word, dont_care_word, start_pos)
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_random_patterns() {
        dbg!(generate_random_patterns(10, 6, 12, 3, 4));
        panic!();
    }
}
