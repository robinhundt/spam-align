use std::collections::HashSet;
use std::convert::TryFrom;
use std::str::FromStr;

use anyhow::{Error, Result};
use fxhash::{hash, FxHashMap, FxHashSet};
use itertools::Itertools;
use rand::prelude::*;

use crate::data_loaders::Sequence;
use crate::score::score_prot_pairwise;
use std::fs::File;
use std::io::Read;
use std::ops::Not;
use std::path::Path;

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
    type Error = Error;

    fn try_from(value: &u8) -> Result<Self, Self::Error> {
        let converted = match value {
            0 => Position::DontCare,
            1 => Position::Match,
            _ => Err(anyhow!(
                "Invalid token in provided pattern, only 0 and 1 allowed"
            ))?,
        };
        Ok(converted)
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct Pattern {
    weight: usize,
    positions: Vec<Position>,
}

pub type Word = Vec<u8>;

#[allow(clippy::len_without_is_empty)]
impl Pattern {
    const MAX_WEIGHT: usize = 12;

    pub fn new(binary_pattern: &[u8]) -> Result<Self> {
        let positions = binary_pattern
            .iter()
            .map(Position::try_from)
            .collect::<Result<Vec<Position>, _>>()?;
        let weight = positions
            .iter()
            .filter(|&&pos| pos == Position::Match)
            .count();
        if weight > Self::MAX_WEIGHT {
            Err(anyhow!(
                "Pattern max weight is {}, but pattern \"{:#?}\" has weight: {}",
                Self::MAX_WEIGHT,
                binary_pattern,
                weight,
            ))?;
        }
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
    pub fn _match_slice(&self, seq_slice: &[u8]) -> Result<(Word, Word)> {
        if seq_slice.len() < self.len() {
            Err(anyhow!("Slice must be longer than pattern"))?;
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

    pub fn match_slice(&self, seq_slice: &[u8]) -> MatchWord {
        if seq_slice.len() != self.len() {
            panic!("Slice must be as exactly as long as pattern");
        }

        let match_buf = self
            .positions()
            .iter()
            .zip(seq_slice.iter())
            .fold(
                ([std::u8::MAX; Self::MAX_WEIGHT], 0_usize),
                |(mut match_buf, mut match_cnt), (pattern_pos, &seq_pos)| {
                    match pattern_pos {
                        Position::Match => {
                            match_buf[match_cnt] = seq_pos;
                            match_cnt += 1;
                        }
                        Position::DontCare => (),
                    };
                    (match_buf, match_cnt)
                },
            )
            .0;
        MatchWord::from(match_buf)
    }
}

impl FromStr for Pattern {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Pattern::new(
            &s.chars()
                .map(|ch| ch.to_digit(10).expect("Pattern must contain only 0 and 1") as u8)
                .collect::<Vec<_>>(),
        )
    }
}

pub fn read_patterns_from_file(path: impl AsRef<Path>) -> Result<Vec<Pattern>> {
    let mut file = File::open(path)?;
    let mut buf = String::new();
    file.read_to_string(&mut buf)?;
    buf.split("\n")
        .filter(|s| s.len() != 0 && s.trim_start().starts_with("#").not())
        .map(|s| s.parse())
        .collect()
}

#[derive(Debug, Ord, PartialOrd, Eq, PartialEq, Copy, Clone)]
pub struct MatchWord(u64);

// the [u8; 12] must contain protein char (e.g. 'A', 'F', 'M') as u8
// if the number of match positions is less than 12, the unused postions
// should be filled with std::u8::MAX
impl From<[u8; 12]> for MatchWord {
    fn from(arr: [u8; 12]) -> Self {
        let mut res = 0_u64;

        for (count, protein) in arr.iter().enumerate() {
            let encoded = encode_protein_u8(*protein);
            res |= (encoded as u64) << count * 5;
        }
        MatchWord(res)
    }
}

impl From<MatchWord> for [u8; 12] {
    fn from(word: MatchWord) -> Self {
        let a = word.0;
        let mut res = [std::u8::MAX; 12];
        for count in 0..12 {
            let mask = 0b11111_u64 << count * 5;
            let mut encoded = (a & mask) >> count * 5;
            if encoded == 0b11111_u64 {
                encoded = std::u8::MAX as u64;
            }
            let decoded = decode_protein_u64(encoded);
            res[count] = decoded;
        }
        res
    }
}

fn encode_protein_u8(a: u8) -> u64 {
    let mut a = a as char;
    a.make_ascii_uppercase();
    let a = a as u8;
    if a == b'Y' {
        23
    } else if a == b'Z' {
        24
    } else if a == b'X' {
        25
    } else if a == b'*' {
        26
    } else if a == std::u8::MAX {
        a as u64
    } else {
        (a - 65) as u64
    }
}

fn decode_protein_u64(a: u64) -> u8 {
    if a == 23 {
        b'Y'
    } else if a == 24 {
        b'Z'
    } else if a == 25 {
        b'X'
    } else if a == 26 {
        b'*'
    } else if a == std::u8::MAX as u64 {
        a as u8
    } else {
        (a + 65) as u8
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
            let mut match_positions: FxHashSet<usize> = HashSet::default();
            match_positions.reserve(weight);
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

#[derive(Clone, Debug, Hash)]
pub struct SpacedWord {
    pub seq: usize,
    pub pos_in_seq: usize,
    pub dont_care_word: Word,
}

#[derive(Clone, Debug, Hash)]
pub struct SpacedWordMatch {
    pub fst_word: SpacedWord,
    pub snd_word: SpacedWord,
    pub score: i32,
}

//TODO provide a way to find matches for mutliple patterns at the same time,
// should be much faster
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

// TODO evaluate if  this is really what is supposed to be happening
fn best_match(spaced_words: Vec<SpacedWord>) -> Option<SpacedWordMatch> {
    let spaced_word_match = spaced_words
        .into_iter()
        .tuple_combinations::<(SpacedWord, SpacedWord)>()
        .filter(|(word1, word2)| word1.seq != word2.seq)
        .fold((vec![], std::i32::MIN), |mut acc, (word1, word2)| {
            let score = score_prot_pairwise(&word1.dont_care_word, &word2.dont_care_word);
            if score > acc.1 {
                acc.0.clear();
                acc.0.push((word1, word2));
                acc.1 = score;
            } else if score == acc.1 {
                acc.0.push((word1, word2));
            }
            acc
        })
        .0
        .into_iter()
        .max_by_key(|word_tuple| hash(word_tuple));

    spaced_word_match.map(|(fst_word, snd_word)| {
        let score = score_prot_pairwise(&fst_word.dont_care_word, &snd_word.dont_care_word);
        SpacedWordMatch {
            fst_word,
            snd_word,
            score,
        }
    })
}
//TODO provide a way to find matches for mutliple patterns at the same time,
// should be much faster
pub fn find_word_match_buckets(
    pattern: &Pattern,
    sequences: &[Sequence],
) -> FxHashMap<Word, Vec<SpacedWord>> {
    let mut spaced_words_in_seqs: FxHashMap<Word, Vec<SpacedWord>> = FxHashMap::default();

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

//TODO provide a way to find matches for mutliple patterns at the same time,
// should be much faster
fn word_matches_in_single_sequence<'b, 'a: 'b>(
    pattern: &'a Pattern,
    sequence: &'b [u8],
) -> impl Iterator<Item = (Word, Word, usize)> + 'b {
    sequence
        .windows(pattern.len())
        .enumerate()
        .map(move |(start_pos, slice)| {
            let (match_word, dont_care_word) = pattern
                ._match_slice(slice)
                .expect("Bug in word_matches_in_single_sequence");
            (match_word, dont_care_word, start_pos)
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::mem::size_of;

    #[test]
    fn test_generate_random_patterns() {
        dbg!(generate_random_patterns(10, 6, 12, 3, 4));
    }

    #[test]
    fn match_word_encoding() {
        let max = std::u8::MAX;
        let a = [
            b'L', b'N', b'A', b'F', b'M', b'L', b'Y', b'M', b'K', b'E', max, max,
        ];
        let b = [
            b'K', b'K', b'K', b'R', b'K', b'R', b'E', b'K', max, max, max, max,
        ];
        let c = [
            b'K', b'K', b'K', b'R', b'K', b'R', b'E', b'K', b'K', b'R', b'E', b'K',
        ];

        assert_eq!(<[u8; 12]>::from(MatchWord::from(a)), a);
        assert_eq!(<[u8; 12]>::from(MatchWord::from(b)), b);
        assert_eq!(<[u8; 12]>::from(MatchWord::from(c)), c);
    }
}
