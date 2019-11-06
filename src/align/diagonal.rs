use crate::spaced_word::{MatchWord, Pattern};
use crate::Sequences;
use fxhash::{FxHashMap, FxHashSet};
use itertools::Itertools;
use ndarray::Array2;
use rayon::prelude::*;
use smallvec::SmallVec;
use std::collections::VecDeque;
use std::iter::FromIterator;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
/// Represents a multidimensional diagonal over n sequences
/// if <= 15 sequences are part of the diagonal the end sites are stored
/// inline, otherwise they spill onto the heap
pub struct Diagonal {
    /// end points (inclusive) in sequences of the diagonal
    /// uses a small vec which stores the end sites inline
    /// if there are less than 15
    pub end_sites: SmallVec<[Site; 15]>, // TODO benchmark this parameter
    /// diagonal length
    k: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ScoredDiagonal {
    diag: Diagonal,
    score: i32,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub struct Site {
    pub seq: usize,
    pub pos: usize,
}

#[derive(Clone, Debug, Copy)]
pub struct Match {
    key: MatchWord,
    end_site: Site,
}

pub fn construct_diagonals_from_patterns(
    patterns: &[Pattern],
    sequences: &Sequences,
    score_fn: fn(&[&[u8]]) -> i32,
    max_n: Option<usize>,
) -> Vec<ScoredDiagonal> {
    let max_seq_len = sequences
        .iter()
        .map(|seq| seq.len())
        .max()
        .expect("Can not construct diagonals from empty sequences");
    //    let mut buf = Vec::with_capacity(max_seq_len * sequences.len() * patterns.len());
    patterns
        .par_iter()
        .map(|pattern| {
            let mut pattern_matches: Vec<Match> = Vec::with_capacity(max_seq_len * sequences.len());
            for (idx, sequence) in sequences.iter().enumerate() {
                for (slice_count, slice) in sequence.data.windows(pattern.len()).enumerate() {
                    let match_word = pattern.match_slice(slice);
                    let end_site = Site {
                        seq: idx,
                        // - 1 because end sites are inclusive
                        pos: slice_count + pattern.len() - 1,
                    };
                    pattern_matches.push(Match {
                        key: match_word,
                        end_site,
                    });
                }
            }
            pattern_matches.sort_by_cached_key(|pattern_match| pattern_match.key);
            (pattern, pattern_matches)
        })
        .flat_map(|(pattern, matches)| {
            matches
                .into_iter()
                .group_by(|pattern_match| pattern_match.key)
                .into_iter()
                .flat_map(|(key, match_group)| {
                    let match_group: SmallVec<[Match; 32]> = SmallVec::from_iter(match_group);
                    let max_n = max_n.unwrap_or(match_group.len());
                    (2..=max_n)
                        .flat_map(|match_dim| {
                            // TODO combinations creates A LOT of vectors
                            // this is bad in such a hot path
                            match_group
                                .iter()
                                .combinations(max_n)
                                .filter(|combinations| {
                                    // TODO this uses a hashset internally and should probably be
                                    // optimized
                                    // combinations.0.end_site.seq != combinations.1.end_site.seq
                                    combinations
                                        .iter()
                                        .unique_by(|pattern_match| pattern_match.end_site.seq)
                                        .count()
                                        == combinations.len()
                                })
                                .map(|combinations| {
                                    //                                    let combinations = tuple_to_arr(combinations);
                                    let msa: SmallVec<[&[u8]; 15]> = SmallVec::from_iter(
                                        combinations.iter().map(|pattern_match| {
                                            sequences.get_site_slice(
                                                pattern_match.end_site,
                                                pattern.len(),
                                            )
                                        }),
                                    );
                                    let msa_score = score_fn(&msa[..]);
                                    let end_sites: SmallVec<[Site; 15]> = SmallVec::from_iter(
                                        combinations
                                            .into_iter()
                                            .map(|pattern_match| pattern_match.end_site),
                                    );
                                    ScoredDiagonal {
                                        diag: Diagonal {
                                            end_sites,
                                            k: pattern.len(),
                                        },
                                        score: msa_score,
                                    }
                                })
                        })
                        .collect_vec()
                        .into_iter()
                })
                .collect_vec()
                .into_par_iter()
        })
        .collect()
}

fn tuple_to_arr((a, b): (&Match, &Match)) -> [Match; 2] {
    [*a, *b]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_loaders::balibase::parse_xml_file;
    use crate::score::score_prot_msa;
    use crate::spaced_word::read_patterns_from_file;
    use std::error::Error;

    #[test]
    fn test_construct_diagonals_from_patterns() -> Result<(), Box<dyn Error>> {
        let alignment = parse_xml_file("data/bb3_release/RV30/BBS30004.xml").unwrap();
        let sequences = Sequences::new(&alignment);
        let patterns = read_patterns_from_file("sample.pat")?;
        let diags =
            construct_diagonals_from_patterns(&patterns[..], &sequences, score_prot_msa, Some(2));
        dbg!(diags.len());
        dbg!(
            diags
                .iter()
                .fold(0, |acc, diag| acc + diag.diag.end_sites.len()) as f64
                / diags.len() as f64
        );
        panic!()
    }
}
