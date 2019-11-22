use crate::spaced_word::{MatchWord, Pattern};
use crate::Sequences;
use fxhash::{FxHashMap, FxHashSet};
use itertools::{Group, GroupBy, Itertools};
use ndarray::Array2;
use rayon::prelude::*;
use smallvec::SmallVec;
use std::cmp::min;
use std::collections::VecDeque;
use std::iter::FromIterator;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
/// Represents a multidimensional diagonal over n sequences
/// if <= 15 sequences are part of the diagonal the end sites are stored
/// inline, otherwise they spill onto the heap
pub struct Diagonal {
    /// start points in sequences of the diagonal
    /// uses a small vec which stores the end sites inline
    /// if there are less than 15
    pub start_sites: SmallVec<[Site; 15]>, // TODO benchmark this parameter
    /// diagonal length
    pub k: usize,
}

impl Diagonal {
    pub fn site_pair_iter(&self) -> impl Iterator<Item = (Site, Site)> + '_ {
        self.start_sites
            .iter()
            .zip(self.start_sites.iter().skip(1))
            .flat_map(move |(site_a, site_b)| {
                (site_a.pos..site_a.pos + self.k)
                    .zip((site_b.pos..site_b.pos + self.k))
                    .map(move |(pos_a, pos_b)| {
                        (
                            Site {
                                seq: site_a.seq,
                                pos: pos_a,
                            },
                            Site {
                                seq: site_b.seq,
                                pos: pos_b,
                            },
                        )
                    })
            })
    }

    pub fn site_pair_par_iter(&self) -> impl ParallelIterator<Item = (Site, Site)> + '_ {
        self.start_sites
            .par_iter()
            .zip(self.start_sites.par_iter().skip(1))
            .flat_map(move |(site_a, site_b)| {
                (site_a.pos..site_a.pos + self.k)
                    .into_par_iter()
                    .zip((site_b.pos..site_b.pos + self.k).into_par_iter())
                    .map(move |(pos_a, pos_b)| {
                        (
                            Site {
                                seq: site_a.seq,
                                pos: pos_a,
                            },
                            Site {
                                seq: site_b.seq,
                                pos: pos_b,
                            },
                        )
                    })
            })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ScoredDiagonal {
    pub diag: Diagonal,
    pub score: i32,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub struct Site {
    pub seq: usize,
    pub pos: usize,
}

#[derive(Clone, Debug, Copy, PartialEq)]
pub struct Match {
    key: MatchWord,
    start_site: Site,
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
    patterns
        .par_iter()
        .map(|pattern| generate_sorted_matches(pattern, sequences, max_seq_len))
        .flat_map(|(pattern, matches)| {
            matches
                .into_iter()
                .group_by(|pattern_match| pattern_match.key)
                .into_iter()
                .flat_map(|(key, match_group)| {
                    let mut match_group: SmallVec<[Match; 8]> = SmallVec::from_iter(match_group);
                    generate_combinations(match_group).map(|combination| {
                        score_match_combination(score_fn, combination, sequences, pattern)
                    })
                })
                .collect_vec()
                .into_par_iter()
        })
        .collect()
}

fn generate_sorted_matches<'a, 'b>(
    pattern: &'a Pattern,
    sequences: &'b Sequences,
    max_seq_len: usize,
) -> (&'a Pattern, Vec<Match>) {
    let mut pattern_matches: Vec<Match> = Vec::with_capacity(max_seq_len * sequences.len());
    for (idx, sequence) in sequences.iter().enumerate() {
        for (slice_count, slice) in sequence.data.windows(pattern.len()).enumerate() {
            let match_word = pattern.match_slice(slice);
            let start_site = Site {
                seq: idx,
                pos: slice_count,
            };
            pattern_matches.push(Match {
                key: match_word,
                start_site,
            });
        }
    }
    pattern_matches.sort_by_cached_key(|pattern_match| pattern_match.key);
    (pattern, pattern_matches)
}

fn generate_combinations(
    mut match_group: SmallVec<[Match; 8]>,
) -> impl Iterator<Item = Vec<Match>> {
    match_group.sort_by_cached_key(|word_match| word_match.start_site.seq);
    //    let match_group = match_group;
    //    if match_group
    //        .iter()
    //        .group_by(|word_match| word_match.start_site.seq)
    //        .into_iter()
    //        .map(|a| a.1)
    //        .fold(1, |acc, x| acc * x.count())
    //        > 100
    {
        return Box::new(
            match_group
                .into_iter()
                .combinations(2)
                .filter(|combination| {
                    combination[0].start_site.seq != combination[1].start_site.seq
                }),
        );
    }
    //    Box::new(
    //        match_group
    //            .into_iter()
    //            .group_by(|word_match| word_match.start_site.seq)
    //            .into_iter()
    //            .map(|(_, group)| group.into_iter().collect_vec())
    //            .multi_cartesian_product()
    //            .filter(|combination| {
    //                let unique_seq_cnt = combination
    //                    .iter()
    //                    .unique_by(|match_word| match_word.start_site.seq)
    //                    .count();
    //                assert_eq!(
    //                    unique_seq_cnt,
    //                    combination.len(),
    //                    "Combination has sites with duplicate seq: {:?}",
    //                    &combination
    //                );
    //                combination.len() > 1
    //            }),
    //    )
}

fn score_match_combination(
    score_fn: fn(&[&[u8]]) -> i32,
    combination: Vec<Match>,
    sequences: &Sequences,
    pattern: &Pattern,
) -> ScoredDiagonal {
    let msa: SmallVec<[&[u8]; 15]> =
        SmallVec::from_iter(combination.iter().map(|pattern_match| {
            sequences.get_site_slice(pattern_match.start_site, pattern.len())
        }));
    let msa_score = score_fn(&msa[..]);
    let start_sites: SmallVec<[Site; 15]> = SmallVec::from_iter(
        combination
            .into_iter()
            .map(|pattern_match| pattern_match.start_site),
    );
    ScoredDiagonal {
        diag: Diagonal {
            start_sites,
            k: pattern.len(),
        },
        score: msa_score,
    }
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
        dbg!("Running test");
        let alignment = parse_xml_file("data/bb3_release/RV30/BBS30004.xml").unwrap();
        let sequences = Sequences::new(&alignment);
        let patterns = read_patterns_from_file("sample.pat")?;
        let diags =
            construct_diagonals_from_patterns(&patterns[..], &sequences, score_prot_msa, Some(2));
        dbg!(diags.len());
        dbg!(
            diags
                .iter()
                .fold(0, |acc, diag| acc + diag.diag.start_sites.len()) as f64
                / diags.len() as f64
        );
        assert!(diags.len() > 0);
        Ok(())
    }

    #[test]
    fn dim_2_site_pair_iter() {
        let diag = Diagonal {
            start_sites: SmallVec::from_vec(vec![Site { seq: 1, pos: 3 }, Site { seq: 2, pos: 5 }]),
            k: 3,
        };

        let expected = vec![
            (Site { seq: 1, pos: 3 }, Site { seq: 2, pos: 5 }),
            (Site { seq: 1, pos: 4 }, Site { seq: 2, pos: 6 }),
            (Site { seq: 1, pos: 5 }, Site { seq: 2, pos: 7 }),
        ];
        assert_eq!(diag.site_pair_iter().collect_vec(), expected);
    }

    #[test]
    fn dim_3_site_pair_iter() {
        let diag = Diagonal {
            start_sites: SmallVec::from_vec(vec![
                Site { seq: 1, pos: 3 },
                Site { seq: 2, pos: 5 },
                Site { seq: 4, pos: 4 },
            ]),
            k: 2,
        };

        let expected = vec![
            (Site { seq: 1, pos: 3 }, Site { seq: 2, pos: 5 }),
            (Site { seq: 1, pos: 4 }, Site { seq: 2, pos: 6 }),
            (Site { seq: 2, pos: 5 }, Site { seq: 4, pos: 4 }),
            (Site { seq: 2, pos: 6 }, Site { seq: 4, pos: 5 }),
        ];
        assert_eq!(diag.site_pair_iter().collect_vec(), expected);
    }

    #[test]
    fn generate_combinations_test() {
        macro_rules! m {
            ($seq:expr) => {
                Match {
                    key: MatchWord::from([100; 12]),
                    start_site: Site { seq: $seq, pos: 0 },
                }
            };
        }
        let match_group: SmallVec<[Match; 8]> =
            SmallVec::from_vec(vec![m!(1), m!(2), m!(2), m!(3)]);
        let expected = vec![vec![m!(1), m!(2), m!(3)], vec![m!(1), m!(2), m!(3)]];
        assert_eq!(generate_combinations(match_group).collect_vec(), expected)
    }
}
