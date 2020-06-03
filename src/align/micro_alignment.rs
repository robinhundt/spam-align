use std::iter::FromIterator;

use itertools::Itertools;
use smallvec::SmallVec;

use crate::align::Strategy;
use crate::spaced_word::{MatchWord, Pattern};
use crate::Sequence;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
/// Represents a multidimensional diagonal over n sequences
/// if <= 15 sequences are part of the diagonal the end sites are stored
/// inline, otherwise they spi2 onto the heap
pub struct MicroAlignment {
    /// start points in sequences of the diagonal
    /// uses a small vec which stores the end sites inline
    /// if there are less than 15
    pub start_sites: SmallVec<[Site; 2]>, // TODO benchmark this parameter
    /// diagonal length
    pub k: usize,
}

impl MicroAlignment {
    pub fn site_iter(&self) -> impl Iterator<Item = Site> + '_ {
        self.start_sites.iter().flat_map(move |start_site| {
            (start_site.pos..start_site.pos + self.k).map(move |pos| Site {
                seq: start_site.seq,
                pos,
            })
        })
    }

    pub fn site_pair_iter(&self) -> impl Iterator<Item = (Site, Site)> + '_ {
        self.start_sites
            .iter()
            .tuple_combinations()
            .flat_map(move |(site_a, site_b)| {
                (site_a.pos..site_a.pos + self.k)
                    .zip(site_b.pos..site_b.pos + self.k)
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
pub struct ScoredMicroAlignment {
    pub micro_alignment: MicroAlignment,
    pub score: i32,
}

impl ScoredMicroAlignment {
    pub fn site_pair_iter(&self) -> impl Iterator<Item = (Site, Site)> + '_ {
        self.micro_alignment.site_pair_iter()
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash, Ord, PartialOrd, Default)]
pub struct Site {
    pub seq: usize,
    pub pos: usize,
}

#[derive(Clone, Debug, Copy, PartialEq)]
pub struct Match {
    key: MatchWord,
    start_site: Site,
}

pub fn construct_2dim_micro_alignments<'a>(
    patterns: &'a [Pattern],
    sequences: &'a [Sequence],
    score_fn: fn(&[&[u8]]) -> i32,
) -> impl Iterator<Item = ScoredMicroAlignment> + 'a {
    let max_seq_len = sequences
        .iter()
        .map(|seq| seq.len())
        .max()
        .expect("Can not construct diagonals from empty sequences");
    patterns
        .iter()
        .map(move |pattern| {
            (
                pattern,
                generate_sorted_matches(pattern, sequences, max_seq_len),
            )
        })
        .flat_map(move |(pattern, matches)| {
            matches
                .into_iter()
                .group_by(|pattern_match| pattern_match.key)
                .into_iter()
                .flat_map(move |(_, match_group)| {
                    let match_group: SmallVec<[Match; 8]> = SmallVec::from_iter(match_group);
                    generate_2_dim_combinations(match_group).map(move |combination| {
                        score_match_combination(score_fn, combination, sequences, pattern)
                    })
                })
                .collect_vec()
                .into_iter()
        })
}

pub fn construct_dyn_dim_micro_alignments<'a>(
    patterns: &'a [Pattern],
    sequences: &'a [Sequence],
    score_fn: fn(&[&[u8]]) -> i32,
) -> impl Iterator<Item = ScoredMicroAlignment> + 'a {
    let max_seq_len = sequences
        .iter()
        .map(|seq| seq.len())
        .max()
        .expect("Can not construct diagonals from empty sequences");
    patterns
        .iter()
        .map(move |pattern| {
            (
                pattern,
                generate_sorted_matches(pattern, sequences, max_seq_len),
            )
        })
        .flat_map(move |(pattern, matches)| {
            matches
                .into_iter()
                .group_by(|pattern_match| pattern_match.key)
                .into_iter()
                .flat_map(move |(_, match_group)| {
                    let match_group: SmallVec<[Match; 8]> = SmallVec::from_iter(match_group);
                    generate_dyn_dim_combinations(match_group).map(move |combination| {
                        score_match_combination(score_fn, combination, sequences, pattern)
                    })
                })
                .collect_vec()
                .into_iter()
        })
}

fn generate_sorted_matches(
    pattern: &Pattern,
    sequences: &[Sequence],
    max_seq_len: usize,
) -> Vec<Match> {
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
    pattern_matches.sort_unstable_by_key(|pattern_match| pattern_match.key);
    pattern_matches
}

fn generate_2_dim_combinations(
    match_group: SmallVec<[Match; 8]>,
) -> impl Iterator<Item = Vec<Match>> {
    match_group
        .into_iter()
        .combinations(2)
        .filter(|combination| combination[0].start_site.seq != combination[1].start_site.seq)
}

fn generate_dyn_dim_combinations(
    mut match_group: SmallVec<[Match; 8]>,
) -> Box<dyn Iterator<Item = Vec<Match>>> {
    match_group.sort_unstable_by_key(|word_match| word_match.start_site.seq);
    let match_group = match_group;
    if match_group
        .iter()
        .group_by(|word_match| word_match.start_site.seq)
        .into_iter()
        .map(|a| a.1.count() as u128)
        .product::<u128>()
        > 100
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
    Box::new(
        match_group
            .clone()
            .into_iter()
            .group_by(|word_match| word_match.start_site.seq)
            .into_iter()
            .map(|(_, group)| group.collect_vec())
            .multi_cartesian_product()
            .filter(|combination| combination.len() > 1)
            .chain(
                match_group
                    .into_iter()
                    .combinations(2)
                    .filter(|combination| {
                        combination[0].start_site.seq != combination[1].start_site.seq
                    }),
            ),
    )
}

fn score_match_combination(
    score_fn: fn(&[&[u8]]) -> i32,
    combination: Vec<Match>,
    sequences: &[Sequence],
    pattern: &Pattern,
) -> ScoredMicroAlignment {
    let msa: SmallVec<[&[u8]; 2]> = SmallVec::from_iter(combination.iter().map(|pattern_match| {
        let start_site = pattern_match.start_site;
        &sequences[start_site.seq].data[start_site.pos..start_site.pos + pattern.len()]
    }));
    let msa_score = score_fn(&msa[..]);
    let start_sites: SmallVec<[Site; 2]> = SmallVec::from_iter(
        combination
            .into_iter()
            .map(|pattern_match| pattern_match.start_site),
    );
    ScoredMicroAlignment {
        micro_alignment: MicroAlignment {
            start_sites,
            k: pattern.len(),
        },
        score: msa_score,
    }
}

impl From<bool> for Strategy {
    fn from(val: bool) -> Self {
        if val {
            Strategy::DynDim
        } else {
            Strategy::TwoDim
        }
    }
}

#[cfg(test)]
mod tests {
    use std::error::Error;

    use super::*;

    #[test]
    fn test_construct_micro_alignments_from_patterns() -> Result<(), Box<dyn Error>> {
        //        dbg!("Running test");
        //        let alignment = parse_xml_file("data/bb3_release/RV30/BBS30004.xml").unwrap();
        //        let sequences = Sequences::new(&alignment);
        //        let patterns = read_patterns_from_file("sample.pat")?;
        //        let diags = construct_micro_alignments_from_patterns(
        //            &patterns[..],
        //            &sequences,
        //            score_prot_msa,
        //            false,
        //        );
        //        dbg!(diags.len());
        //        dbg!(
        //            diags
        //                .iter()
        //                .fold(0, |acc, diag| acc + diag.micro_alignment.start_sites.len())
        //                as f64
        //                / diags.len() as f64
        //        );
        //        assert!(diags.len() > 0);
        Ok(())
    }

    #[test]
    fn dim_2_site_pair_iter() {
        let diag = MicroAlignment {
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
        let diag = MicroAlignment {
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
            (Site { seq: 1, pos: 3 }, Site { seq: 4, pos: 4 }),
            (Site { seq: 1, pos: 4 }, Site { seq: 4, pos: 5 }),
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
        //        let expected = vec![vec![m!(1), m!(2), vec![m!(1), m!(2), m!(3)]];
        // if dim == 2
        let expected = vec![
            vec![m!(1), m!(2)],
            vec![m!(1), m!(2)],
            vec![m!(1), m!(3)],
            vec![m!(2), m!(3)],
            vec![m!(2), m!(3)],
        ];

        // if dynaimc max dim
        // let expected = vec![
        //     vec![m!(1), m!(2)],
        //     vec![m!(1), m!(2)],
        //     vec![m!(1), m!(3)],
        //     vec![m!(2), m!(3)],
        //     vec![m!(2), m!(3)],
        // ];

        assert_eq!(
            generate_2_dim_combinations(match_group).collect_vec(),
            expected
        )
    }
}
