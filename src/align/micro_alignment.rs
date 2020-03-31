use crate::spaced_word::{MatchWord, Pattern};
use crate::Sequences;
use fxhash::{FxHashMap, FxHashSet};
use itertools::Itertools;
use rayon::prelude::*;
use smallvec::SmallVec;
use std::collections::BTreeSet;
use std::iter::FromIterator;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
/// Represents a multidimensional diagonal over n sequences
/// if <= 15 sequences are part of the diagonal the end sites are stored
/// inline, otherwise they spill onto the heap
pub struct MicroAlignment {
    /// start points in sequences of the diagonal
    /// uses a small vec which stores the end sites inline
    /// if there are less than 15
    pub start_sites: SmallVec<[Site; 15]>, // TODO benchmark this parameter
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
            .zip(self.start_sites.iter().skip(1))
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

pub fn compute_multi_dim_micro_alignment_information(
    patterns: &[Pattern],
    sequences: &Sequences,
    //    score_fn: fn(&[&[u8]]) -> i32,
) -> FxHashMap<usize, usize> {
    let max_seq_len = sequences
        .iter()
        .map(|seq| seq.len())
        .max()
        .expect("Can not construct diagonals from empty sequences");
    //    let (sum, count) =
    let mut dim_counts = FxHashMap::default();
    for pattern in patterns {
        generate_sorted_matches(pattern, sequences, max_seq_len)
            .into_iter()
            .group_by(|pattern_match| pattern_match.key)
            .into_iter()
            .for_each(|(_, matches)| {
                let mut matches: SmallVec<[Match; 8]> = SmallVec::from_iter(matches);
                matches.sort_unstable_by_key(|word_match| word_match.start_site.seq);
                let dim = matches
                    .into_iter()
                    .group_by(|pattern_match| pattern_match.start_site.seq)
                    .into_iter()
                    .count();
                if dim == 1 {
                    return;
                }
                *dim_counts.entry(dim).or_insert(0) += 1;
            });
    }
    dim_counts
    //        .fold((0, 0), |(sum, count), match_dim| {
    //            (sum + match_dim, count + 1)
    //        });
    //    sum as f64 / count as f64
}

pub fn construct_micro_alignments_from_patterns<'a>(
    patterns: &'a [Pattern],
    sequences: &'a Sequences,
    score_fn: fn(&[&[u8]]) -> i32,
    one_to_one_mapping: bool,
    //    max_n: Option<usize>,
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
                .flat_map(
                    move |(_, match_group)| -> Box<dyn Iterator<Item = ScoredMicroAlignment>> {
                        let match_group: SmallVec<[Match; 8]> = SmallVec::from_iter(match_group);
                        let mut scored_combinations = generate_combinations(match_group)
                            .map(|combination| {
                                score_match_combination(score_fn, combination, sequences, pattern)
                            })
                            .collect_vec();
                        if one_to_one_mapping {
                            scored_combinations = generate_one_to_one_mapping(scored_combinations);
                            Box::new(scored_combinations.into_iter())
                        } else {
                            Box::new(scored_combinations.into_iter())
                        }
                    },
                )
                .collect_vec()
                .into_iter()
        })
}

fn generate_sorted_matches(
    pattern: &Pattern,
    sequences: &Sequences,
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

fn generate_combinations(match_group: SmallVec<[Match; 8]>) -> impl Iterator<Item = Vec<Match>> {
    //    match_group.sort_unstable_by_key(|word_match| word_match.start_site.seq);
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
) -> ScoredMicroAlignment {
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
    ScoredMicroAlignment {
        micro_alignment: MicroAlignment {
            start_sites,
            k: pattern.len(),
        },
        score: msa_score,
    }
}

fn generate_one_to_one_mapping(mut data: Vec<ScoredMicroAlignment>) -> Vec<ScoredMicroAlignment> {
    // TODO this method seems overly complex and could likely be improved
    data.sort_unstable_by_key(|ma| {
        BTreeSet::from_iter(
            ma.micro_alignment
                .start_sites
                .iter()
                .map(|start_site| start_site.seq),
        )
    });
    let ma_in_same_seqs = data.into_iter().group_by(|ma| {
        FxHashSet::from_iter(
            ma.micro_alignment
                .start_sites
                .iter()
                .map(|start_site| start_site.seq),
        )
    });
    let mut one_to_one_mapping = vec![];
    for (_, micro_alignments) in ma_in_same_seqs.into_iter() {
        let mut micro_alignments = micro_alignments.collect_vec();
        micro_alignments.sort_unstable_by_key(|ma: &ScoredMicroAlignment| ma.score);
        while let Some(micro_alignment) = micro_alignments.pop() {
            let diag_start_sites =
                FxHashSet::from_iter(micro_alignment.micro_alignment.start_sites.iter());
            // TODO here the retain method could be used
            micro_alignments = micro_alignments
                .into_iter()
                .filter(|el| {
                    let el_start_sites =
                        FxHashSet::from_iter(el.micro_alignment.start_sites.iter());
                    diag_start_sites.is_disjoint(&el_start_sites)
                })
                .collect();
            one_to_one_mapping.push(micro_alignment);
        }
    }

    one_to_one_mapping
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_loaders::balibase::parse_xml_file;
    use crate::score::score_prot_msa;
    use crate::spaced_word::read_patterns_from_file;
    use std::error::Error;

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

        assert_eq!(generate_combinations(match_group).collect_vec(), expected)
    }

    #[test]
    fn one_to_one_mapping() {
        macro_rules! s {
            ($seq:expr;$pos:expr) => {
                Site {
                    seq: $seq,
                    pos: $pos
                }
            };
            ($len:expr, $score:expr, [$($seq:expr;$pos:expr),+]) => {
                {
                    let start_sites = SmallVec::from_vec(
                        vec![$(s!($seq;$pos),)+]
                    );
                    ScoredMicroAlignment {
                        micro_alignment: MicroAlignment {
                            start_sites: start_sites,
                            k: $len
                        },
                        score: $score
                    }
                }
            };
        }

        let scored_micro_alignment = vec![
            s!(5, 10, [1;1, 2;3, 4;5]),
            s!(5, 15, [1;1, 2;30, 4;50]),
            s!(5, 5, [6;1, 7;3]),
            s!(5, 8, [6;10, 7;10]),
        ];

        let expected_mapping = FxHashSet::from_iter(
            vec![
                s!(5, 15, [1;1, 2;30, 4;50]),
                s!(5, 5, [6;1, 7;3]),
                s!(5, 8, [6;10, 7;10]),
            ]
            .into_iter(),
        );
        assert_eq!(
            FxHashSet::from_iter(generate_one_to_one_mapping(scored_micro_alignment).into_iter()),
            expected_mapping
        );

        let scored_micro_alignment = vec![s!(1, 10, [1;1, 2;1]), s!(1, 8, [2;1, 3;1])];

        let expected_mapping = scored_micro_alignment.clone();

        assert_eq!(
            generate_one_to_one_mapping(scored_micro_alignment),
            expected_mapping
        );
    }
}
