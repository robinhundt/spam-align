use std::collections::HashMap;
use std::error::Error;

use bio::io::fasta::Reader;
use itertools::Itertools;
use rayon::prelude::*;

use spam_align::data_loaders::balibase::{BBAlignment, FilterXmlFile};
use spam_align::score::{score_prot_msa, score_prot_pairwise};
use spam_align::spaced_word::{find_word_match_buckets, Pattern};

fn main() -> Result<(), Box<dyn Error>> {
    //    for folder_path in &[
    //        "./data/bb3_release/RV11/",
    //        "./data/bb3_release/RV12/",
    //        "./data/bb3_release/RV20/",
    //        "./data/bb3_release/RV30/",
    //        // "./data/bb3_release/RV40/", // has no BBS files for some reason
    //        "./data/bb3_release/RV50/",
    //    ] {
    //        println!("{}", folder_path);
    //        let alignments =
    //            BBAlignment::from_xml_files_in_dir(folder_path, FilterXmlFile::CoreBlockDataOnly)?;
    //        let mut core_match_sum = 0;
    //        let mut mismatch_sum = 0;
    //        let mut partial_match_sum = 0;
    //        let mut cum_len = 0;
    //        let mut seq_cnt = 0;
    //        for alignment in alignments {
    //            // println!("Evaluating: {}", alignment.name);
    //            let unaligned_data: Vec<Vec<u8>> = alignment
    //                .sequences
    //                .iter()
    //                .map(|seq| {
    //                    seq.data
    //                        .chars()
    //                        .filter(|ch| *ch != '-' && ch.is_ascii() && !ch.is_whitespace())
    //                        .map(|ch| ch as u8)
    //                        .collect()
    //                })
    //                .collect();
    //            cum_len += unaligned_data.iter().fold(0, |acc, seq| acc + seq.len());
    //            seq_cnt += unaligned_data.len();
    //
    //            let mut covered_by_fswm: Vec<_> = unaligned_data
    //                .iter()
    //                .map(|seq| vec![false; seq.len()])
    //                .collect();
    //            let pos_in_aligned_data: Vec<Vec<usize>> = alignment
    //                .sequences
    //                .iter()
    //                .map(|seq| {
    //                    seq.data
    //                        .chars()
    //                        .enumerate()
    //                        .filter(|(pos, ch)| *ch != '-' && ch.is_ascii() && !ch.is_whitespace())
    //                        .map(|(pos, _)| pos)
    //                        .collect()
    //                })
    //                .collect();
    //            let patterns = ["0010011001"]
    //                .iter()
    //                .map(|str_pattern| str_pattern.parse().unwrap());
    //            let matches: Vec<_> = patterns
    //                .map(|pattern| {
    //                    let matches = find_word_match_buckets(&pattern, &unaligned_data);
    //                    (pattern, matches)
    //                })
    //                .collect();
    //
    //            let core_blocks: Vec<bool> = alignment
    //                .core_block
    //                .data
    //                .split(" ")
    //                .map(|s| s == "1")
    //                .collect();
    //
    //            // let match_count_threshold = 0.3;
    //            let match_count_threshold = 2;
    //            let mut core_match = 0;
    //            let mut partial_match = 0;
    //            let mut mismatch = 0;
    //
    //            for (pattern, matches) in matches {
    //                for (spaced_word, positions) in matches {
    //                    // remove matches in same sequence
    //                    // TODO use score!. Should be fine at the moment since it occurs
    //                    // rarely that a pattern occurs twice in the same sequence
    //                    let pair_matches = positions
    //                        .into_iter()
    //                        .map(|word_match| {
    //                            (
    //                                word_match.seq,
    //                                word_match.pos_in_seq,
    //                                pos_in_aligned_data[word_match.seq][word_match.pos_in_seq],
    //                                word_match.dont_care_word,
    //                            )
    //                        })
    //                        .tuple_combinations::<(_, _)>()
    //                        .filter(|(pos1, pos2)| pos1.0 != pos2.0)
    //                        .filter(|(pos1, pos2)| score_prot_pairwise(&pos1.3, &pos2.3) > 0);
    //
    //                    for (pos1, pos2) in pair_matches {
    //                        for pos in
    //                            covered_by_fswm[pos1.0][pos1.1..pos1.1 + pattern.len()].iter_mut()
    //                        {
    //                            *pos = true;
    //                        }
    //
    //                        for pos in
    //                            covered_by_fswm[pos2.0][pos2.1..pos2.1 + pattern.len()].iter_mut()
    //                        {
    //                            *pos = true;
    //                        }
    //
    //                        if (pos1.1..pos1.1 + pattern.len())
    //                            .zip(pos2.1..pos2.1 + pattern.len())
    //                            .all(|(p1, p2)| {
    //                                pos_in_aligned_data[pos1.0][p1] == pos_in_aligned_data[pos2.0][p2]
    //                            })
    //                        {
    //                            core_match += 1;
    //                        } else if (pos1.1..pos1.1 + pattern.len())
    //                            .zip(pos2.1..pos2.1 + pattern.len())
    //                            .any(|(p1, p2)| {
    //                                pos_in_aligned_data[pos1.0][p1] == pos_in_aligned_data[pos2.0][p2]
    //                            })
    //                        {
    //                            partial_match_sum += 1;
    //                        } else {
    //                            mismatch += 1;
    //                        }
    //                    }
    //
    //                    // let dont_care_alignment: Vec<_> = positions
    //                    //     .iter()
    //                    //     .map(|(_, _, dont_care_word)| dont_care_word)
    //                    //     .collect();
    //                    // if positions.len() < match_count_threshold || score_prot_msa(&dont_care_alignment) <= 0. {
    //                    //     continue;
    //                    // }
    //
    //                    // if positions.iter().all(|(_, pos, _)| *pos == positions[0].1)
    //                    //     && core_blocks[positions[0].1]
    //                    // {
    //                    //     core_match += 1;
    //                    // } else {
    //                    //     mismatch += 1;
    //                    // }
    //
    //                    // println!(
    //                    //     "{}: {:?}",
    //                    //     String::from_utf8(spaced_word).unwrap(),
    //                    //     positions
    //                    //         .iter()
    //                    //         .map(|(seq, pos)| (*seq, pos_in_aligned_data[*seq][*pos]))
    //                    //         .collect::<Vec<(usize, usize)>>()
    //                    // );
    //                }
    //            }
    //            let (hit, miss) = covered_by_fswm
    //                .iter()
    //                .fold((0, 0), |(mut hit, mut miss), seq| {
    //                    let (seq_hit, seq_miss) =
    //                        seq.iter().fold((0, 0), |(mut hit, mut miss), pos| {
    //                            if *pos {
    //                                hit += 1;
    //                            } else {
    //                                miss += 1;
    //                            }
    //                            (hit, miss)
    //                        });
    //                    hit += seq_hit;
    //                    miss += seq_miss;
    //                    (hit, miss)
    //                });
    //            println!("{}, {}", alignment.name, hit as f64 / (hit + miss) as f64);
    //            core_match_sum += core_match;
    //            partial_match_sum += partial_match;
    //            mismatch_sum += mismatch;
    //        }
    //        println!(
    //            "core_match: {}, mismatch: {}, partial_match: {}, ratio: {}",
    //            core_match_sum,
    //            mismatch_sum,
    //            partial_match_sum,
    //            core_match_sum as f64 / (mismatch_sum + core_match_sum) as f64
    //        );
    //        println!("Average len: {}", cum_len as f64 / seq_cnt as f64);
    //    }

    Ok(())
}
