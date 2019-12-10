use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use bio::io::fasta::Reader;
use itertools::Itertools;
use rayon::prelude::*;
use structopt::StructOpt;

use fxhash::{FxHashMap, FxHashSet};
use num_integer::binomial;
use spam_align::align::micro_alignment::{
    compute_multi_dim_micro_alignment_information, construct_micro_alignments_from_patterns, Site,
};
use spam_align::data_loaders::balibase::{BBAlignment, FilterXmlFile};
use spam_align::data_loaders::{balibase, MicroAlignmentCheck, PositionAlignment};
use spam_align::score::{score_prot_msa, score_prot_pairwise};
use spam_align::spaced_word::{
    find_word_match_buckets, find_word_matches, read_patterns_from_file, Pattern,
};
use spam_align::Sequences;
use std::iter::FromIterator;
use std::ops::Not;

#[derive(StructOpt, Debug)]
struct Opt {
    #[structopt(name = "IN_FOLDER", parse(from_os_str))]
    input_folder: PathBuf,
    #[structopt(name = "PATTERN", long = "pattern")]
    pattern: Option<Pattern>,
    #[structopt(name = "PATTERN_SET", long = "pattern-set", parse(from_os_str))]
    pattern_set: Option<PathBuf>,
    #[structopt(long = "multi-dim")]
    multi_dim: bool,
}

type BoxResult<T> = Result<T, Box<dyn Error>>;

fn compute_precision(
    path: &PathBuf,
    patterns: &[Pattern],
) -> BoxResult<FxHashMap<PositionAlignment, usize>> {
    let alignments_all_data = balibase::parse_xml_files_in_dir(path, FilterXmlFile::AllData)?;
    let mut precision_data = FxHashMap::default();
    for alignment in alignments_all_data {
        let sequences = Sequences::new(&alignment);
        let mut correct_site_pairs = FxHashSet::default();
        correct_site_pairs.reserve(sequences.total_len() * 10);
        let mut incorrect_site_pairs = FxHashSet::default();
        incorrect_site_pairs.reserve(sequences.total_len() * 10);
        let mut unknown_site_pairs = FxHashSet::default();
        unknown_site_pairs.reserve(sequences.total_len() * 10);
        let mut micro_alignments =
            construct_micro_alignments_from_patterns(patterns, &sequences, score_prot_msa, true)
                .collect_vec();
        micro_alignments.sort_by_cached_key(|ma| -ma.score);
        println!("Counting precision");
        for micro_alignment in micro_alignments {
            for (s1, s2) in micro_alignment.micro_alignment.site_pair_iter() {
                let consistent = alignment.pos_aligned(s1, s2);
                match consistent {
                    PositionAlignment::Correct => {
                        correct_site_pairs.insert((s1, s2));
                        correct_site_pairs.insert((s2, s1));
                    }
                    PositionAlignment::Incorrect => {
                        incorrect_site_pairs.insert((s1, s2));
                        incorrect_site_pairs.insert((s2, s1));
                    }
                    PositionAlignment::Unknown => {
                        unknown_site_pairs.insert((s1, s2));
                        unknown_site_pairs.insert((s2, s1));
                    }
                }
            }
        }
        *precision_data
            .entry(PositionAlignment::Correct)
            .or_insert(0) += correct_site_pairs.len();
        *precision_data
            .entry(PositionAlignment::Incorrect)
            .or_insert(0) += incorrect_site_pairs.len();
        *precision_data
            .entry(PositionAlignment::Unknown)
            .or_insert(0) += unknown_site_pairs.len();
    }
    Ok(precision_data)
}

fn compute_recall(path: &PathBuf, patterns: &[Pattern]) -> BoxResult<()> {
    let alignments_core_data =
        balibase::parse_xml_files_in_dir(path, FilterXmlFile::CoreBlockDataOnly)?;

    alignments_core_data.into_par_iter().for_each(|alignment| {
        let mut matched_site_pairs: FxHashSet<(Site, Site)> = FxHashSet::default();
        let sequences = Sequences::new(&alignment);
        let aligned_seq_len = alignment.aligned_data[0].data.len();
        let aligned_pos_per_col =
            (0..aligned_seq_len).map(|pos| {
                alignment.aligned_data.iter().fold(0, |acc, seq| {
                    if seq.data[pos] == b'-' {
                        acc
                    } else {
                        acc + 1
                    }
                })
            });
        let max_site_pair_count =
            aligned_pos_per_col.fold(0, |acc, aligned_count| acc + binomial(aligned_count, 2) * 2);
        matched_site_pairs.reserve(max_site_pair_count);
        for micro_alignment in
            construct_micro_alignments_from_patterns(patterns, &sequences, score_prot_msa, true)
        {
            for (s1, s2) in micro_alignment.micro_alignment.site_pair_iter() {
                match alignment.pos_aligned(s1, s2) {
                    PositionAlignment::Correct => {
                        matched_site_pairs.insert((s1, s2));
                        matched_site_pairs.insert((s2, s1));
                    }
                    PositionAlignment::Incorrect | PositionAlignment::Unknown => {}
                }
            }
        }
        println!(
            "{}: {}",
            alignment.name,
            matched_site_pairs.len() as f64 / max_site_pair_count as f64
        )
    });
    Ok(())
}

fn recall_from_coverage(coverage: Vec<Vec<bool>>) -> (u64, u64) {
    let mut match_and_total_per_column: Vec<(u64, u64)> = vec![];
    for seq in coverage {
        for (pos, aligned) in seq.into_iter().enumerate() {
            match match_and_total_per_column.get_mut(pos) {
                Some(counts) => {
                    if aligned {
                        counts.0 += 1
                    }
                    counts.1 += 1
                }
                None => match_and_total_per_column.push({
                    if aligned {
                        (1, 1)
                    } else {
                        (0, 1)
                    }
                }),
            }
        }
    }
    match_and_total_per_column.into_iter().fold(
        (0, 0),
        |(match_acc, total_acc), (match_col, total_col)| {
            (match_acc + match_col, total_acc + total_col)
        },
    )
}

fn compute_multi_dim_stats(path: &PathBuf, patterns: &[Pattern]) -> BoxResult<()> {
    let alignments_all_data = balibase::parse_xml_files_in_dir(path, FilterXmlFile::AllData)?;
    for alignment in alignments_all_data {
        let sequences = Sequences::new(&alignment);
        let dim_data =
            compute_multi_dim_micro_alignment_information(&patterns, &sequences, score_prot_msa);
        let mut dim_data = dim_data.iter().collect_vec();
        dim_data.sort_by_cached_key(|el| el.0);
        println!("{} : {:?}", alignment.name, dim_data);
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let opt = Opt::from_args();

    let data_path = opt.input_folder;

    let pattern = opt.pattern;
    let pattern_set = opt.pattern_set;
    let multi_dim = opt.multi_dim;

    if pattern.is_some() && pattern_set.is_some()
        || pattern.is_some().not() && pattern_set.is_some().not()
    {
        panic!("You need to either set prec_pattern or rec_pattern_set (not both)");
    }

    if let Some(pattern_set_path) = pattern_set {
        let patterns = read_patterns_from_file(pattern_set_path)?;
        if multi_dim {
            compute_multi_dim_stats(&data_path, &patterns)?;
        }
        //        let precision_data = compute_precision(&data_path, &patterns)?;
        //        println!("{:#?}", precision_data);
        let recall = compute_recall(&data_path, &patterns)?;
        //        println!("Recall: {}", recall)
    }

    //    buffered_file.write(format!("{};{}\n", spaced_word_match.score, consistent).as_bytes());

    Ok(())
}
