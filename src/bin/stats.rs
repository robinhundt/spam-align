use std::error::Error;
use std::path::PathBuf;

use itertools::Itertools;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use structopt::StructOpt;

use fxhash::{FxHashMap, FxHashSet};
use num_integer::binomial;
use spam_align::align::micro_alignment::{
    compute_multi_dim_micro_alignment_information, construct_micro_alignments_from_patterns, Site,
};
use spam_align::data_loaders::balibase::FilterXmlFile;
use spam_align::data_loaders::{balibase, PositionAlignment};
use spam_align::score::score_prot_msa;
use spam_align::spaced_word::{read_patterns_from_file, Pattern};
use spam_align::Sequences;
use std::ops::Not;

type BoxResult<T> = Result<T, Box<dyn Error>>;

#[derive(StructOpt, Debug)]
struct Opt {
    #[structopt(name = "IN_FOLDER", parse(from_os_str))]
    input_folder: PathBuf,
    #[structopt(name = "PATTERN_SET", long = "pattern-set", parse(from_os_str))]
    pattern_set: PathBuf,
}

#[derive(Serialize, Deserialize)]
struct EvaluationResult {
    path: PathBuf,
    results: Vec<AlignmentResult>,
}

#[derive(Serialize, Deserialize)]
struct AlignmentResult {
    name: String,
    true_positive_site_pair_count: usize,
    false_positive_site_pair_count: usize,
    true_site_pair_count: usize,
    precision: f64,
    recall: f64,
}

fn compute_precision(
    path: &PathBuf,
    patterns: &[Pattern],
) -> BoxResult<FxHashMap<PositionAlignment, usize>> {
    let alignments_all_data = balibase::parse_xml_files_in_dir(path, FilterXmlFile::AllData)?;
    let precision_data = FxHashMap::default();
    for alignment in alignments_all_data {
        let sequences = Sequences::new(&alignment);
        let mut correct_site_pairs = FxHashSet::default();
        correct_site_pairs.reserve(sequences.total_len() * 10);
        let mut incorrect_site_pairs = FxHashSet::default();
        incorrect_site_pairs.reserve(sequences.total_len() * 10);
        let micro_alignments =
            construct_micro_alignments_from_patterns(patterns, &sequences, score_prot_msa, true);
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
                    PositionAlignment::Unknown => {}
                }
            }
        }
        let prec = correct_site_pairs.len() as f64
            / (correct_site_pairs.len() + incorrect_site_pairs.len()) as f64;
        println!("{}: {}", alignment.name, prec)
    }
    Ok(precision_data)
}

fn compute_recall(path: &PathBuf, patterns: &[Pattern]) -> BoxResult<()> {
    let alignments_core_data = balibase::parse_xml_files_in_dir(path, FilterXmlFile::AllData)?;

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

//fn compute_multi_dim_stats(path: &PathBuf, patterns: &[Pattern]) -> BoxResult<()> {
//    let alignments_all_data = balibase::parse_xml_files_in_dir(path, FilterXmlFile::AllData)?;
//    for alignment in alignments_all_data {
//        let sequences = Sequences::new(&alignment);
//        let dim_data = compute_multi_dim_micro_alignment_information(&patterns, &sequences);
//        let mut dim_data = dim_data.iter().collect_vec();
//        dim_data.sort_by_cached_key(|el| el.0);
//        println!("{} : {:?}", alignment.name, dim_data);
//    }
//    Ok(())
//}

fn compute_results(path: PathBuf, pattern_set: &[Pattern]) -> EvaluationResult {
    todo!()
}

fn main() -> Result<(), Box<dyn Error>> {
    let opt = Opt::from_args();

    let data_path = opt.input_folder;

    let pattern_set = opt.pattern_set;

    //    let alignments_all_data =
    //        balibase::parse_xml_files_in_dir(data_path, FilterXmlFile::CoreBlockDataOnly)?;
    //    let mut match_count = 0;
    //    let mut total_count = 0;
    //    for alignment in alignments_all_data {
    //        for (seq1, seq2) in alignment.aligned_data.iter().tuple_combinations() {
    //            for (pos1, pos2) in seq1.data.iter().zip(seq2.data.iter()) {
    //                if *pos1 == b'-' || *pos2 == b'-' {
    //                    continue;
    //                }
    //                if pos1 != pos2 {
    //                    match_count += 1;
    //                }
    //                total_count += 1;
    //            }
    //        }
    //    }
    //    println!("{}", match_count as f64 / total_count as f64);

    //    if let Some(pattern_set_path) = pattern_set {
    //        let patterns = read_patterns_from_file(pattern_set_path)?;
    //        if multi_dim {
    //            compute_multi_dim_stats(&data_path, &patterns)?;
    //        }
    //        let _precision_data = compute_precision(&data_path, &patterns)?;
    //        //        println!("{:#?}", precision_data);
    //        println!();
    //        let _recall = compute_recall(&data_path, &patterns)?;
    //        //        println!("Recall: {}", recall)
    //    }

    //    buffered_file.write(format!("{};{}\n", spaced_word_match.score, consistent).as_bytes());

    Ok(())
}
