use std::error::Error;
use std::path::{Path, PathBuf};

use indicatif::ParallelProgressIterator;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use structopt::StructOpt;

use fxhash::{FxHashMap, FxHashSet};
use num_integer::binomial;
use spam_align::align::micro_alignment::construct_micro_alignments_from_patterns;
use spam_align::data_loaders::balibase::FilterXmlFile;
use spam_align::data_loaders::{balibase, Alignment, PositionAlignment};
use spam_align::score::score_prot_msa;
use spam_align::spaced_word::{read_patterns_from_file, Pattern};
use spam_align::Sequences;
use std::fs;
use std::fs::File;

type BoxResult<T> = Result<T, Box<dyn Error>>;

#[derive(StructOpt, Debug)]
enum Opt {
    Evaluate {
        #[structopt(name = "IN_FOLDER", parse(from_os_str))]
        data_path: PathBuf,
        #[structopt(name = "PATTERN_SET", long = "pattern-set", parse(from_os_str))]
        pattern_set_path: PathBuf,
        #[structopt(name = "OUT_FILE", long = "out", short = "o", parse(from_os_str))]
        out_path: PathBuf,
        #[structopt(long = "histogram", short = "h")]
        histogram_data: bool,
    },
    Balibase {
        #[structopt(
            name = "BB_RELEASE",
            help = "Path to the bb3_release folder",
            parse(from_os_str),
            default_value = "./data/bb3_release"
        )]
        data_path: PathBuf,
    },
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

#[derive(Serialize, Deserialize, Default)]
struct HistogramData {
    path: Option<PathBuf>,
    true_positive: FxHashMap<i32, usize>,
    false_positive: FxHashMap<i32, usize>,
    unknown: FxHashMap<i32, usize>,
}

fn main() -> BoxResult<()> {
    let opt = Opt::from_args();

    match opt {
        Opt::Evaluate {
            data_path,
            pattern_set_path,
            out_path,
            histogram_data,
        } => {
            let patterns = read_patterns_from_file(pattern_set_path)?;

            if histogram_data {
                let mut results = compute_histogram_data(&data_path, &patterns)?;
                results.path = Some(data_path);
                let out_file = File::create(out_path)?;
                serde_json::to_writer_pretty(out_file, &results)?;
            } else {
                let results = compute_results(data_path, &patterns)?;

                let out_file = File::create(out_path)?;
                serde_json::to_writer_pretty(out_file, &results)?;
            }
        }
        Opt::Balibase { data_path } => {
            let balibase_folders = ["RV11", "RV12", "RV20", "RV30", "RV50"];
            let balibase_folders = balibase_folders.iter().map(|folder| {
                let mut path = data_path.clone();
                path.push(folder);
                path
            });
            for folder in balibase_folders {
                let out_folder_path = format!(
                    "./stats-out-modified-cb/{}",
                    folder.file_name().unwrap().to_str().unwrap()
                );
                fs::create_dir_all(&out_folder_path)?;
                for pattern_set_path in fs::read_dir("./pattern_sets/data")? {
                    let pattern_set_path = pattern_set_path?;
                    let pattern_set = read_patterns_from_file(pattern_set_path.path())?;
                    let results = compute_results(&folder, &pattern_set)?;
                    let out_file = File::create(format!(
                        "{}/{}",
                        &out_folder_path,
                        pattern_set_path.file_name().to_str().unwrap()
                    ))?;
                    serde_json::to_writer_pretty(out_file, &results)?;
                }
            }
        }
    }

    Ok(())
}

fn compute_histogram_data(
    path: impl AsRef<Path>,
    pattern_set: &[Pattern],
) -> BoxResult<HistogramData> {
    let alignments = balibase::parse_xml_files_in_dir(&path, FilterXmlFile::AllData)?;
    let results = alignments
        .par_iter()
        .map(|alignment| compute_histogram_data_for_alignment(alignment, pattern_set))
        .progress_count(alignments.len() as u64)
        .reduce(
            || HistogramData::default(),
            |mut acc, data_for_alig| {
                for (score, count) in data_for_alig.unknown {
                    *acc.unknown.entry(score).or_insert(0) += count;
                }
                for (score, count) in data_for_alig.true_positive {
                    *acc.true_positive.entry(score).or_insert(0) += count;
                }
                for (score, count) in data_for_alig.false_positive {
                    *acc.false_positive.entry(score).or_insert(0) += count;
                }

                acc
            },
        );
    Ok(results)
}

fn compute_histogram_data_for_alignment(
    alignment: &Alignment,
    patterns: &[Pattern],
) -> HistogramData {
    let mut data = HistogramData::default();
    let sequences = Sequences::new(&alignment);
    let micro_alignments =
        construct_micro_alignments_from_patterns(patterns, &sequences, score_prot_msa, true);
    for micro_alignment in micro_alignments {
        for (s1, s2) in micro_alignment.micro_alignment.site_pair_iter() {
            let consistent = alignment.pos_aligned(s1, s2);
            match consistent {
                PositionAlignment::Correct => {
                    *data.true_positive.entry(micro_alignment.score).or_insert(0) += 1
                }
                PositionAlignment::Incorrect => {
                    *data
                        .false_positive
                        .entry(micro_alignment.score)
                        .or_insert(0) += 1
                }
                PositionAlignment::Unknown => {
                    *data.unknown.entry(micro_alignment.score).or_insert(0) += 1
                }
            }
        }
    }
    data
}

fn compute_results(path: impl AsRef<Path>, pattern_set: &[Pattern]) -> BoxResult<EvaluationResult> {
    let alignments = balibase::parse_xml_files_in_dir(&path, FilterXmlFile::AllData)?;
    let results = alignments
        .par_iter()
        .map(|alignment| compute_results_for_alignment(alignment, pattern_set))
        .progress_count(alignments.len() as u64)
        .collect();
    Ok(EvaluationResult {
        path: PathBuf::from(path.as_ref()),
        results,
    })
}

fn compute_results_for_alignment(alignment: &Alignment, patterns: &[Pattern]) -> AlignmentResult {
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

    let true_site_pair_count = true_site_pair_count(alignment);
    let tp = correct_site_pairs.len();
    let fp = incorrect_site_pairs.len();
    AlignmentResult {
        name: alignment.name.clone(),
        true_positive_site_pair_count: tp,
        false_positive_site_pair_count: fp,
        true_site_pair_count,
        precision: tp as f64 / (tp + fp) as f64,
        recall: tp as f64 / true_site_pair_count as f64,
    }
}

fn true_site_pair_count(alignment: &Alignment) -> usize {
    let core_block_data = alignment.core_block_data();
    let aligned_seq_len = core_block_data[0].data.len();
    let aligned_pos_per_col = (0..aligned_seq_len).map(|pos| {
        core_block_data.iter().fold(
            0,
            |acc, seq| {
                if seq.data[pos] == b'-' {
                    acc
                } else {
                    acc + 1
                }
            },
        )
    });
    let max_site_pair_count =
        aligned_pos_per_col.fold(0, |acc, aligned_count| acc + binomial(aligned_count, 2) * 2);

    max_site_pair_count
}
