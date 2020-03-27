use anyhow::Result;
use fxhash::{FxBuildHasher, FxHashSet};
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use num_integer::binomial;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use spam_align::align::eq_class::EqClasses;
use spam_align::align::micro_alignment::Site;
use spam_align::align::{align, Align};
use spam_align::data_loaders::balibase::FilterXmlFile;
use spam_align::data_loaders::{
    balibase, format_as_fasta, write_as_fasta, Alignment, PositionAlignment, Sequence,
};
use spam_align::spaced_word::{read_patterns_from_file, Pattern};
use spam_align::Sequences;
use std::collections::HashSet;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::iter::FromIterator;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
enum Opt {
    Single {
        #[structopt(long, name = "FILE", parse(from_os_str))]
        in_file: PathBuf,
        #[structopt(long, name = "PATTERN_SET", parse(from_os_str))]
        pattern_set_path: PathBuf,
        #[structopt(long, name = "OUT", parse(from_os_str))]
        out_file: PathBuf,
    },
    Balibase {
        #[structopt(
            long,
            name = "BALIBASE PATH",
            parse(from_os_str),
            default_value = "data/bb3_release"
        )]
        balibase_path: PathBuf,
        #[structopt(
            long,
            name = "PATTERNS PATH",
            parse(from_os_str),
            default_value = "pattern_sets/data"
        )]
        patterns_path: PathBuf,
        #[structopt(short, long, name = "OUT", parse(from_os_str))]
        out_path: PathBuf,
    },
}

#[derive(Serialize, Deserialize)]
struct EvaluationResult {
    path: PathBuf,
    results: Vec<AlignmentResult>,
}

#[derive(Serialize, Deserialize, Debug)]
struct AlignmentResult {
    name: String,
    true_positive_site_pair_count: usize,
    false_positive_site_pair_count: usize,
    true_site_pair_count: usize,
    precision: f64,
    recall: f64,
    aligned_sequences_fasta: String,
    alignment_execution_time_ms: u128,
}

fn main() -> Result<()> {
    let opt: Opt = Opt::from_args();
    match opt {
        Opt::Single {
            in_file,
            pattern_set_path,
            out_file,
        } => {
            let alignment = balibase::parse_xml_file(in_file)?;
            let patterns = read_patterns_from_file(pattern_set_path)?;
            let results = compute_results_for_alignment(&alignment, &patterns, Align::ShowProgress);
            let out_file = File::create(out_file)?;
            serde_json::to_writer_pretty(out_file, &results)?;
        }
        Opt::Balibase {
            balibase_path,
            patterns_path,
            out_path,
        } => {
            compute_results_for_balibase(balibase_path, patterns_path, out_path)?;
        }
    }

    Ok(())
}

fn compute_results_for_balibase(
    balibase_path: PathBuf,
    patterns_path: PathBuf,
    out_path: PathBuf,
) -> Result<()> {
    let balibase_folders = ["RV11", "RV12", "RV20", "RV30", "RV50"];
    let balibase_folders = balibase_folders.iter().map(|folder| {
        let mut path = balibase_path.clone();
        path.push(folder);
        path
    });
    for folder in balibase_folders {
        let mut out_folder_path = out_path.clone();
        out_folder_path.push(folder.file_name().unwrap());

        fs::create_dir_all(&out_folder_path)?;
        for pattern_set_path in fs::read_dir("./pattern_sets/data")? {
            let pattern_set_path = pattern_set_path?;
            let pattern_set = read_patterns_from_file(pattern_set_path.path())?;
            let results = compute_results_for_folder(&folder, &pattern_set)?;
            let out_path = PathBuf::from_iter(
                [
                    out_folder_path.as_path(),
                    &pattern_set_path.file_name().as_ref(),
                ]
                .iter(),
            );
            let out_file = File::create(out_path)?;
            serde_json::to_writer_pretty(out_file, &results)?;
        }
    }
    Ok(())
}

fn compute_results_for_folder(
    path: impl AsRef<Path>,
    pattern_set: &[Pattern],
) -> Result<EvaluationResult> {
    let alignments = balibase::parse_xml_files_in_dir(&path, FilterXmlFile::AllData)?;
    let results = alignments
        .par_iter()
        .map(|alignment| compute_results_for_alignment(alignment, pattern_set, Align::HideProgress))
        .progress_count(alignments.len() as u64)
        .collect();
    Ok(EvaluationResult {
        path: PathBuf::from(path.as_ref()),
        results,
    })
}

fn compute_results_for_alignment(
    alignment: &Alignment,
    patterns: &[Pattern],
    progress: Align,
) -> AlignmentResult {
    let sequences = Sequences::new(&alignment);
    let mut correct_site_pairs: HashSet<(Site, Site), FxBuildHasher> = FxHashSet::default();
    correct_site_pairs.reserve(sequences.total_len() * 10);
    let mut incorrect_site_pairs: HashSet<(Site, Site), FxBuildHasher> = FxHashSet::default();
    incorrect_site_pairs.reserve(sequences.total_len() * 10);
    let now = Instant::now();
    let (scored_diagonals, closure) = align(&sequences, &patterns, progress);
    let eq_classes = EqClasses::new(&scored_diagonals, &closure);
    let mut orig_sequences = alignment.unaligned_data.clone();
    eq_classes.align_sequences(&mut orig_sequences);
    let alignment_execution_time = now.elapsed();

    for eq_class in eq_classes.iter() {
        for (&s1, &s2) in eq_class.iter().tuple_combinations() {
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

    let aligned_sequences_fasta = format_as_fasta(&orig_sequences);

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
        aligned_sequences_fasta,
        alignment_execution_time_ms: alignment_execution_time.as_millis(),
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
