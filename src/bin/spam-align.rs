use fxhash::FxHashSet;
use itertools::Itertools;
use num_integer::binomial;
use serde::{Deserialize, Serialize};
use spam_align::align::align;
use spam_align::align::eq_class::EqClasses;
use spam_align::data_loaders::{balibase, write_as_fasta, Alignment, PositionAlignment};
use spam_align::spaced_word::{read_patterns_from_file, Pattern};
use spam_align::Sequences;
use std::error::Error;
use std::path::PathBuf;
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
struct Opt {
    #[structopt(name = "FILE", parse(from_os_str))]
    in_file: PathBuf,
    #[structopt(name = "PATTERN_SET", parse(from_os_str))]
    pattern_set_path: PathBuf,
    #[structopt(name = "OUT", parse(from_os_str))]
    out_file: PathBuf,
}

#[derive(Serialize, Deserialize, Debug)]
struct AlignmentResult {
    name: String,
    true_positive_site_pair_count: usize,
    false_positive_site_pair_count: usize,
    true_site_pair_count: usize,
    precision: f64,
    recall: f64,
}

fn main() -> Result<(), Box<dyn Error>> {
    let opt: Opt = Opt::from_args();
    let alignment = balibase::parse_xml_file(opt.in_file)?;
    let patterns = read_patterns_from_file(opt.pattern_set_path)?;
    //    let sequences = Sequences::new(&alignment);
    //    let (scored_diagonals, closure) = align(&sequences, &patterns);
    //    eprintln!("Added {} diagonals", scored_diagonals.len());
    //    let mut orig_sequences = alignment.unaligned_data;
    //    let eq_classes = EqClasses::new(&scored_diagonals, &closure);
    //    eq_classes.align_sequences(&mut orig_sequences);
    //    write_as_fasta(opt.out_file, &orig_sequences)?;
    let results = compute_results_for_alignment(&alignment, &patterns);
    dbg!(results);
    Ok(())
}

fn compute_results_for_alignment(alignment: &Alignment, patterns: &[Pattern]) -> AlignmentResult {
    let sequences = Sequences::new(&alignment);
    let mut correct_site_pairs = FxHashSet::default();
    correct_site_pairs.reserve(sequences.total_len() * 10);
    let mut incorrect_site_pairs = FxHashSet::default();
    incorrect_site_pairs.reserve(sequences.total_len() * 10);
    let (scored_diagonals, closure) = align(&sequences, &patterns);
    let eq_classes = EqClasses::new(&scored_diagonals, &closure);

    for eq_class in eq_classes {
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
