use crate::align::gabios::TransitiveClosure;
use crate::align::micro_alignment::{
    construct_micro_alignments_from_patterns, ScoredMicroAlignment,
};
use crate::score::score_prot_msa;
use crate::spaced_word::Pattern;
use crate::Sequences;
use indicatif::{ProgressBar, ProgressIterator};
use itertools::Itertools;

pub mod eq_class;
pub mod gabios;
pub mod micro_alignment;

pub fn align(
    sequences: &Sequences,
    patterns: &[Pattern],
) -> (Vec<ScoredMicroAlignment>, TransitiveClosure) {
    let diagonals =
        construct_micro_alignments_from_patterns(patterns, sequences, score_prot_msa, false);
    let mut diagonals = diagonals.collect_vec();
    diagonals.sort_by_cached_key(|diag| -diag.score);

    println!("Found {} diagonals", diagonals.len());
    let max_seq_len = sequences
        .iter()
        .max_by_key(|seq| seq.len())
        .expect("There are no Sequences to align!")
        .len();

    let mut transitive_closure = TransitiveClosure::new(max_seq_len, sequences.len());
    let num_diagonals = diagonals.len();
    let progress_bar = ProgressBar::new(num_diagonals as u64);
    progress_bar.set_draw_delta(num_diagonals as u64 / 100);
    let added_diagonals = diagonals
        .into_iter()
        .progress_with(progress_bar)
        .filter_map(|mut scored_diag| {
            if transitive_closure.add_diagonal(&mut scored_diag.micro_alignment) {
                Some(scored_diag)
            } else {
                None
            }
        })
        .collect_vec();

    (added_diagonals, transitive_closure)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_loaders::balibase;
    use crate::spaced_word::read_patterns_from_file;

    #[test]
    fn test_align() {
        //        let alignment = balibase::parse_xml_file("data/bb3_release/RV30/BBS30004.xml").unwrap();
        //        let mut sequences = Sequences::new(&alignment);
        //        let patterns = read_patterns_from_file("./sample.pat").unwrap();
        //        let (diagonals, _) = align(&sequences, &patterns);
        //        dbg!(diagonals.len());
        //        panic!()
    }
}
