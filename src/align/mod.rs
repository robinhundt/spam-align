use crate::align::diagonal::{construct_diagonals_from_patterns, ScoredDiagonal};
use crate::align::gabios::TransitiveClosure;
use crate::score::score_prot_msa;
use crate::spaced_word::Pattern;
use crate::Sequences;
use itertools::Itertools;

pub mod diagonal;
pub mod eq_class;
pub mod gabios;

pub fn align(
    sequences: &Sequences,
    patterns: &[Pattern],
) -> (Vec<ScoredDiagonal>, TransitiveClosure) {
    let mut diagonals =
        construct_diagonals_from_patterns(patterns, sequences, score_prot_msa, Some(2));
    diagonals.sort_by_cached_key(|diag| -diag.score);

    println!("Found {} diagonals", diagonals.len());
    let max_seq_len = sequences
        .iter()
        .max_by_key(|seq| seq.len())
        .expect("There are no Sequences to align!")
        .len();

    let mut transitive_closure = TransitiveClosure::new(max_seq_len, sequences.len());

    let added_diagonals = diagonals
        .into_iter()
        .enumerate()
        .filter_map(|(id, mut scored_diag)| {
            if id % 1000 == 0 {
                println!("Adding no {}", id);
            }
            if transitive_closure.add_diagonal(&mut scored_diag.diag) {
                Some(scored_diag)
            } else {
                None
            }
        })
        .collect_vec();

    for x in &added_diagonals {
        for y in &added_diagonals {}
    }

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
