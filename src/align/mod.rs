use std::ops::{Index, IndexMut};

use indicatif::{ProgressBar, ProgressIterator};
use itertools::Itertools;

use crate::align::gabios::Closure as TransitiveClosure;
use crate::align::micro_alignment::{
    construct_micro_alignments_from_patterns, ScoredMicroAlignment,
};
use crate::score::score_prot_msa;
use crate::Sequence;
use crate::spaced_word::Pattern;

pub mod gabios;
pub mod micro_alignment;

#[derive(Debug)]
pub enum AlignProgress {
    Show,
    Hide,
}

pub fn align(sequences: &mut [Sequence], patterns: &[Pattern], progress: AlignProgress) {
    let diagonals =
        construct_micro_alignments_from_patterns(patterns, sequences, score_prot_msa, false);
    let mut diagonals = diagonals.collect_vec();
    diagonals.sort_by_cached_key(|diag| -diag.score);

    let seq_lengths = sequences.iter().map(|seq| seq.len()).collect_vec();

    let mut transitive_closure = TransitiveClosure::new(&seq_lengths);
    let num_diagonals = diagonals.len();
    let mut diagonals: Box<dyn Iterator<Item = ScoredMicroAlignment>> =
        Box::new(diagonals.into_iter());

    if let AlignProgress::Show = progress {
        let progress_bar = ProgressBar::new(num_diagonals as u64);
        progress_bar.set_draw_delta(num_diagonals as u64 / 100);
        diagonals = Box::new(diagonals.progress_with(progress_bar));
    }

    diagonals.for_each(|scored_diag| {
        transitive_closure.try_add_micro_alignment(&scored_diag.micro_alignment);
    });

    transitive_closure.align(sequences);
}

#[derive(Default)]
pub struct Matrix<E> {
    rows: usize,
    cols: usize,
    data: Vec<E>,
}

impl<E> Matrix<E> {
    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn row(&self, row: usize) -> &[E] {
        let start = row * self.cols;
        let end = (row + 1) * self.cols;
        &self.data[start..end]
    }
}

impl Matrix<usize> {
    pub fn zeros(shape: [usize; 2]) -> Self {
        let data = vec![0; shape[0] * shape[1]];
        Self {
            rows: shape[0],
            cols: shape[1],
            data,
        }
    }

    pub fn resize_rows(&mut self, rows: usize) {
        self.rows = rows;
        self.data.resize(rows * self.cols, 0)
    }
}

impl<E> Index<[usize; 2]> for Matrix<E> {
    type Output = E;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        let idx = index[0] * self.cols + index[1];
        &self.data[idx]
    }
}

impl<E> IndexMut<[usize; 2]> for Matrix<E> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        let idx = index[0] * self.cols + index[1];
        &mut self.data[idx]
    }
}

impl From<bool> for AlignProgress {
    fn from(flag: bool) -> Self {
        if flag {
            AlignProgress::Show
        } else {
            AlignProgress::Hide
        }
    }
}

#[cfg(test)]
mod tests {

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
