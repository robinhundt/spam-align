use std::ops::{Index, IndexMut};

use indicatif::{ProgressBar, ProgressIterator};
use itertools::Itertools;

use crate::align::gabios::Closure as TransitiveClosure;
use crate::align::micro_alignment::{
    construct_micro_alignments_from_patterns, ScoredMicroAlignment,
};
use crate::score::score_prot_msa;
use crate::spaced_word::Pattern;
use crate::{log_elapsed_time, Sequence};

pub mod gabios;
pub mod micro_alignment;

#[derive(Debug)]
pub enum AlignProgress {
    Show,
    Hide,
}

pub fn align(sequences: &mut [Sequence], patterns: &[Pattern], progress: AlignProgress) {
    let diagonals = log_elapsed_time("Constructing diagonals", || {
        let diagonals =
            construct_micro_alignments_from_patterns(patterns, sequences, score_prot_msa, false);
        let mut diagonals = diagonals.collect_vec();
        diagonals.sort_by_cached_key(|diag| -diag.score);
        diagonals
    });

    let seq_lengths = sequences.iter().map(|seq| seq.len()).collect_vec();

    let mut transitive_closure = TransitiveClosure::new(&seq_lengths);
    let num_diagonals = diagonals.len();
    let mut diagonals: Box<dyn Iterator<Item = ScoredMicroAlignment>> =
        Box::new(diagonals.into_iter());

    if let AlignProgress::Show = progress {
        let progress_bar = ProgressBar::new(num_diagonals as u64);
        progress_bar.set_draw_delta(num_diagonals as u64 / 100);
        progress_bar.println("Adding micro alignments...");
        diagonals = Box::new(diagonals.progress_with(progress_bar));
    }
    log_elapsed_time("Adding diagonals", || {
        diagonals.for_each(|scored_diag| {
            transitive_closure.try_add_micro_alignment(&scored_diag.micro_alignment);
        });
    });

    log_elapsed_time("Aligning sequences", || {
        transitive_closure.align(sequences);
    });
}

#[derive(Default, Clone)]
pub struct Matrix<E: Clone> {
    rows: usize,
    cols: usize,
    elem: E,
    data: Vec<E>,
}

pub struct MatrixMutView<'a, E> {
    rows: usize,
    cols: usize,
    data: &'a mut [E],
}

impl<E: Clone> Matrix<E> {
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

    pub fn row_iter(&self) -> impl Iterator<Item = &[E]> + '_ {
        (0..self.rows).map(move |row| self.row(row))
    }

    pub fn row_mut(&mut self, row: usize) -> &mut [E] {
        let start = row * self.cols;
        let end = (row + 1) * self.cols;
        &mut self.data[start..end]
    }

    pub fn split_at_mut(&mut self, row_idx: usize) -> (MatrixMutView<'_, E>, &mut [E]) {
        let start = row_idx * self.cols;
        let (before, after) = self.data.split_at_mut(start);
        let row = &mut after[0..self.cols];
        let matrix_view = MatrixMutView {
            rows: row_idx - 1,
            cols: self.cols,
            data: before,
        };
        (matrix_view, row)
    }
}

impl Matrix<usize> {
    pub fn zeros(shape: [usize; 2]) -> Self {
        Self::from_elem(shape, 0)
    }

    pub fn from_elem(shape: [usize; 2], elem: usize) -> Self {
        let rows = shape[0];
        let cols = shape[1];
        let data = vec![elem; rows * cols];
        Self {
            rows,
            cols,
            elem,
            data,
        }
    }

    pub fn resize_rows(&mut self, rows: usize) {
        self.rows = rows;
        self.data.resize(rows * self.cols, 0)
    }
}

impl<'a, E> MatrixMutView<'a, E> {
    pub fn row(&'a self, row: usize) -> &[E] {
        let start = row * self.cols;
        let end = (row + 1) * self.cols;
        &self.data[start..end]
    }
}

impl<E: Clone> Index<[usize; 2]> for Matrix<E> {
    type Output = E;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        let idx = index[0] * self.cols + index[1];
        // &self.data[idx]
        // TODO this is really naughty thing and should not be done this way....
        // but it brings a nice performance boost. Change to better unsafe encapsulation
        // in future
        unsafe { self.data.get_unchecked(idx) }
    }
}

impl<E: Clone> IndexMut<[usize; 2]> for Matrix<E> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        let idx = index[0] * self.cols + index[1];
        // &mut self.data[idx]
        // TODO this is really naughty thing and should not be done this way....
        // but it brings a nice performance boost. Change to better unsafe encapsulation
        // in future
        unsafe { self.data.get_unchecked_mut(idx) }
    }
}

impl<E> Index<[usize; 2]> for MatrixMutView<'_, E> {
    type Output = E;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        let idx = index[0] * self.cols + index[1];
        // &self.data[idx]
        // TODO this is really naughty thing and should not be done this way....
        // but it brings a nice performance boost. Change to better unsafe encapsulation
        // in future
        unsafe { self.data.get_unchecked(idx) }
    }
}

impl<E> IndexMut<[usize; 2]> for MatrixMutView<'_, E> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        let idx = index[0] * self.cols + index[1];
        // &mut self.data[idx]
        // TODO this is really naughty thing and should not be done this way....
        // but it brings a nice performance boost. Change to better unsafe encapsulation
        // in future
        unsafe { self.data.get_unchecked_mut(idx) }
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
