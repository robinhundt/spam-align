#[macro_use]
extern crate anyhow;

use crate::align::micro_alignment::Site;
use crate::data_loaders::Alignment;
use itertools::{EitherOrBoth, Itertools};
use std::iter::FromIterator;
use std::ops::Deref;

pub mod align;
pub mod data_loaders;
pub mod score;
pub mod spaced_word;

#[derive(Debug)]
pub struct Sequences {
    seq_data: Vec<u8>,
    seq_start_indices: Vec<usize>,
}

#[derive(Debug)]
pub struct Sequence<'a> {
    data: &'a [u8],
}

impl Sequences {
    pub fn new(alignment: &Alignment) -> Self {
        let mut seq_indices = vec![];
        let seq_data = Vec::from_iter(
            alignment
                .unaligned_data
                .iter()
                .map(|seq| {
                    let s_len = seq.data.len();
                    let next_idx = match seq_indices.last() {
                        None => (0, s_len),
                        Some(&(_start, end)) => (end, end + s_len),
                    };
                    seq_indices.push(next_idx);
                    seq.data.iter().copied()
                })
                .flatten(),
        );
        let seq_start_indices = seq_indices.into_iter().map(|(start, _end)| start).collect();
        Self {
            seq_data,
            seq_start_indices,
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = Sequence> {
        self.seq_start_indices
            .iter()
            .zip_longest(self.seq_start_indices.iter().skip(1))
            .map(move |bounds| {
                let data = match bounds {
                    EitherOrBoth::Both(idx1, idx2) => &self.seq_data[*idx1..*idx2],
                    EitherOrBoth::Left(idx) => &self.seq_data[*idx..],
                    EitherOrBoth::Right(_) => unreachable!(
                        "The left iterator should never be exhausted before the right one"
                    ),
                };
                Sequence { data }
            })
    }

    pub fn len(&self) -> usize {
        self.seq_start_indices.len()
    }

    pub fn total_len(&self) -> usize {
        self.seq_data.len()
    }

    // TODO possible optimization: store the actual bounds of the seq in Sequences to reduce
    // branches in this code
    pub fn get(&self, index: usize) -> Sequence<'_> {
        let start_pos = self.seq_start_indices[index];
        let end_pos = self.seq_start_indices.get(index + 1);
        let data = match end_pos {
            Some(&end_pos) => &self.seq_data[start_pos..end_pos],
            None => &self.seq_data[start_pos..],
        };
        Sequence { data }
    }

    pub fn get_site_slice(&self, start_site: Site, slice_len: usize) -> &[u8] {
        let start_pos = self.seq_start_indices[start_site.seq];
        let end_pos = self.seq_start_indices.get(start_site.seq + 1);
        let data = match end_pos {
            Some(&end_pos) => &self.seq_data[start_pos..end_pos],
            None => &self.seq_data[start_pos..],
        };
        &data[start_site.pos..start_site.pos + slice_len]
    }
}

impl Deref for Sequence<'_> {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        self.data
    }
}

#[cfg(test)]
mod tests {
    use crate::Sequences;

    #[test]
    fn sequence_iter() {
        let data = vec![1, 2, 3, 1, 2, 3, 4, 4, 4, 4, 4];
        let indices = vec![0, 3, 6];
        let seqs = Sequences {
            seq_data: data,
            seq_start_indices: indices,
        };
        let iter_seqs: Vec<_> = seqs.iter().map(|seq| seq.data.clone()).collect();
        assert_eq!(
            vec![vec![1, 2, 3], vec![1, 2, 3], vec![4, 4, 4, 4, 4]],
            iter_seqs
        )
    }
}
