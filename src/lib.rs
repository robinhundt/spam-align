use crate::data_loaders::Alignment;
use itertools::{EitherOrBoth, Itertools};
use std::iter::FromIterator;
use std::ops::Index;

pub mod align;
pub mod data_loaders;
pub mod score;
pub mod spaced_word;

pub struct Sequences {
    seq_data: Vec<u8>,
    seq_end_indices: Vec<usize>,
}

pub struct Sequence<'a> {
    data: &'a [u8],
}

impl Sequences {
    pub fn new(alignment: &Alignment) -> Self {
        let mut seq_end_indices = vec![];
        let seq_data = Vec::from_iter(
            alignment
                .unaligned_data
                .iter()
                .map(|seq| {
                    let end_idx = seq_end_indices.last().unwrap_or(&0) + seq.data.len();
                    seq_end_indices.push(end_idx);
                    seq.data.iter().copied()
                })
                .flatten(),
        );
        Self {
            seq_data,
            seq_end_indices,
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = Sequence> {
        [0].iter()
            .chain(self.seq_end_indices.iter())
            .zip_longest(self.seq_end_indices.iter())
            .map(move |bounds| {
                let data = match bounds {
                    EitherOrBoth::Both(idx1, idx2) => &self.seq_data[*idx1..*idx2],
                    EitherOrBoth::Left(idx) => &self.seq_data[*idx..],
                    EitherOrBoth::Right(idx) => unreachable!(
                        "The left iterator should never be exhausted before the right one"
                    ),
                };
                Sequence { data }
            })
    }

    //    pub fn get(&self, index: usize) -> Sequence<'_> {
    //        let start_pos = self.seq_start_indices[index];
    //        let end_pos = self.seq_start_indices.get(index + 1);
    //        let data = match end_pos {
    //            Some(&end_pos) => &self.seq_data[start_pos..end_pos],
    //            None => &self.seq_data[start_pos..],
    //        };
    //        Sequence { data }
    //    }
}

impl Index<usize> for Sequence<'_> {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        Index::index(self.data, index)
    }
}

#[cfg(test)]
mod tests {
    use crate::Sequences;

    #[test]
    fn sequence_iter() {
        let data = vec![1, 2, 3, 1, 2, 3, 4, 4, 4, 4, 4];
        let indices = vec![3, 6, data.len()];
        let seqs = Sequences {
            seq_data: data,
            seq_end_indices: indices,
        };
        let iter_seqs: Vec<_> = seqs.iter().map(|seq| seq.data.clone()).collect();
        assert_eq!(
            iter_seqs,
            vec![vec![1, 2, 3], vec![1, 2, 3], vec![4, 4, 4, 4, 4]]
        )
    }
}
