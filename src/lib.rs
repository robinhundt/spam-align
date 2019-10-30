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
    seq_start_indices: Vec<usize>,
}

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
                        Some((&start, &end)) => (end, end + s_len),
                        None => (0, s_len),
                    };
                    seq_indices.push(next_idx);
                    seq.data.iter().copied()
                })
                .flatten(),
        );
        let seq_start_indices = seq_indices.into_iter().map(|(start, end)| start).collect();
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
        let indices = vec![0, 3, 6];
        let seqs = Sequences {
            seq_data: data,
            seq_start_indices: indices,
        };
        let iter_seqs: Vec<_> = seqs.iter().map(|seq| seq.data.clone()).collect();
        assert_eq!(
            iter_seqs,
            vec![vec![1, 2, 3], vec![1, 2, 3], vec![4, 4, 4, 4, 4]]
        )
    }
}
