use core::borrow::BorrowMut;
use std::ops::{Deref, Index, IndexMut, Not};

use ndarray::{s, Array3};

use crate::align::{Diagonal, Site};
use std::cell::RefCell;
use std::cmp::{max, min, Ordering};
use std::slice::{Iter, IterMut};
use std::vec::IntoIter;

#[derive(Debug, Clone)]
pub struct PartialAlignment {
    max_seq_len: usize,
    seq_cnt: usize,
    succ: TransFrontiers,
    pred: TransFrontiers,
    site_aligned_to_seq: Array3<bool>,
}

#[derive(Debug, Clone)]
struct TransFrontiers {
    data: Array3<Option<usize>>,
}

impl TransFrontiers {
    /// Creates a new transitivity frontier with dimensions $$seq_cnt x seq_cnt x max_seq_len$$
    /// It is initialized such that for every sequence, the consistency bounds with
    /// itself are equal to the corresponding position. For example: after initialization the
    /// consistency bound for the Site(seq: 3, pos: 4) to seq: 3 will be 4, since a position
    /// in a sequence can only be aligned with itself in the same sequence
    ///
    /// The first dimension corresponds to the target sequence, the second to the origin sequence
    /// and the last dimension to the position in the origin sequence
    fn new(max_seq_len: usize, seq_cnt: usize) -> Self {
        let mut data = Array3::from_elem((seq_cnt, seq_cnt, max_seq_len), None);
        for i in 0..seq_cnt {
            data.slice_mut(s!(i, i, ..))
                .iter_mut()
                .enumerate()
                .for_each(|(idx, data)| *data = Some(idx))
        }

        Self { data }
    }

    /// Given two sites and a target sequence, determine the minimal transitivity bound position
    /// in that sequence
    ///
    /// Returns None if neither of the Sites is aligned with the target sequence
    fn min(&self, sites: [Site; 2], seq: usize) -> Option<usize> {
        match (self[(sites[0], seq)], self[(sites[1], seq)]) {
            (Some(pos1), Some(pos2)) => Some(min(pos1, pos2)),
            (None, Some(pos)) => Some(pos),
            (Some(pos), None) => Some(pos),
            _ => None,
        }
    }
    /// Given two sites and a target sequence, determine the maximal transitivity bound position
    /// in that sequence
    ///
    /// Returns None if neither of the Sites is aligned with the target sequence
    fn max(&self, sites: [Site; 2], seq: usize) -> Option<usize> {
        match (self[(sites[0], seq)], self[(sites[1], seq)]) {
            (Some(pos1), Some(pos2)) => Some(max(pos1, pos2)),
            (None, Some(pos)) => Some(pos),
            (Some(pos), None) => Some(pos),
            _ => None,
        }
    }
}

impl Index<(Site, usize)> for TransFrontiers {
    type Output = Option<usize>;

    fn index(&self, index: (Site, usize)) -> &Self::Output {
        &self.data[[index.1, index.0.seq, index.0.pos]]
    }
}

impl IndexMut<(Site, usize)> for TransFrontiers {
    fn index_mut(&mut self, index: (Site, usize)) -> &mut Self::Output {
        &mut self.data[[index.1, index.0.seq, index.0.pos]]
    }
}

impl PartialAlignment {
    /// Constructs a new empty partial alignment with initialized transitivity bounds
    pub fn new(max_seq_len: usize, seq_cnt: usize) -> Self {
        let succ = TransFrontiers::new(max_seq_len, seq_cnt);
        let pred = TransFrontiers::new(max_seq_len, seq_cnt);
        let site_aligned_to_seq = Array3::from_elem((seq_cnt, seq_cnt, max_seq_len), false);

        Self {
            succ,
            pred,
            max_seq_len,
            seq_cnt,
            site_aligned_to_seq,
        }
    }

    pub fn add_diagonal(&mut self, diag: &Diagonal) -> bool {
        if !self.is_consistent(diag) {
            return false;
        }

        for site_pair in diag.site_iter() {
            self.add_site_pair(site_pair)
        }
        true
    }

    // TODO what happens if the sub sequences described by the diagonal are already
    // aligned? Are they aligned again?
    pub fn is_consistent(&self, diag: &Diagonal) -> bool {
        diag.site_iter().all(|(site_i, site_j)| {
            let i_consistent_with_j = self.lower_bound(&site_i, site_j.seq) <= site_j.pos
                && site_j.pos as i64
                    <= self
                        .upper_bound(&site_i, site_j.seq)
                        .map_or(-1, |bound| bound as i64);
            let j_consistent_with_i = self.lower_bound(&site_j, site_i.seq) <= site_i.pos
                && site_i.pos as i64
                    <= self
                        .upper_bound(&site_j, site_i.seq)
                        .map_or(-1, |bound| bound as i64);
            i_consistent_with_j && j_consistent_with_i
        })
    }

    fn add_site_pair(&mut self, sites: (Site, Site)) {
        if self.is_aligned_with_seq(&sites.0, sites.1.seq) {
            return;
        }
        self.site_aligned_to_seq[[sites.0.seq, sites.1.seq, sites.1.pos]] = true;
        self.site_aligned_to_seq[[sites.1.seq, sites.0.seq, sites.0.pos]] = true;

        let mut succ = self.succ.clone();
        let mut pred = self.pred.clone();

        //TODO this code can probably be parallelized
        for j in 0..self.seq_cnt {
            for i in 0..self.seq_cnt {
                for k in 0..self.max_seq_len {
                    // update the successor bounds
                    let site = Site { seq: i, pos: k };
                    succ[(site, j)] = if self.le(site, sites.0) {
                        self.succ.min([site, sites.1], j)
                    } else if self.le(site, sites.1) {
                        self.succ.min([site, sites.0], j)
                    } else {
                        self.succ[(site, j)]
                    };

                    //update the predecessor bounds

                    pred[(site, j)] = if self.ge(site, sites.0) {
                        self.pred.max([site, sites.1], j)
                    } else if self.ge(site, sites.1) {
                        self.pred.max([site, sites.0], j)
                    } else {
                        self.pred[(site, j)]
                    };
                }
            }
        }

        self.succ = succ;
        self.pred = pred;
    }

    fn upper_bound(&self, site: &Site, seq: usize) -> Option<usize> {
        if self.site_aligned_to_seq[[seq, site.seq, site.pos]] {
            let succ = self.succ[(*site, seq)].unwrap_or_else(|| {
                panic!(
                    "succ[({:?}, {})] is None! Invalid call of upper_cound",
                    site, seq
                )
            });
            return Some(succ);
        }

        match self.succ[(*site, seq)] {
            None => Some(self.max_seq_len - 1),
            Some(0) => None,
            Some(succ) => Some(succ - 1),
        }
    }

    fn lower_bound(&self, site: &Site, seq: usize) -> usize {
        if self.is_aligned_with_seq(site, seq) {
            return self.pred[(*site, seq)].expect("pred[(Site, i)] is None!");
        }
        match self.pred[(*site, seq)] {
            None => 0,
            Some(pred) => pred + 1,
        }
    }

    pub fn le(&self, left: Site, right: Site) -> bool {
        match self.succ[(left, right.seq)] {
            Some(pos) => pos <= right.pos,
            None => false,
        }
    }

    pub fn ge(&self, left: Site, right: Site) -> bool {
        match self.pred[(left, right.seq)] {
            Some(pos) => pos >= right.pos,
            None => false,
        }
    }

    pub fn sites_are_aligned(&self, fst: Site, snd: Site) -> bool {
        match [self.succ[(fst, snd.seq)], self.pred[(fst, snd.seq)]] {
            [Some(pos_succ), Some(pos_pred)] => pos_succ == snd.pos && pos_pred == snd.pos,
            _ => false,
        }
    }

    fn is_aligned_with_seq(&self, site: &Site, seq: usize) -> bool {
        self.site_aligned_to_seq[[seq, site.seq, site.pos]]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::arr3;

    #[test]
    fn test_new_trans_frontier() {
        let actual = TransFrontiers::new(2, 3);
        let expected = arr3(&[
            [[Some(0), Some(1)], [None; 2], [None; 2]],
            [[None; 2], [Some(0), Some(1)], [None; 2]],
            [[None; 2], [None; 2], [Some(0), Some(1)]],
        ]);
        assert_eq!(actual.data, expected);
    }
}
