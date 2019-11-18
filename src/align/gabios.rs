use core::borrow::BorrowMut;
use std::ops::{Deref, Index, IndexMut, Not};

use ndarray::{s, Array3, ArrayView2, ArrayViewMut2, Axis};

use crate::align::diagonal::{Diagonal, Site};
use itertools::Itertools;
use rayon::iter::IntoParallelRefIterator;
use std::cell::RefCell;
use std::cmp::{max, min, Ordering};
use std::convert::{TryFrom, TryInto};
use std::slice::{Iter, IterMut};
use std::sync::{Arc, Mutex, RwLock};
use std::vec::IntoIter;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TransitiveClosure {
    succ: TransitiveFrontier<usize>,
    pred: TransitiveFrontier<usize>,

    site_aligned_to_seq: SiteAlignment,

    succ_buf: TransitiveFrontier<usize>,
    pred_buf: TransitiveFrontier<usize>,
}

/// TransitiveFrontiers store information
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TransitiveFrontier<E: Eq + Ord + Clone + Copy> {
    data: Array3<E>,
}

impl<E> TransitiveFrontier<E>
where
    E: Eq + Ord + Clone + Copy,
    E: From<usize>,
{
    /// Creates a new transitivity frontier with dimensions $$seq_cnt x seq_cnt x max_seq_len$$
    /// It is initialized such that for every sequence, the consistency bounds with
    /// itself are equal to the corresponding position. For example: after initialization the
    /// consistency bound for the Site(seq: 3, pos: 4) to seq: 3 will be 4, since a position
    /// in a sequence can only be aligned with itself in the same sequence
    ///
    /// The first dimension corresponds to the target sequence, the second to the origin sequence
    /// and the last dimension to the position in the origin sequence
    pub fn new(max_seq_len: usize, seq_cnt: usize, default: E) -> Self {
        let mut data = Array3::from_elem((seq_cnt, seq_cnt, max_seq_len + 2), default);
        for i in 0..seq_cnt {
            data.slice_mut(s!(i, i, ..))
                .iter_mut()
                .enumerate()
                .for_each(|(idx, data)| *data = idx.into())
        }

        Self { data }
    }

    pub fn iter_origin_sequences(&self) -> impl Iterator<Item = ArrayView2<'_, E>> {
        self.data.axis_iter(Axis(0))
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct SiteAlignment {
    data: Array3<bool>,
}

impl SiteAlignment {
    pub fn new(max_seq_len: usize, seq_cnt: usize) -> Self {
        let mut data = Array3::from_elem((seq_cnt, seq_cnt, max_seq_len + 2), false);
        for i in 0..seq_cnt {
            // every position in a sequence is aligned with itself
            data.slice_mut(s!(i, i, ..))
                .iter_mut()
                .for_each(|data| *data = true)
        }

        Self { data }
    }
}

impl TransitiveClosure {
    pub fn new(max_seq_len: usize, seq_cnt: usize) -> Self {
        let pred = TransitiveFrontier::new(max_seq_len, seq_cnt, 0);
        let succ = TransitiveFrontier::new(max_seq_len, seq_cnt, std::usize::MAX);
        let site_aligned_to_seq = SiteAlignment::new(max_seq_len, seq_cnt);
        let pred_buf = pred.clone();
        let succ_buf = succ.clone();

        Self {
            succ,
            pred,
            site_aligned_to_seq,
            pred_buf,
            succ_buf,
        }
    }

    pub fn add_diagonal(&mut self, diagonal: &mut Diagonal) -> bool {
        let mut added = false;
        // we shift the Diagonal because Sequences starting at pos 1 is easier for the algo
        shift_diagonal(diagonal, 1);
        // TODO Maybe we could parallelize at this place
        let alignable_site_pairs = self.alignable_site_pairs(diagonal);
        for site_pair in alignable_site_pairs {
            self.add_site_pair(site_pair);
            added = true;
        }
        shift_diagonal(diagonal, -1);
        added
    }

    pub fn sites_are_aligned(&self, x: Site, y: Site) -> bool {
        let x = shift_site(x, 1);
        let y = shift_site(y, 1);
        self.shifted_sites_are_aligned(x, y)
    }

    fn alignable_site_pairs<'a>(&'a self, diag: &'a Diagonal) -> Vec<(Site, Site)> {
        let consistent = diag.site_pair_iter().all(|(Site { seq, pos }, x)| {
            self.lower_bound(x, seq) <= pos && pos <= self.upper_bound(x, seq)
        });
        if !consistent {
            return vec![];
        }

        diag.site_pair_iter()
            .filter(|(a, b)| !self.shifted_sites_are_aligned(*a, *b))
            .collect_vec()
    }

    fn shifted_sites_are_aligned(&self, x: Site, y: Site) -> bool {
        let x_to_y_frontiers = self.pred[(x, y.seq)] == y.pos && y.pos == self.succ[(x, y.seq)];
        let y_to_x_frontiers = self.pred[(y, x.seq)] == x.pos && x.pos == self.succ[(y, x.seq)];
        let equal_frontiers = x_to_y_frontiers && y_to_x_frontiers;

        //        let x_to_y = self.site_aligned_to_seq[(x, y.seq)];
        //        let y_to_x = self.site_aligned_to_seq[(y, x.seq)];
        //        let sites_are_aligned = x_to_y && y_to_x;

        equal_frontiers //&& sites_are_aligned
    }

    fn add_site_pair(&mut self, (site_a, site_b): (Site, Site)) {
        for (origin_seq, mut target_pos_view) in self.succ.iter_origin_sequences().enumerate() {
            'next_target: for (target_seq, pos_view) in
                target_pos_view.axis_iter(Axis(0)).enumerate()
            {
                for pos in 1..pos_view.len() {
                    let origin_site = Site {
                        seq: origin_seq,
                        pos,
                    };
                    let mut no_further_changes = false;
                    self.succ_buf[(origin_site, target_seq)] =
                        if self.shifted_le(origin_site, site_a) {
                            min(
                                self.succ[(origin_site, target_seq)],
                                self.succ[(site_b, target_seq)],
                            )
                        } else if self.shifted_le(origin_site, site_b) {
                            min(
                                self.succ[(origin_site, target_seq)],
                                self.succ[(site_a, target_seq)],
                            )
                        } else {
                            no_further_changes = true;
                            self.succ[(origin_site, target_seq)]
                        };
                    if no_further_changes {
                        continue 'next_target;
                    }
                }
            }
        }

        for (origin_seq, mut target_pos_view) in self.pred.iter_origin_sequences().enumerate() {
            'next_target: for (target_seq, pos_view) in
                target_pos_view.axis_iter(Axis(0)).enumerate()
            {
                for pos in (1..pos_view.len()).rev() {
                    let origin_site = Site {
                        seq: origin_seq,
                        pos,
                    };
                    let mut no_further_changes = false;
                    self.pred_buf[(origin_site, target_seq)] =
                        if self.shifted_ge(origin_site, site_a) {
                            max(
                                self.pred[(origin_site, target_seq)],
                                self.pred[(site_b, target_seq)],
                            )
                        } else if self.shifted_ge(origin_site, site_b) {
                            max(
                                self.pred[(origin_site, target_seq)],
                                self.pred[(site_a, target_seq)],
                            )
                        } else {
                            no_further_changes = true;
                            self.pred[(origin_site, target_seq)]
                            // TODO when this happens, it's probably correct to skip to the
                            // next iteration of the target_pos_view for loop
                        };
                    if no_further_changes {
                        continue 'next_target;
                    }
                }
            }
        }

        self.succ.clone_from(&self.succ_buf);
        self.pred.clone_from(&self.pred_buf);

        self.site_aligned_to_seq[(site_a, site_b.seq)] = true;
        self.site_aligned_to_seq[(site_b, site_a.seq)] = true;
    }

    fn lower_bound(&self, origin_site: Site, target_seq: usize) -> usize {
        let pred = self.pred[(origin_site, target_seq)];
        let succ = self.succ[(origin_site, target_seq)];

        if pred == succ {
            pred
        } else {
            pred + 1
        }
    }

    fn upper_bound(&self, origin_site: Site, target_seq: usize) -> usize {
        let pred = self.pred[(origin_site, target_seq)];
        let succ = self.succ[(origin_site, target_seq)];

        if pred == succ {
            pred
        } else {
            succ - 1
        }
    }

    pub fn le(&self, left: Site, right: Site) -> bool {
        let left = shift_site(left, 1);
        let right = shift_site(right, 1);
        self.shifted_le(left, right)
    }

    fn shifted_le(&self, left: Site, right: Site) -> bool {
        self.succ[(left, right.seq)] <= right.pos
    }

    pub fn ge(&self, left: Site, right: Site) -> bool {
        let left = shift_site(left, 1);
        let right = shift_site(right, 1);
        self.shifted_ge(left, right)
    }

    fn shifted_ge(&self, left: Site, right: Site) -> bool {
        self.pred[(left, right.seq)] >= right.pos
    }
}

fn shift_diagonal(diagonal: &mut Diagonal, shift: i64) {
    for start_site in diagonal.start_sites.iter_mut() {
        start_site.pos = (i64::try_from(start_site.pos).unwrap() + shift)
            .try_into()
            .unwrap();
    }
}

fn shift_site(a: Site, shift: i64) -> Site {
    Site {
        seq: a.seq,
        pos: (i64::try_from(a.pos).unwrap() + shift).try_into().unwrap(),
    }
}

impl<E> Index<(Site, usize)> for TransitiveFrontier<E>
where
    E: Eq + Ord + Clone + Copy,
{
    type Output = E;
    /// Index the Transitive Frontier with a tuple of a site and a target sequence
    /// the underlying array will be indexed with like this:
    /// [origin sequence, target sequence, origin position]
    fn index(&self, index: (Site, usize)) -> &Self::Output {
        &self.data[[index.1, index.0.seq, index.0.pos]]
    }
}

impl<E> IndexMut<(Site, usize)> for TransitiveFrontier<E>
where
    E: Eq + Ord + Clone + Copy,
{
    fn index_mut(&mut self, index: (Site, usize)) -> &mut Self::Output {
        &mut self.data[[index.1, index.0.seq, index.0.pos]]
    }
}

impl Index<(Site, usize)> for SiteAlignment {
    type Output = bool;
    /// Index the Transitive Frontier with a tuple of a site and a target sequence
    /// the underlying array will be indexed with like this:
    /// [origin sequence, target sequence, origin position]
    fn index(&self, index: (Site, usize)) -> &Self::Output {
        &self.data[[index.1, index.0.seq, index.0.pos]]
    }
}

impl IndexMut<(Site, usize)> for SiteAlignment {
    fn index_mut(&mut self, index: (Site, usize)) -> &mut Self::Output {
        &mut self.data[[index.1, index.0.seq, index.0.pos]]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;
    use ndarray::{arr2, arr3};
    use smallvec::SmallVec;

    macro_rules! s {
        ($seq:expr, $pos:expr) => {
            Site {
                seq: $seq,
                pos: $pos,
            }
        };
    }

    #[test]
    fn test_new_trans_frontier() {
        let actual: TransitiveFrontier<Option<usize>> = TransitiveFrontier::new(2, 3, None);
        let expected = arr3(&[
            [[Some(0), Some(1), Some(2), Some(3)], [None; 4], [None; 4]],
            [[None; 4], [Some(0), Some(1), Some(2), Some(3)], [None; 4]],
            [[None; 4], [None; 4], [Some(0), Some(1), Some(2), Some(3)]],
        ]);
        assert_eq!(actual.data, expected);
    }

    #[test]
    fn iter_origin_sequences() {
        let mut actual: TransitiveFrontier<Option<usize>> = TransitiveFrontier::new(2, 3, None);
        let expected = vec![
            arr2(&[[Some(0), Some(1), Some(2), Some(3)], [None; 4], [None; 4]]),
            arr2(&[[None; 4], [Some(0), Some(1), Some(2), Some(3)], [None; 4]]),
            arr2(&[[None; 4], [None; 4], [Some(0), Some(1), Some(2), Some(3)]]),
        ];

        assert_eq!(actual.iter_origin_sequences().collect_vec(), expected);
    }

    #[test]
    fn shift_diagonal_test() {
        let mut diag = Diagonal {
            start_sites: SmallVec::from_vec(vec![Site { seq: 2, pos: 3 }, Site { seq: 4, pos: 6 }]),
            k: 2,
        };

        let mut expected_diag = Diagonal {
            start_sites: SmallVec::from_vec(vec![Site { seq: 2, pos: 4 }, Site { seq: 4, pos: 7 }]),
            k: 2,
        };

        shift_diagonal(&mut diag, 1);

        assert_eq!(diag, expected_diag)
    }

    #[test]
    fn add_site_pair() {
        let mut closure = TransitiveClosure::new(15, 10);

        let a = Site { seq: 0, pos: 1 };
        let b = Site { seq: 1, pos: 6 };
        let c = Site { seq: 2, pos: 10 };
        let d = Site { seq: 3, pos: 8 };

        closure.add_site_pair((a, b));
        assert!(closure.shifted_sites_are_aligned(a, b));
        closure.add_site_pair((d, c));
        assert!(closure.shifted_sites_are_aligned(d, c));
        assert!(closure.shifted_sites_are_aligned(b, c).not());
        closure.add_site_pair((b, d));
        assert!(closure.shifted_sites_are_aligned(b, c));
        assert!(closure.shifted_sites_are_aligned(a, c));
        assert!(closure.shifted_sites_are_aligned(a, d));
    }

    #[test]
    fn add_site_pair_reversed_sites() {
        let mut closures = [
            TransitiveClosure::new(15, 10),
            TransitiveClosure::new(15, 10),
        ];

        let a = Site { seq: 0, pos: 1 };
        let b = Site { seq: 1, pos: 6 };

        closures[0].add_site_pair((a, b));
        closures[1].add_site_pair((b, a));

        assert_eq!(closures[0], closures[1])
    }

    #[test]
    fn add_site_cycle() {
        let mut closure = TransitiveClosure::new(15, 10);

        let a = Site { seq: 0, pos: 1 };
        let b = Site { seq: 1, pos: 6 };
        let c = Site { seq: 2, pos: 8 };

        closure.add_site_pair((a, b));
        closure.add_site_pair((b, c));
        let tmp_closure = closure.clone();
        closure.add_site_pair((c, a));

        assert_eq!(closure, tmp_closure);
        assert!(closure.shifted_sites_are_aligned(a, b));
        assert!(closure.shifted_sites_are_aligned(b, c));
        assert!(closure.shifted_sites_are_aligned(c, a));
    }

    #[test]
    fn corner_case() {
        macro_rules! s {
            ($seq:expr, $pos:expr) => {
                Site {
                    seq: $seq,
                    pos: $pos,
                }
            };
        }
        let mut closure = TransitiveClosure::new(15, 10);
        let site_pairs = vec![
            (s!(0, 1), s!(1, 1)),
            (s!(0, 3), s!(1, 3)),
            (s!(1, 2), s!(3, 2)),
        ];

        for site_pair in site_pairs.into_iter() {
            closure.add_site_pair(site_pair);
        }
        assert!(closure.shifted_sites_are_aligned(s!(0, 2), s!(1, 2)).not());
        assert!(closure.shifted_sites_are_aligned(s!(0, 2), s!(3, 2)).not());
    }

    #[test]
    fn corner_case_2() {
        let mut closure = TransitiveClosure::new(15, 10);
        let site_pairs = vec![
            (s!(0, 1), s!(1, 1)),
            (s!(0, 1), s!(2, 1)),
            (s!(0, 1), s!(3, 1)),
        ];

        for site_pair in site_pairs.into_iter() {
            closure.add_site_pair(site_pair);
        }
        assert!(closure.shifted_sites_are_aligned(s!(1, 1), s!(3, 1)));
    }
}
