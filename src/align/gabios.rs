use std::cmp::{max, min};
use std::collections::HashMap;
use std::convert::{TryFrom, TryInto};
use std::ops::{Index, IndexMut};

use fxhash::FxHashMap;
use itertools::Itertools;
use ndarray::{s, Array2, Array3, ArrayView2, Axis};
use num_integer::Integer;
use petgraph::graph::IndexType;
use petgraph::unionfind::UnionFind;

use crate::align::micro_alignment::{MicroAlignment, Site};
use std::cell::Cell;

#[derive(Clone, Debug, PartialEq)]
pub struct TransitiveClosure {
    pub succ: TransitiveFrontier,
    pub pred: TransitiveFrontier,

    succ_buf: TransitiveFrontier,
    pred_buf: TransitiveFrontier,
}

impl TransitiveClosure {
    pub fn new(max_seq_len: usize, seq_cnt: usize) -> Self {
        let pred = TransitiveFrontier::new(max_seq_len + 2, seq_cnt, FrontierKind::Pred);
        let succ = TransitiveFrontier::new(max_seq_len + 2, seq_cnt, FrontierKind::Succ);

        let pred_buf = pred.clone();
        let succ_buf = succ.clone();

        Self {
            succ,
            pred,
            pred_buf,
            succ_buf,
        }
    }

    pub fn add_diagonal(&mut self, diagonal: &mut MicroAlignment) -> bool {
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

    fn alignable_site_pairs<'a>(&'a self, diag: &'a MicroAlignment) -> Vec<(Site, Site)> {
        let consistent = diag.site_pair_iter().all(|(Site { seq, pos }, x)| {
            if seq == 2 && pos == 40 {
                let a = self.lower_bound(x, seq);
                let b = self.upper_bound(x, seq);
            }
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
        //        dbg!(&x, &y);
        let x_to_y_frontiers =
            self.pred.get((x, y.seq)) == y.pos && y.pos == self.succ.get((x, y.seq));
        let y_to_x_frontiers =
            self.pred.get((y, x.seq)) == x.pos && x.pos == self.succ.get((y, x.seq));
        let equal_frontiers = x_to_y_frontiers && y_to_x_frontiers;

        let ret = equal_frontiers;
        if ret != self.pred.eq_classes.equiv(x, y) {
            //            panic!("{:?}, {:?}, {}", x, y, ret);
        }
        ret
    }

    fn add_site_pair(&mut self, (site_a, site_b): (Site, Site)) {
        //        eprintln!("{:?}, {:?}", site_a, site_b);
        self.pred_buf.merge_eq_classes(site_a, site_b);
        self.succ_buf.merge_eq_classes(site_a, site_b);

        self.succ_buf.next_class.update(site_a, FrontierKind::Succ);
        self.succ_buf.next_class.update(site_b, FrontierKind::Succ);

        self.pred_buf.next_class.update(site_a, FrontierKind::Pred);
        self.pred_buf.next_class.update(site_b, FrontierKind::Pred);

        for eq_class_repr in self.succ_buf.iter_eq_classes() {
            for target_seq in 0..self.succ_buf.seq_cnt {
                let val = if self.shifted_le(eq_class_repr, site_a) {
                    min(
                        //                        succ_a, succ_b
                        self.succ.get((eq_class_repr, target_seq)),
                        self.succ.get((site_b, target_seq)),
                    )
                } else if self.shifted_le(eq_class_repr, site_b) {
                    min(
                        //                        succ_c,
                        //                        succ_d
                        self.succ.get((eq_class_repr, target_seq)),
                        self.succ.get((site_a, target_seq)),
                    )
                } else {
                    //                    succ_a
                    self.succ.get((eq_class_repr, target_seq))
                };
                self.succ_buf.set((eq_class_repr, target_seq), val);
            }
        }

        for eq_class_repr in self.pred_buf.iter_eq_classes() {
            for target_seq in 0..self.pred_buf.seq_cnt {
                let val = if self.shifted_ge(eq_class_repr, site_a) {
                    max(
                        self.pred.get((eq_class_repr, target_seq)),
                        self.pred.get((site_b, target_seq)),
                    )
                } else if self.shifted_ge(eq_class_repr, site_b) {
                    max(
                        self.pred.get((eq_class_repr, target_seq)),
                        self.pred.get((site_a, target_seq)),
                    )
                } else {
                    self.pred.get((eq_class_repr, target_seq))
                };
                self.pred_buf.set((eq_class_repr, target_seq), val);
            }
        }

        self.succ.clone_from(&self.succ_buf);
        self.pred.clone_from(&self.pred_buf);

        assert!(self.shifted_sites_are_aligned(site_a, site_b));
    }

    fn lower_bound(&self, origin_site: Site, target_seq: usize) -> usize {
        let pred = self.pred.get((origin_site, target_seq));
        let succ = self.succ.get((origin_site, target_seq));

        if pred == succ {
            pred
        } else {
            pred + 1
        }
    }

    fn upper_bound(&self, origin_site: Site, target_seq: usize) -> usize {
        let pred = self.pred.get((origin_site, target_seq));
        let succ = self.succ.get((origin_site, target_seq));

        if pred == succ {
            pred
        } else {
            succ - 1
        }
    }

    pub fn less(&self, left: Site, right: Site) -> bool {
        let left = shift_site(left, 1);
        let right = shift_site(right, 1);
        self.shifted_less(left, right)
    }

    fn shifted_less(&self, left: Site, right: Site) -> bool {
        let a = match self.pred.get((right, left.seq)) {
            val if val == self.pred.default_val => false,
            val => left.pos < val,
        };
        let b = match self.succ.get((left, right.seq)) {
            val if val == self.succ.default_val => false,
            val => val < right.pos,
        };

        a || b
    }

    fn shifted_le(&self, left: Site, right: Site) -> bool {
        match self.pred.get((right, left.seq)) {
            //            val if val == self.pred.default_val => false,
            val => left.pos <= val,
        }
    }

    pub fn greater(&self, left: Site, right: Site) -> bool {
        let left = shift_site(left, 1);
        let right = shift_site(right, 1);
        self.shifted_greater(left, right)
    }

    fn shifted_greater(&self, left: Site, right: Site) -> bool {
        let a = match self.succ.get((right, left.seq)) {
            val if val == self.succ.default_val => false,
            val => left.pos > val,
        };
        let b = match self.pred.get((left, right.seq)) {
            val if val == self.pred.default_val => false,
            val => val > right.pos,
        };
        a || b
    }

    fn shifted_ge(&self, left: Site, right: Site) -> bool {
        match self.succ.get((right, left.seq)) {
            //            val if val == self.succ.default_val => false,
            val => left.pos >= val,
        }
    }
}

fn shift_diagonal(diagonal: &mut MicroAlignment, shift: i64) {
    for start_site in diagonal.start_sites.iter_mut() {
        start_site.pos = (i64::try_from(start_site.pos).unwrap() + shift)
            .try_into()
            .unwrap();
    }
}

pub fn shift_site(a: Site, shift: i64) -> Site {
    Site {
        seq: a.seq,
        pos: (i64::try_from(a.pos).unwrap() + shift).try_into().unwrap(),
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct TransitiveFrontier {
    default_val: usize,
    kind: FrontierKind,
    seq_cnt: usize,
    next_class: NextClass,
    eq_classes: EqClasses,
    frontiers: FxHashMap<Site, Vec<Cell<usize>>>,
}

#[derive(Eq, PartialEq, Clone, Copy, Debug)]
enum FrontierKind {
    Pred,
    Succ,
}

impl FrontierKind {
    fn get_default_val(&self) -> usize {
        match self {
            FrontierKind::Pred => 0,
            FrontierKind::Succ => std::usize::MAX,
        }
    }
}

impl TransitiveFrontier {
    fn new(max_seq_len: usize, seq_cnt: usize, kind: FrontierKind) -> Self {
        let next_class = NextClass::new(max_seq_len, seq_cnt);
        let eq_classes = EqClasses::new(max_seq_len, seq_cnt);

        let frontiers = FxHashMap::default();
        Self {
            default_val: kind.get_default_val(),
            kind: kind,
            seq_cnt,
            next_class,
            eq_classes,
            frontiers,
        }
    }

    fn get(&self, index: (Site, usize)) -> usize {
        if index.0.seq == index.1 {
            return index.0.pos;
        }

        let repr = self.eq_classes.find(index.0);
        match self.frontiers.get(&repr) {
            Some(frontiers) => frontiers[index.1].get(),
            None => {
                let next_class_pos = self.next_class[index.0];
                let next_class = Site {
                    seq: index.0.seq,
                    pos: next_class_pos,
                };

                if next_class_pos == NextClass::DEFAULT_VAL {
                    return self.default_val;
                }

                let repr = self.eq_classes.find(next_class);
                let frontiers = self.frontiers.get(&repr).unwrap_or_else(|| {
                    panic!("Index: Next class repr {:?} not in frontiers", repr)
                });
                frontiers[index.1].get()
            }
        }
    }

    fn set(&self, index: (Site, usize), val: usize) {
        //        let old_val = self.frontiers.get(&index.0).unwrap()[index.1].get();
        //        if index.0.seq == index.1 && index.0.pos != val {
        //            panic!("{} {} {:?}", old_val, val, index)
        //        }
        self.frontiers
            .get(&index.0)
            .expect("set must be called with repr")[index.1]
            .set(val);
    }

    fn merge_eq_classes(&mut self, a: Site, b: Site) {
        let a_repr = self.eq_classes.find(a);
        let b_repr = self.eq_classes.find(b);

        if !self.eq_classes.union(a_repr, b_repr) {
            //            return ();
            panic!(
                            "Called TransitiveFrontier::align on sites [{:?}, {:?}], repr: {:?} that were already aligned. \
                        This is likely a bug.",
                            a, b, a_repr
                        )
        }
        let unified_repr = self.eq_classes.find(a_repr);
        let (remove_repr, update_repr) = if a_repr == unified_repr {
            (b_repr, a_repr)
        } else if b_repr == unified_repr {
            (a_repr, b_repr)
        } else {
            panic!(
                "Broke invariant that after union one of the \
                                former reprs must be the repr of the unified set"
            )
        };

        let (default_val, seq_cnt) = (self.default_val, self.seq_cnt);

        let removed_frontier = self.frontiers.remove(&remove_repr);
        let update_frontier = self.frontiers.get(&update_repr);

        let create_default_front = |site: Site| {
            let mut front = vec![Cell::new(default_val); seq_cnt];
            front[site.seq].set(site.pos);
            front
        };

        let do_front_update = |update_front: &Vec<_>, removed_front: Vec<_>| {
            let _removed_front = removed_front.clone();
            update_front.iter().enumerate().zip(removed_front).for_each(
                |((idx, l), r): ((_, &Cell<usize>), Cell<usize>)| {
                    let r = r.into_inner();

                    if self.kind == FrontierKind::Pred {
                        l.set(max(l.get(), r))
                    } else {
                        l.set(min(l.get(), r))
                    }
                },
            )
        };

        match (update_frontier, removed_frontier) {
            (Some(update_front), Some(removed_front)) => {
                do_front_update(update_front, removed_front);
            }
            (Some(update_front), None) => update_front[remove_repr.seq].set(remove_repr.pos),
            (None, Some(removed_front)) => {
                let mut front = create_default_front(update_repr);
                do_front_update(&front, removed_front);
                self.frontiers.insert(unified_repr, front);
            }
            (None, None) => {
                let mut front = create_default_front(update_repr);
                front[remove_repr.seq].set(remove_repr.pos);
                self.frontiers.insert(unified_repr, front);
            }
        }
    }

    pub fn iter_eq_classes(&self) -> impl Iterator<Item = Site> + '_ {
        self.frontiers.keys().copied()
    }
}

#[derive(Clone, Debug, PartialEq)]
struct NextClass {
    data: Array2<usize>,
}

impl NextClass {
    const DEFAULT_VAL: usize = 0;
    fn new(max_seq_len: usize, seq_cnt: usize) -> Self {
        //TODO this probably need to be initialized with something else than 0
        // also the shape might be off
        let data = Array2::zeros((seq_cnt, max_seq_len));
        Self { data }
    }

    fn update(&mut self, aligned_site: Site, frontier_kind: FrontierKind) {
        let orig_next_class = self.data[[aligned_site.seq, aligned_site.pos]];
        let mut slice = match frontier_kind {
            FrontierKind::Succ => self
                .data
                .slice_mut(s![aligned_site.seq, ..aligned_site.pos+1;-1]),
            FrontierKind::Pred => self
                .data
                .slice_mut(s![aligned_site.seq, aligned_site.pos..]),
        };
        slice
            .iter_mut()
            .take_while(|next_class| **next_class == orig_next_class)
            .for_each(|next_class| *next_class = aligned_site.pos);
    }
}

impl Index<Site> for NextClass {
    type Output = usize;

    fn index(&self, index: Site) -> &Self::Output {
        &self.data[[index.seq, index.pos]]
    }
}

impl IndexMut<Site> for NextClass {
    fn index_mut(&mut self, index: Site) -> &mut Self::Output {
        &mut self.data[[index.seq, index.pos]]
    }
}

#[derive(Clone, Debug)]
struct EqClasses {
    seq_cnt: usize,
    max_seq_len: usize,
    data: UnionFind<usize>,
}

impl EqClasses {
    fn new(max_seq_len: usize, seq_cnt: usize) -> Self {
        // TODO is this the correct number of elements
        let data = UnionFind::new(seq_cnt * max_seq_len);
        Self {
            seq_cnt,
            max_seq_len,
            data,
        }
    }
    fn find(&self, site: Site) -> Site {
        self.idx_to_site(self.data.find(self.site_to_idx(site)))
    }
    fn union(&mut self, a: Site, b: Site) -> bool {
        self.data.union(self.site_to_idx(a), self.site_to_idx(b))
    }

    fn equiv(&self, a: Site, b: Site) -> bool {
        self.data.equiv(self.site_to_idx(a), self.site_to_idx(b))
    }

    fn site_to_idx(&self, site: Site) -> usize {
        site.seq * self.max_seq_len + site.pos
    }

    fn idx_to_site(&self, idx: usize) -> Site {
        let (seq, pos) = idx.div_rem(&self.max_seq_len);
        Site { seq, pos }
    }
}

impl PartialEq for EqClasses {
    fn eq(&self, other: &Self) -> bool {
        self.seq_cnt == other.seq_cnt
            && self.max_seq_len == other.max_seq_len
            && self.data.clone().into_labeling() == other.data.clone().into_labeling()
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Not;

    use itertools::Itertools;
    use ndarray::{arr2, arr3};
    use smallvec::SmallVec;

    use super::*;

    macro_rules! s {
        ($seq:expr, $pos:expr) => {
            Site {
                seq: $seq,
                pos: $pos,
            }
        };
    }

    //    #[test]
    //    fn test_new_trans_frontier() {
    //        let actual: TransitiveFrontier<Option<usize>> = TransitiveFrontier::new(2, 3, None);
    //        let expected = arr3(&[
    //            [[Some(0), Some(1), Some(2), Some(3)], [None; 4], [None; 4]],
    //            [[None; 4], [Some(0), Some(1), Some(2), Some(3)], [None; 4]],
    //            [[None; 4], [None; 4], [Some(0), Some(1), Some(2), Some(3)]],
    //        ]);
    //        assert_eq!(actual.data, expected);
    //    }
    //
    //    #[test]
    //    fn iter_origin_sequences() {
    //        let mut actual: TransitiveFrontier<Option<usize>> = TransitiveFrontier::new(2, 3, None);
    //        let expected = vec![
    //            arr2(&[[Some(0), Some(1), Some(2), Some(3)], [None; 4], [None; 4]]),
    //            arr2(&[[None; 4], [Some(0), Some(1), Some(2), Some(3)], [None; 4]]),
    //            arr2(&[[None; 4], [None; 4], [Some(0), Some(1), Some(2), Some(3)]]),
    //        ];
    //
    //        assert_eq!(actual.iter_origin_sequences().collect_vec(), expected);
    //    }

    #[test]
    fn shift_diagonal_test() {
        let mut diag = MicroAlignment {
            start_sites: SmallVec::from_vec(vec![Site { seq: 2, pos: 3 }, Site { seq: 4, pos: 6 }]),
            k: 2,
        };

        let mut expected_diag = MicroAlignment {
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
        let mut closure = TransitiveClosure::new(2, 3);

        let a = Site { seq: 0, pos: 1 };
        let b = Site { seq: 1, pos: 1 };
        let c = Site { seq: 2, pos: 1 };

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

    #[test]
    fn greater_and_less() {
        let mut closure = TransitiveClosure::new(15, 10);
        let site_pairs = vec![(s!(0, 1), s!(1, 1))];

        for site_pair in site_pairs.into_iter() {
            closure.add_site_pair(site_pair);
        }
        assert!(closure.shifted_less(s!(0, 0), s!(1, 1)));
        assert!(closure.shifted_greater(s!(0, 0), s!(1, 1)).not());

        assert!(closure.shifted_greater(s!(1, 2), s!(0, 1)));
        assert!(closure.shifted_less(s!(1, 2), s!(0, 1)).not());
    }

    #[test]
    fn compare() {
        let mut closure = TransitiveClosure::new(15, 10);
        let site_pairs = vec![
            (s!(0, 1), s!(1, 1)),
            (s!(2, 1), s!(3, 1)),
            (s!(1, 2), s!(2, 2)),
        ];

        for site_pair in site_pairs.into_iter() {
            closure.add_site_pair(site_pair);
        }

        assert!(closure.shifted_greater(s!(1, 1), s!(2, 1)).not());
        assert!(closure.shifted_less(s!(1, 1), s!(2, 1)).not());
        assert!(closure.shifted_greater(s!(2, 2), s!(2, 1)));
        assert!(closure.shifted_less(s!(2, 1), s!(2, 2)));
    }

    #[test]
    fn ge_and_le() {
        let mut closure = TransitiveClosure::new(15, 10);
        let site_pairs = vec![(s!(0, 1), s!(1, 1))];

        for site_pair in site_pairs.into_iter() {
            closure.add_site_pair(site_pair);
        }
        assert!(closure.shifted_le(s!(0, 0), s!(1, 1)));
        assert!(closure.shifted_ge(s!(0, 0), s!(1, 1)).not());

        assert!(closure.shifted_ge(s!(1, 2), s!(0, 1)));
        assert!(closure.shifted_le(s!(1, 2), s!(0, 1)).not());
    }

    #[test]
    fn test_next_class_succ() {
        let mut next_class = NextClass::new(10, 1);
        next_class.update(s!(0, 3), FrontierKind::Succ);
        next_class.update(s!(0, 8), FrontierKind::Succ);

        assert_eq!(arr2(&[[3, 3, 3, 3, 8, 8, 8, 8, 8, 0]]), next_class.data);
    }

    #[test]
    fn test_next_class_pred() {
        let mut next_class = NextClass::new(10, 1);
        next_class.update(s!(0, 3), FrontierKind::Pred);
        next_class.update(s!(0, 8), FrontierKind::Pred);

        assert_eq!(arr2(&[[0, 0, 0, 3, 3, 3, 3, 3, 8, 8]]), next_class.data);
    }
}
