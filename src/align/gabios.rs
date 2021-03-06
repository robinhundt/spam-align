use std::cmp::{max, min, Ordering};
use std::ops::Not;

use itertools::repeat_n;
use smallvec::SmallVec;

use crate::align::micro_alignment::{MicroAlignment, Site};
use crate::align::Matrix;

pub struct Closure {
    /// Information about the positions in the aligned sequences.
    sequences: Sequences,
    /// This stores information about which eq class is aligned with which
    /// positions for each sequence
    ///
    /// Is used to tell if a eq class is aligned with a seq (not aligned if val == 0)
    alig_set: Matrix<usize>,

    /// How many eq classes are there in the current iteration.
    nbr_alig_sets: usize,

    /// these are prob used to save the actual frontiers,
    /// they are indexed with a eq class id and a sequence id
    pred_frontier: Matrix<usize>,
    /// these are prob used to save the actual frontiers,
    /// they are indexed with a eq class id and a sequence id
    succ_frontier: Matrix<usize>,

    /// This temporarily stores information about which pred/succ
    /// frontiers to update. First the necessary changes are determined
    /// and saved in this vec and then they are applied. This is needed
    /// because changing the frontiers while determining which ones need
    /// changing conflicts.
    pred_frontier_ops: Vec<FrontierOp>,
    succ_frontier_ops: Vec<FrontierOp>,
}

// Sequences is a struct of Matrices as opposed to a Vec of Structs of Vecs
// as in the original Gabios Lib implementation to improve cache friendliness
struct Sequences {
    lengths: Vec<usize>,
    /// index with seq pos, stores id of alig set this site belongs to
    alig_set_nbr: Matrix<usize>,
    /// next class information for pred frontier
    pred_alig_set_pos: Matrix<usize>,
    /// next class information for succ frontier
    succ_alig_set_pos: Matrix<usize>,
}

#[derive(Clone, Debug, Ord, PartialOrd, Eq, PartialEq)]
struct FrontierOp {
    seq: usize,
    eq_class: usize,
    new_frontier: usize,
}

enum Alignable {
    True,
    False,
    AlreadyAligned,
}

/// Wrapper for a site whose position is increased by 1
/// since the algorithm operates in sequences with `1` being the first position
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash, Default)]
struct ShiftedSite(Site);

impl Closure {
    pub fn new(seq_lengths: &[usize]) -> Self {
        let &max_length = seq_lengths.iter().max().expect("No sequences");
        let sequences = Sequences::new(seq_lengths);
        let nbr_seqs = sequences.len();

        let pred_frontier = Matrix::zeros([max_length + 2, nbr_seqs]);
        let succ_frontier = Matrix::from_elem([max_length + 2, nbr_seqs], usize::MAX);
        let alig_set = Matrix::zeros([max_length + 2, nbr_seqs]);

        let pred_frontier_ops = vec![];
        let succ_frontier_ops = vec![];
        Self {
            sequences,
            pred_frontier,
            succ_frontier,
            alig_set,
            nbr_alig_sets: 0,
            pred_frontier_ops,
            succ_frontier_ops,
        }
    }

    /// Try to add the passed micro alignment to the partial alignment
    /// returns false if the micro alignment is inconsistent with previously
    /// added ones or if it is already included in earlier ones
    pub fn try_add_micro_alignment(&mut self, micro_alignment: &MicroAlignment) -> bool {
        let mut added = false;
        let mut sites_to_add: SmallVec<[(ShiftedSite, ShiftedSite); 15]> = SmallVec::new();
        for (a, b) in micro_alignment.site_pair_iter() {
            let (a, b) = (a.into(), b.into());
            let alignable = self.alignable(a, b);
            match alignable {
                Alignable::True => sites_to_add.push((a, b)),
                Alignable::False => return false,
                Alignable::AlreadyAligned => continue,
            }
        }
        for (a, b) in sites_to_add {
            self.add_aligned_positions(a, b);
            added = true;
        }
        added
    }

    pub fn less(&self, left: Site, right: Site) -> bool {
        let (left, right) = (left.into(), right.into());
        self.path_exists(left, right)
    }

    /// Checks whether two sites are alignable
    fn alignable(&self, a: ShiftedSite, b: ShiftedSite) -> Alignable {
        // sites a and b are alignable if path(a, b) <=> path(b, a)
        if self.path_exists(a, b) {
            // if there is a path from a to b and b to a they are
            // already aligned and thus alignable
            if self.path_exists(b, a) {
                Alignable::AlreadyAligned
            } else {
                Alignable::False
            }
        } else if self.path_exists(b, a).not() {
            Alignable::True
        } else {
            Alignable::False
        }
    }

    /// Checks if there is a path in the alignment graph from the first
    /// to the second site
    fn path_exists(&self, from: ShiftedSite, to: ShiftedSite) -> bool {
        let (from, to) = (from.0, to.0);
        if from.seq == to.seq {
            return from.pos < to.pos;
        }
        let mut pred_alig_set_to = self.sequences.alig_set_nbr[[to.seq, to.pos]];
        if pred_alig_set_to == 0 {
            let k = self.sequences.pred_alig_set_pos[[to.seq, to.pos]];
            if k > 0 {
                pred_alig_set_to = self.sequences.alig_set_nbr[[to.seq, k]];
            }
        }
        if pred_alig_set_to == 0 {
            false
        } else {
            from.pos <= self.pred_frontier[[pred_alig_set_to, from.seq]]
        }
    }

    /// Uses the current closure to align the passed sequences
    pub fn align(&self, seqs: &mut [crate::Sequence]) {
        // Gather a vec of vec of sites, where each Vec<Site> represents
        // an equivalency class of aligned sites
        let mut classes = vec![vec![]; self.nbr_alig_sets];
        for (idx, seq_alig_set) in self.sequences.alig_set_nbr.row_iter().enumerate() {
            for (pos, alig_set) in seq_alig_set.iter().enumerate().skip(1) {
                if *alig_set == 0 {
                    continue;
                }
                let site: Site = ShiftedSite(Site { seq: idx, pos }).into();
                let alig_set = alig_set - 1;
                classes[alig_set].push(site);
            }
        }

        quicksort_by(&mut classes, &|a, b| {
            let (a, b) = (a[0], b[0]);
            if self.less(a, b) {
                Ordering::Less
            } else if self.less(b, a) {
                Ordering::Greater
            } else {
                Ordering::Equal
            }
        });

        // for seq in seqs.iter_mut() {
        //     for residue in seq.data.iter_mut() {
        //         residue.make_ascii_lowercase();
        //     }
        // }

        // This handles the actual gap insertion by iterating the eq classes
        // from smallest (meaning the left most) to highest,
        // inserting gaps such that all the positions of the class are in the
        // same column in the output. The remaining classes are shifted to
        // the right correspondingly
        classes.reverse();
        let mut shifted_by = vec![0; self.sequences.len()];
        while let Some(class) = classes.pop() {
            let max_pos = class
                .iter()
                .max_by_key(|site| site.pos)
                .expect("Unexpected empty eq class")
                .pos;
            shifted_by.iter_mut().for_each(|el| *el = 0);
            for site in &class {
                seqs[site.seq].data[site.pos].make_ascii_uppercase();
                let shift = max_pos - site.pos;
                if shift == 0 {
                    continue;
                }
                seqs[site.seq]
                    .data
                    .splice(site.pos..site.pos, repeat_n(b'-', shift));
                shifted_by[site.seq] = shift;
            }
            classes.iter_mut().for_each(|shifted_class| {
                shifted_class.iter_mut().for_each(|mut site| {
                    site.pos += shifted_by[site.seq];
                });
            });
        }

        // fill ends of seqs with gap '-'
        let max_len = seqs
            .iter()
            .max_by_key(|seq| seq.data.len())
            .expect("No sequences to align")
            .data
            .len();
        for seq in seqs {
            let missing_gaps = max_len - seq.data.len();
            seq.data.extend(repeat_n(b'-', missing_gaps));
        }
    }

    fn add_aligned_positions(&mut self, s1: ShiftedSite, s2: ShiftedSite) {
        let (s1, s2) = (s1.0, s2.0);

        // this condition is guaranteed by filtering the
        // micro alignment that is added in the `try_add_micro_alignment`
        // method
        // assert!(
        //     n1 != n2 || n1 == 0 || n2 == 0,
        //     "Sites must belong to different eq classes or at least one must not be aligned"
        // );

        // this closure looks up the alignment set n of a site and the alignment sites
        // that are closest and <= and >= than n. If site is aligned (n, n, n) is returned
        // with n being the set of site
        let lookup_alig_set = |site: Site| {
            let n = self.sequences.alig_set_nbr[[site.seq, site.pos]];
            if n != 0 {
                return (n, n, n);
            }
            let mut n_left = n;
            let mut n_right = n;
            let k = self.sequences.pred_alig_set_pos[[site.seq, site.pos]];
            if k > 0 {
                n_left = self.sequences.alig_set_nbr[[site.seq, k]];
            }
            let k = self.sequences.succ_alig_set_pos[[site.seq, site.pos]];
            if k > 0 {
                n_right = self.sequences.alig_set_nbr[[site.seq, k]];
            }
            (n, n_left, n_right)
        };

        let (n1, n_left1, n_right1) = lookup_alig_set(s1);
        let (n2, n_left2, n_right2) = lookup_alig_set(s2);

        if !(n1 != n2 || n1 == 0 || n2 == 0) {
            return;
        }
        // nn is the index of the alignment set created by aligning
        // s1 and s2, it might be moved at the end of the method if
        // existing sets were merged
        let nn = self.nbr_alig_sets + 1;
        let (mut sub_mat_pred, nn_pred) = self.pred_frontier.split_at_mut(nn);
        let (mut sub_mat_succ, nn_succ) = self.succ_frontier.split_at_mut(nn);
        let left1 = sub_mat_pred.row(n_left1);
        let left2 = sub_mat_pred.row(n_left2);

        let right1 = sub_mat_succ.row(n_right1);
        let right2 = sub_mat_succ.row(n_right2);

        // initialize the alignment set nn
        for x in 0..self.sequences.len() {
            self.alig_set[[nn, x]] = 0;
            if n1 > 0 && self.alig_set[[n1, x]] > 0 {
                self.alig_set[[nn, x]] = self.alig_set[[n1, x]];
            } else if n2 > 0 && self.alig_set[[n2, x]] > 0 {
                self.alig_set[[nn, x]] = self.alig_set[[n2, x]];
            }

            if self.alig_set[[nn, x]] == 0 {
                nn_pred[x] = max(left1[x], left2[x]);
                nn_succ[x] = min(right1[x], right2[x]);
            } else {
                let front = self.alig_set[[nn, x]];
                nn_pred[x] = front;
                nn_succ[x] = front;
            }
        }
        nn_pred[s1.seq] = s1.pos;
        nn_succ[s1.seq] = s1.pos;
        self.alig_set[[nn, s1.seq]] = s1.pos;

        nn_pred[s2.seq] = s2.pos;
        nn_succ[s2.seq] = s2.pos;
        self.alig_set[[nn, s2.seq]] = s2.pos;

        for x in 0..self.sequences.len() {
            let k = self.alig_set[[nn, x]];
            if k > 0 {
                self.sequences.alig_set_nbr[[x, k]] = nn;
            }
        }

        // change pred and succ frontiers for all eq classes != nn
        // if influenced by alignment of s1 and s2

        self.pred_frontier_ops.clear();

        for x in 0..self.sequences.len() {
            // SAFETY: do not remove this, because succ frontiers are initialized
            // with usize::MAX, removing this check will lead to buffer overflow
            if right1[x] == right2[x] {
                continue;
            }
            for y in 0..self.sequences.len() {
                let mut k = nn_succ[x];
                if k == self.alig_set[[nn, x]] {
                    // eq class nn is directly aligned with pos k in seq x
                    // so we take the the pos of the next aligned pos after k
                    k = self.sequences.succ_alig_set_pos[[x, k]];
                }

                // k is the nearest succ_frontier position in x from the perspective
                // of eq class nn, but not part of it
                // we then iterate over the eq classes present in x which have position
                // greater than the initial k and as long as their pred_frontier to **y**
                // (so from [x, k] to y) is less than the pred frontier off nn to y
                // we save the need to change the bound from n to y to the value of
                // nn_to_y in self.frontier_ops
                let nn_to_y = nn_pred[y];
                while k > 0 {
                    let n = self.sequences.alig_set_nbr[[x, k]];
                    if sub_mat_pred[[n, y]] < nn_to_y {
                        // frontier ops cannot be applied directly because
                        // that would influence above if statement
                        self.pred_frontier_ops.push(FrontierOp::new(n, y, nn_to_y));
                        k = self.sequences.succ_alig_set_pos[[x, k]];
                    } else {
                        break;
                    }
                }
            }
        }

        self.succ_frontier_ops.clear();

        for x in 0..self.sequences.len() {
            if left1[x] == left2[x] {
                continue;
            }
            for y in 0..self.sequences.len() {
                let mut k = nn_pred[x];
                if k > 0 && k == self.alig_set[[nn, x]] {
                    k = self.sequences.pred_alig_set_pos[[x, k]];
                }

                let nn_to_y = nn_succ[y];
                while k > 0 {
                    let n = self.sequences.alig_set_nbr[[x, k]];
                    if sub_mat_succ[[n, y]] > nn_to_y {
                        self.succ_frontier_ops.push(FrontierOp::new(n, y, nn_to_y));
                        k = self.sequences.pred_alig_set_pos[[x, k]];
                    } else {
                        break;
                    }
                }
            }
        }
        self.pred_frontier_ops.iter().for_each(
            |FrontierOp {
                 eq_class,
                 seq,
                 new_frontier,
             }| {
                sub_mat_pred[[*eq_class, *seq]] = *new_frontier;
            },
        );
        self.succ_frontier_ops.iter().for_each(
            |FrontierOp {
                 eq_class,
                 seq,
                 new_frontier,
             }| {
                sub_mat_succ[[*eq_class, *seq]] = *new_frontier;
            },
        );

        let mut update_sequences = |n: usize, site: Site| {
            if n != 0 {
                return;
            }

            let mut k = site.pos - 1;
            while k > 0 && self.sequences.alig_set_nbr[[site.seq, k]] == 0 {
                self.sequences.succ_alig_set_pos[[site.seq, k]] = site.pos;
                k -= 1;
            }
            if k > 0 {
                self.sequences.succ_alig_set_pos[[site.seq, k]] = site.pos;
            }

            let mut k = site.pos + 1;
            while k <= self.sequences.lengths[site.seq]
                && self.sequences.alig_set_nbr[[site.seq, k]] == 0
            {
                self.sequences.pred_alig_set_pos[[site.seq, k]] = site.pos;
                k += 1;
            }
            if k <= self.sequences.lengths[site.seq] {
                self.sequences.pred_alig_set_pos[[site.seq, k]] = site.pos
            }
        };

        update_sequences(n1, s1);
        update_sequences(n2, s2);

        let (n1, n2) = order(n1, n2);

        if n2 == 0 {
            // implies that n1 == 0
            self.nbr_alig_sets += 1;
            self.resize();
        } else if n1 == 0 {
            // implies that n2 != 0, which means that nn
            // should be move into n2
            self.move_alig_set(n2, nn)
        } else {
            // n1 != 0 && n2 != 0, which means
            // nn is merged from n1 and n2, thus
            // nn can be moved into n1 and the now last set
            // self.nbr_alig_sets can be moved into n2 which
            // is not needed anymore
            // self.nbr_alig_sets is reduced by 1
            self.move_alig_set(n1, nn);
            if n2 < self.nbr_alig_sets {
                self.move_alig_set(n2, self.nbr_alig_sets)
            }
            self.nbr_alig_sets -= 1;
        }
    }

    /// Move the alignment set n2 **into** n1
    fn move_alig_set(&mut self, n1: usize, n2: usize) {
        for x in 0..self.sequences.len() {
            let k = self.alig_set[[n2, x]];
            self.alig_set[[n1, x]] = k;
            if k > 0 {
                self.sequences.alig_set_nbr[[x, k]] = n1;
            }
            self.pred_frontier[[n1, x]] = self.pred_frontier[[n2, x]];
            self.succ_frontier[[n1, x]] = self.succ_frontier[[n2, x]];
        }
    }

    fn resize(&mut self) {
        let new_rows = self.nbr_alig_sets + 2;
        self.pred_frontier.resize_rows(new_rows);
        self.succ_frontier.resize_rows(new_rows);
        self.alig_set.resize_rows(new_rows);
    }
}

impl Sequences {
    fn new(seq_lengths: &[usize]) -> Self {
        let max_len = seq_lengths.iter().max().expect("No sequences") + 2;
        let mat = Matrix::zeros([seq_lengths.len(), max_len]);
        Sequences {
            lengths: seq_lengths.to_vec(),
            alig_set_nbr: mat.clone(),
            pred_alig_set_pos: mat.clone(),
            succ_alig_set_pos: mat,
        }
    }

    fn len(&self) -> usize {
        self.lengths.len()
    }
}

impl FrontierOp {
    fn new(eq_class: usize, seq: usize, new_frontier: usize) -> Self {
        Self {
            eq_class,
            seq,
            new_frontier,
        }
    }
}

impl From<Site> for ShiftedSite {
    fn from(mut site: Site) -> Self {
        site.pos += 1;
        ShiftedSite(site)
    }
}

impl From<ShiftedSite> for Site {
    fn from(mut site: ShiftedSite) -> Self {
        site.0.pos -= 1;
        site.0
    }
}

fn quicksort_by<T, F: Fn(&T, &T) -> Ordering>(data: &mut [T], cmp: &F) {
    if data.is_empty() {
        return;
    }
    let p = partition(data, cmp);
    let (data1, data2) = data.split_at_mut(p);
    quicksort_by(data1, cmp);
    quicksort_by(&mut data2[1..], cmp);
}

fn partition<T, F: Fn(&T, &T) -> Ordering>(data: &mut [T], cmp: F) -> usize {
    let pivot_idx = data.len() - 1;
    let mut i = 0;

    for j in 0..pivot_idx {
        if cmp(&data[j], &data[pivot_idx]) == Ordering::Less {
            data.swap(i, j);
            i += 1;
        }
    }
    data.swap(i, pivot_idx);
    i
}

/// Return input swapped such that ret.0 is <= than ret.1
fn order(a: usize, b: usize) -> (usize, usize) {
    if a > b {
        (b, a)
    } else {
        (a, b)
    }
}

#[cfg(test)]
mod tests {
    use rand::prelude::*;

    use super::quicksort_by;

    #[test]
    fn test_quicksort() {
        let mut rng = rand::thread_rng();
        let mut nums: Vec<usize> = (1..1000).collect();
        let nums_expected = nums.clone();
        for _ in 0..10 {
            nums.shuffle(&mut rng);
            quicksort_by(&mut nums, &|a, b| a.cmp(b));
            assert_eq!(nums_expected, nums)
        }
    }
}
