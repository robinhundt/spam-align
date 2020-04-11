use crate::align::micro_alignment::{MicroAlignment, Site};
use crate::align::Matrix;

use std::ops::Not;

pub struct Closure {
    sequences: Vec<Sequence>,
    alig_set: Vec<PositionSet>,

    nbr_alig_sets: usize,
    old_nbr_alig_sets: usize,

    pred_frontier: Matrix<usize>,
    succ_frontier: Matrix<usize>,

    left1: Vec<usize>,
    left2: Vec<usize>,
    right1: Vec<usize>,
    right2: Vec<usize>,

    pos: Matrix<usize>,
}

struct Sequence {
    length: usize,
    alig_set_nbr: Vec<usize>,
    pred_alig_set_pos: Vec<usize>,
    succ_alig_set_pos: Vec<usize>,
}

impl Sequence {
    fn len(&self) -> usize {
        self.length
    }
}

#[derive(Clone)]
struct PositionSet {
    pos: Vec<usize>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash, Default)]
struct ShiftedSite(Site);

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

impl Closure {
    pub fn new(seq_lengths: &[usize]) -> Self {
        let &max_length = seq_lengths.iter().max().expect("No sequences");
        let sequences = Closure::init_sequences(seq_lengths);
        let nbr_seqs = sequences.len();

        let pred_frontier = Matrix::zeros([max_length + 2, nbr_seqs]);
        let succ_frontier = Matrix::zeros([max_length + 2, nbr_seqs]);
        let alig_set = vec![
            PositionSet {
                pos: vec![0; nbr_seqs]
            };
            max_length + 2
        ];

        let left1 = vec![0; nbr_seqs];
        let left2 = vec![0; nbr_seqs];
        let right1 = vec![0; nbr_seqs];
        let right2 = vec![0; nbr_seqs];

        let pos = Matrix::zeros([nbr_seqs, nbr_seqs]);
        Self {
            sequences,
            pred_frontier,
            succ_frontier,
            alig_set,
            left1,
            nbr_alig_sets: 0,
            old_nbr_alig_sets: 0,
            left2,
            right1,
            right2,
            pos,
        }
    }

    fn init_sequences(long_seq: &[usize]) -> Vec<Sequence> {
        long_seq
            .iter()
            .map(|len| {
                let padded_len = len + 2;
                Sequence {
                    length: *len,
                    alig_set_nbr: vec![0; padded_len],
                    pred_alig_set_pos: vec![0; padded_len],
                    succ_alig_set_pos: vec![0; padded_len],
                }
            })
            .collect()
    }

    pub fn try_add_micro_alignment(&mut self, micro_alignment: &MicroAlignment) -> bool {
        // dbg!(micro_alignment);
        let mut added = false;
        let mut sites_to_add = vec![];
        for (a, b) in micro_alignment.site_pair_iter() {
            let (a, b) = (a.into(), b.into());
            if !self.alignable(a, b) {
                return false;
            }
            if !self.aligned_shifted(a, b) {
                sites_to_add.push((a, b));
            }
        }
        for (a, b) in sites_to_add {
            self.add_aligned_positions(a, b);
            added = true;
        }
        added
    }

    pub fn aligned(&self, a: Site, b: Site) -> bool {
        let (a, b) = (a.into(), b.into());
        self.aligned_shifted(a, b)
    }

    fn aligned_shifted(&self, a: ShiftedSite, b: ShiftedSite) -> bool {
        let (a, b) = (a.0, b.0);
        a.seq == b.seq && a.pos == b.pos
            || self.sequences[a.seq].alig_set_nbr[a.pos] != 0
                && self.sequences[a.seq].alig_set_nbr[a.pos]
                    == self.sequences[b.seq].alig_set_nbr[b.pos]
    }

    pub fn less(&self, left: Site, right: Site) -> bool {
        let (left, right) = (left.into(), right.into());
        self.path_exists(left, right)
    }

    fn alignable(&self, a: ShiftedSite, b: ShiftedSite) -> bool {
        if self.path_exists(a, b) {
            self.path_exists(b, a)
        } else {
            self.path_exists(b, a).not()
        }
    }

    fn path_exists(&self, a: ShiftedSite, b: ShiftedSite) -> bool {
        let (a, b) = (a.0, b.0);
        if a.seq == b.seq {
            return a.pos < b.pos;
        }
        let mut n2 = self.sequences[b.seq].alig_set_nbr[b.pos];
        if n2 == 0 {
            let k = self.sequences[b.seq].pred_alig_set_pos[b.pos];
            if k > 0 {
                n2 = self.sequences[b.seq].alig_set_nbr[k];
            }
        }
        if n2 == 0 {
            false
        } else {
            a.pos <= self.pred_frontier[[n2, a.seq]]
        }
    }

    pub fn eq_classes(&self) -> Vec<Vec<Site>> {
        let mut classes = vec![vec![]; self.nbr_alig_sets];
        for (idx, seq) in self.sequences.iter().enumerate() {
            for (pos, alig_set) in seq.alig_set_nbr.iter().enumerate().skip(1) {
                if *alig_set == 0 {
                    continue;
                }
                let site: Site = ShiftedSite(Site { seq: idx, pos }).into();
                let alig_set = alig_set - 1;
                classes[alig_set].push(site);
            }
        }
        classes
    }

    // fn pred_frontier(&self, origin_site: ShiftedSite, target_seq: usize) -> usize {
    //     let s = origin_site.0;
    //     let mut n = self.sequences[s.seq].alig_set_nbr[s.pos];
    //     if n == 0 {
    //         let k = self.sequences[s.seq].pred_alig_set_pos[s.pos];
    //         if k > 0 {
    //             n = self.sequences[s.seq].alig_set_nbr[k];
    //         }
    //     }
    //     if n > 0 {
    //         self.pred_frontier[[n, target_seq]]
    //     } else {
    //         0
    //     }
    // }
    //
    // fn succ_frontier(&self, origin_site: ShiftedSite, target_seq: usize) -> usize {
    //     let s = origin_site.0;
    //     let mut n = self.sequences[s.seq].alig_set_nbr[s.pos];
    //     if n == 0 {
    //         let k = self.sequences[s.seq].succ_alig_set_pos[s.pos];
    //         if k > 0 {
    //             n = self.sequences[s.seq].alig_set_nbr[k];
    //         }
    //     }
    //     if n > 0 {
    //         self.succ_frontier[[n, target_seq]]
    //     } else {
    //         self.sequences[target_seq].len() + 1
    //     }
    // }

    fn add_aligned_positions(&mut self, a: ShiftedSite, b: ShiftedSite) {
        let (a, b) = (a.0, b.0);
        let n1 = self.sequences[a.seq].alig_set_nbr[a.pos];
        let mut ng1 = n1;
        let mut nd1 = n1;

        let n2 = self.sequences[b.seq].alig_set_nbr[b.pos];
        let mut ng2 = n2;
        let mut nd2 = n2;

        if !(n1 == 0 || n2 == 0 || n1 != n2) {
            return;
        }
        if ng1 == 0 {
            let mut k = self.sequences[a.seq].pred_alig_set_pos[a.pos];
            if k > 0 {
                ng1 = self.sequences[a.seq].alig_set_nbr[k];
            }
            k = self.sequences[a.seq].succ_alig_set_pos[a.pos];
            if k > 0 {
                nd1 = self.sequences[a.seq].alig_set_nbr[k];
            }
        }

        if ng2 == 0 {
            let mut k = self.sequences[b.seq].pred_alig_set_pos[b.pos];
            if k > 0 {
                ng2 = self.sequences[b.seq].alig_set_nbr[k];
            }
            k = self.sequences[b.seq].succ_alig_set_pos[b.pos];
            if k > 0 {
                nd2 = self.sequences[b.seq].alig_set_nbr[k];
            }
        }

        if ng1 == 0 {
            self.left1.iter_mut().for_each(|x| *x = 0);
        } else {
            assert!(
                ng1 < self.pred_frontier.rows(),
                format!("{}, {}", ng1, self.pred_frontier.rows())
            );
            let pred_frontier = self.pred_frontier.row(ng1);
            self.left1
                .iter_mut()
                .zip(pred_frontier)
                .for_each(|(g1, pred)| *g1 = *pred);
        }

        if nd1 == 0 {
            self.right1
                .iter_mut()
                .zip(&self.sequences)
                .for_each(|(d1, seq)| *d1 = seq.len() + 1);
        } else {
            let succ_frontier = self.succ_frontier.row(nd1);
            self.right1
                .iter_mut()
                .zip(succ_frontier)
                .for_each(|(d1, succ)| *d1 = *succ);
        }

        if ng2 == 0 {
            self.left2.iter_mut().for_each(|g2| *g2 = 0);
        } else {
            let pred_frontier = self.pred_frontier.row(ng2);
            self.left2
                .iter_mut()
                .zip(pred_frontier)
                .for_each(|(g2, pred)| *g2 = *pred);
        }

        if nd2 == 0 {
            self.right2
                .iter_mut()
                .zip(&self.sequences)
                .for_each(|(d1, seq)| *d1 = seq.len() + 1);
        } else {
            let succ_frontier = self.succ_frontier.row(nd2);
            self.right2
                .iter_mut()
                .zip(succ_frontier)
                .for_each(|(d2, succ)| *d2 = *succ);
        }

        self.left1[a.seq] = a.pos;
        self.right1[a.seq] = a.pos;

        self.left2[b.seq] = b.pos;
        self.right2[b.seq] = b.pos;

        let nn = self.nbr_alig_sets + 1;
        for x in 0..self.sequences.len() {
            self.alig_set[nn].pos[x] = 0;
            if n1 > 0 && self.alig_set[n1].pos[x] > 0 {
                self.alig_set[nn].pos[x] = self.alig_set[n1].pos[x];
            } else if n2 > 0 && self.alig_set[n2].pos[x] > 0 {
                self.alig_set[nn].pos[x] = self.alig_set[n2].pos[x];
            }

            if self.alig_set[nn].pos[x] == 0 {
                self.pred_frontier[[nn, x]] = std::cmp::max(self.left1[x], self.left2[x]);
                self.succ_frontier[[nn, x]] = std::cmp::min(self.right1[x], self.right2[x]);
            } else {
                self.pred_frontier[[nn, x]] = self.alig_set[nn].pos[x];
                self.succ_frontier[[nn, x]] = self.alig_set[nn].pos[x];
            }
        }
        self.pred_frontier[[nn, a.seq]] = a.pos;
        self.succ_frontier[[nn, a.seq]] = a.pos;
        self.alig_set[nn].pos[a.seq] = a.pos;

        self.pred_frontier[[nn, b.seq]] = b.pos;
        self.succ_frontier[[nn, b.seq]] = b.pos;
        self.alig_set[nn].pos[b.seq] = b.pos;

        for x in 0..self.sequences.len() {
            let k = self.alig_set[nn].pos[x];
            if k > 0 {
                self.sequences[x].alig_set_nbr[k] = nn;
            }
        }

        for x in 0..self.sequences.len() {
            if self.right1[x] == self.right2[x] {
                continue;
            }
            for y in 0..self.sequences.len() {
                self.pos[[x, y]] = 0;
                let mut k = self.succ_frontier[[nn, x]];
                if k == self.alig_set[nn].pos[x] {
                    k = self.sequences[x].succ_alig_set_pos[k];
                }
                if k <= self.sequences[x].len() {
                    while k > 0 {
                        let n = self.sequences[x].alig_set_nbr[k];
                        if self.pred_frontier[[n, y]] < self.pred_frontier[[nn, y]] {
                            self.pos[[x, y]] = k;
                            k = self.sequences[x].succ_alig_set_pos[k];
                        } else {
                            break;
                        }
                    }
                }
            }
        }

        for x in 0..self.sequences.len() {
            if self.right1[x] == self.right2[x] {
                continue;
            }
            for y in 0..self.sequences.len() {
                let mut k = self.succ_frontier[[nn, x]];
                if k == self.alig_set[nn].pos[x] {
                    k = self.sequences[x].succ_alig_set_pos[k];
                }
                if self.pos[[x, y]] > 0 {
                    while k > 0 && k <= self.pos[[x, y]] {
                        let n = self.sequences[x].alig_set_nbr[k];
                        self.pred_frontier[[n, y]] = self.pred_frontier[[nn, y]];
                        k = self.sequences[x].succ_alig_set_pos[k];
                    }
                }
            }
        }

        for x in 0..self.sequences.len() {
            if self.left1[x] == self.left2[x] {
                continue;
            }
            for y in 0..self.sequences.len() {
                self.pos[[x, y]] = 0;
                let mut k = self.pred_frontier[[nn, x]];
                if k > 0 && k == self.alig_set[nn].pos[x] {
                    k = self.sequences[x].pred_alig_set_pos[k];
                }
                while k > 0 {
                    let n = self.sequences[x].alig_set_nbr[k];
                    if self.succ_frontier[[n, y]] > self.succ_frontier[[nn, y]] {
                        self.pos[[x, y]] = k;
                        k = self.sequences[x].pred_alig_set_pos[k];
                    } else {
                        break;
                    }
                }
            }
        }

        for x in 0..self.sequences.len() {
            if self.left1[x] == self.left2[x] {
                continue;
            }
            for y in 0..self.sequences.len() {
                let mut k = self.pred_frontier[[nn, x]];
                if k > 0 && k == self.alig_set[nn].pos[x] {
                    k = self.sequences[x].pred_alig_set_pos[k];
                }
                if self.pos[[x, y]] > 0 {
                    while k >= self.pos[[x, y]] {
                        let n = self.sequences[x].alig_set_nbr[k];
                        self.succ_frontier[[n, y]] = self.succ_frontier[[nn, y]];
                        k = self.sequences[x].pred_alig_set_pos[k];
                    }
                }
            }
        }

        if n1 == 0 {
            let mut k = a.pos - 1;
            while k > 0 && self.sequences[a.seq].alig_set_nbr[k] == 0 {
                self.sequences[a.seq].succ_alig_set_pos[k] = a.pos;
                k -= 1;
            }
            if k > 0 {
                self.sequences[a.seq].succ_alig_set_pos[k] = a.pos;
            }

            let mut k = a.pos + 1;
            while k <= self.sequences[a.seq].len() && self.sequences[a.seq].alig_set_nbr[k] == 0 {
                self.sequences[a.seq].pred_alig_set_pos[k] = a.pos;
                k += 1;
            }
            if k <= self.sequences[a.seq].len() {
                self.sequences[a.seq].pred_alig_set_pos[k] = a.pos
            }
        }

        if n2 == 0 {
            let mut k = b.pos - 1;
            while k > 0 && self.sequences[b.seq].alig_set_nbr[k] == 0 {
                self.sequences[b.seq].succ_alig_set_pos[k] = b.pos;
                k -= 1;
            }
            if k > 0 {
                self.sequences[b.seq].succ_alig_set_pos[k] = b.pos;
            }

            let mut k = b.pos + 1;
            while k <= self.sequences[b.seq].len() && self.sequences[b.seq].alig_set_nbr[k] == 0 {
                self.sequences[b.seq].pred_alig_set_pos[k] = b.pos;
                k += 1;
            }
            if k <= self.sequences[b.seq].len() {
                self.sequences[b.seq].pred_alig_set_pos[k] = b.pos
            }
        }

        let (n1, n2) = if n1 > n2 { (n2, n1) } else { (n1, n2) };

        if n2 == 0 {
            self.nbr_alig_sets += 1;
            self.realloc();
        } else if n1 == 0 {
            self.move_alig_set(n2, nn)
        } else {
            self.move_alig_set(n1, nn);
            if n2 < self.nbr_alig_sets {
                self.move_alig_set(n2, self.nbr_alig_sets)
            }
            self.nbr_alig_sets -= 1;
            self.realloc();
        }
    }

    fn move_alig_set(&mut self, n1: usize, n2: usize) {
        for x in 0..self.sequences.len() {
            let k = self.alig_set[n2].pos[x];
            self.alig_set[n1].pos[x] = k;
            if k > 0 {
                self.sequences[x].alig_set_nbr[k] = n1;
            }
            self.pred_frontier[[n1, x]] = self.pred_frontier[[n2, x]];
            self.succ_frontier[[n1, x]] = self.succ_frontier[[n2, x]];
        }
    }

    fn realloc(&mut self) {
        if self.nbr_alig_sets <= self.old_nbr_alig_sets {
            return;
        }
        let new_rows = self.nbr_alig_sets + 2;
        self.pred_frontier.resize_rows(new_rows);
        self.succ_frontier.resize_rows(new_rows);
        self.alig_set.resize(
            self.nbr_alig_sets + 2,
            PositionSet {
                pos: vec![0; self.sequences.len()],
            },
        );
        self.old_nbr_alig_sets = self.nbr_alig_sets;
    }
}
