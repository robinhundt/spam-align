use std::ops::{Not, Deref};

use crate::align::Matrix;
use crate::align::micro_alignment::{MicroAlignment, Site};
use itertools::{repeat_n, Itertools};
use fxhash::FxHashMap;
use petgraph::{Graph, algo::toposort};


pub struct Closure {
    /// information about the positions in the aligned sequences
    sequences: Vec<Sequence>,
    /// this stores information about which eq class is aligned with which
    /// positions for each sequence
    /// is used to tell if a eq class is aligned with a seq (not aligned if val == 0)
    alig_set: Vec<PositionSet>,

    /// how many eq classes are there in the current iteration
    nbr_alig_sets: usize,
    /// how many eq classes are there in the previous iteration
    old_nbr_alig_sets: usize,

    /// these are prob used to save the actual frontiers,
    /// they are indexed with a eq class id and a sequence id
    pred_frontier: Matrix<usize>,
    /// these are prob used to save the actual frontiers,
    /// they are indexed with a eq class id and a sequence id
    succ_frontier: Matrix<usize>,

    /// these might be used to as buffers storing frontiers
    /// which are potentially change by adding the current positions
    left1: Vec<usize>,
    left2: Vec<usize>,
    right1: Vec<usize>,
    right2: Vec<usize>,

    /// indexed with 2 seq numbers
    pos: Matrix<usize>,
}

struct Sequence {
    length: usize,
    /// index with seq pos, stores id of alig set this site belongs to
    alig_set_nbr: Vec<usize>,
    /// next class information for pred frontier
    pred_alig_set_pos: Vec<usize>,
    /// next class information for succ frontier
    succ_alig_set_pos: Vec<usize>,
}

impl Sequence {
    fn len(&self) -> usize {
        self.length
    }
}

// TODO this can prob be replaced by a matrix which would reduce number of
// bounds checks
#[derive(Clone)]
struct PositionSet {
    pos: Vec<usize>,
}


/// Wrapper for a site whose position is increased by 1
/// since the algorithm operates in sequences with `1` being the first position
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

    /// Try to add the passed micro alignment to the partial alignment
    /// returns false if the micro alignment is inconsistent with previously
    /// added ones or if it is already included in earlier ones
    pub fn try_add_micro_alignment(&mut self, micro_alignment: &MicroAlignment) -> bool {
        // dbg!(micro_alignment);
        let mut added = false;
        let mut sites_to_add = vec![];
        for (a, b) in micro_alignment.site_pair_iter() {
            let (a, b) = (a.into(), b.into());
            if !self.alignable(a, b) {
                return false;
            }
            if !self.aligned(a, b) {
                sites_to_add.push((a, b));
            }
        }
        for (a, b) in sites_to_add {
            self.add_aligned_positions(a, b);
            added = true;
        }
        added
    }


    fn aligned(&self, a: ShiftedSite, b: ShiftedSite) -> bool {
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

    pub fn align(&self, seqs: &mut [crate::Sequence])  {
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

        let mut graph = Graph::<&Vec<Site>, ()>::with_capacity(classes.len(), 0);
        let mut node_indices = vec![];
        for eq_class in &classes {
            node_indices.push(graph.add_node(eq_class));
        }
        for idx1 in &node_indices {
            for idx2 in &node_indices {
                let class1 = graph.node_weight(*idx1).unwrap();
                let class2 = graph.node_weight(*idx2).unwrap();
                let repr1 = class1.first().expect("Empty eq class");
                let repr2 = class2.first().expect("Empty eq class");
                if self.less(*repr1, *repr2) {
                    graph.add_edge(*idx1, *idx2, ());
                }
            }
        }
        let sorted_indices = toposort(&graph, None).unwrap();
        let sorted_classes = sorted_indices
            .into_iter()
            .map(|idx| graph.node_weight(idx).unwrap().deref().clone())
            .collect_vec();
        let mut classes = sorted_classes;

        classes.reverse();
        let mut shifted_by: FxHashMap<usize, usize> = FxHashMap::default();
        while let Some(class) = classes.pop() {
            let max_pos = class
                .iter()
                .max_by_key(|site| site.pos)
                .expect("Unexpected empty eq class")
                .pos;
            shifted_by.clear();
            //            println!("Handling class {:?}", class);
            for site in &class {
                seqs[site.seq].data[site.pos].make_ascii_uppercase();
                let shift = max_pos - site.pos;
                if shift == 0 {
                    continue;
                }
                seqs[site.seq]
                    .data
                    .splice(site.pos..site.pos, repeat_n(b'-', shift));
                shifted_by.insert(site.seq, shift);
            }
            classes.iter_mut().for_each(|shifted_class| {
                shifted_class.iter_mut().for_each(|mut site| {
                    match shifted_by.get(&site.seq) {
                        Some(shift) => {
                            site.pos += *shift;
                            site
                        }
                        None => site,
                    };
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

    fn add_aligned_positions(&mut self, a: ShiftedSite, b: ShiftedSite) {
        let (a, b) = (a.0, b.0);
        let n1 = self.sequences[a.seq].alig_set_nbr[a.pos];
        let n2 = self.sequences[b.seq].alig_set_nbr[b.pos];

        // this condition is guaranteed by filtering the
        // micro alignment that is added in the `try_add_micro_alignment`
        // method
        assert!(
            n1 != n2 || n1 == 0 || n2 == 0,
            "Sites must belong to different eq classes or at least one must not be aligned"
        );

        let lookup_alig_set = |n: usize, site: Site| {
            let mut ng = n;
            let mut nd = n;
            if n == 0 {
                let mut k = self.sequences[site.seq].pred_alig_set_pos[site.pos];
                if k > 0 {
                    ng = self.sequences[site.seq].alig_set_nbr[k];
                }
                k = self.sequences[site.seq].succ_alig_set_pos[site.pos];
                if k > 0 {
                    nd = self.sequences[site.seq].alig_set_nbr[k];
                }
            }
            (ng, nd)
        };

        let (ng1, nd1) = lookup_alig_set(n1, a);
        let (ng2, nd2) = lookup_alig_set(n2, b);

        let init_left = |ng: usize, left: &mut [usize], pred_frontier: &Matrix<usize>| {
            if ng == 0 {
                left.iter_mut().for_each(|x| *x = 0);
            } else {
                let pred_frontier = pred_frontier.row(ng);
                left
                    .iter_mut()
                    .zip(pred_frontier)
                    .for_each(|(g1, pred)| *g1 = *pred);
            }

        };

        init_left(ng1, &mut self.left1, &self.pred_frontier);
        init_left(ng2, &mut self.left2, &self.pred_frontier);

        let init_right = |nd: usize, right: &mut [usize], sequences: &[Sequence], succ_frontier: &Matrix<usize>| {
            if nd == 0 {
                right
                    .iter_mut()
                    .zip(sequences)
                    .for_each(|(d1, seq)| *d1 = seq.len() + 1);
            } else {
                let succ_frontier = succ_frontier.row(nd);
                right
                    .iter_mut()
                    .zip(succ_frontier)
                    .for_each(|(d1, succ)| *d1 = *succ);
            }

        };

        init_right(nd1, &mut self.right1, &self.sequences, &self.succ_frontier);
        init_right(nd2, &mut self.right2, &self.sequences, &self.succ_frontier);


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
                let front = self.alig_set[nn].pos[x];
                self.pred_frontier[[nn, x]] = front;
                self.succ_frontier[[nn, x]] = front;
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

        // TODO dry out the following code

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

        let mut update_sequences = |n: usize, site: Site| {
            if n != 0 {
                return;
            }

            let mut k = site.pos - 1;
            while k > 0 && self.sequences[site.seq].alig_set_nbr[k] == 0 {
                self.sequences[site.seq].succ_alig_set_pos[k] = site.pos;
                k -= 1;
            }
            if k > 0 {
                self.sequences[site.seq].succ_alig_set_pos[k] = site.pos;
            }

            let mut k = site.pos + 1;
            while k <= self.sequences[site.seq].len() && self.sequences[site.seq].alig_set_nbr[k] == 0 {
                self.sequences[site.seq].pred_alig_set_pos[k] = site.pos;
                k += 1;
            }
            if k <= self.sequences[site.seq].len() {
                self.sequences[site.seq].pred_alig_set_pos[k] = site.pos
            }
        };

        update_sequences(n1, a);
        update_sequences(n2, b);

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
