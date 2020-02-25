use crate::align::micro_alignment::Site;
use bio::data_structures::bwt::less;
use ndarray::{Array2, NdIndex};
use std::cmp::Ordering;

#[derive(Debug, Clone, PartialEq)]
pub(super) struct UnionFind {
    parent: Array2<Site>,

    rank: Array2<u8>,
}

impl UnionFind {
    pub(super) fn new(max_row: usize, max_col: usize) -> Self {
        let rank = Array2::zeros((max_row, max_col));
        let mut parent = Array2::default((max_row, max_col));
        parent
            .indexed_iter_mut()
            .for_each(|((seq, pos), el): (_, &mut Site)| {
                el.seq = seq;
                el.pos = pos;
            });
        Self { rank, parent }
    }

    pub(super) fn find(&self, mut site: Site) -> Site {
        unsafe {
            loop {
                let site_parent = self.parent[<[usize; 2]>::from(site)];
                if site_parent == site {
                    break;
                }
                site = site_parent;
            }
            site
        }
    }

    pub(super) fn find_mut(&mut self, mut site: Site) -> Site {
        let mut parent = self.parent[<[usize; 2]>::from(site)];
        while parent != site {
            let grandparent = self.parent[<[usize; 2]>::from(parent)];
            self.parent[<[usize; 2]>::from(site)] = grandparent;
            site = parent;
            parent = grandparent;
        }
        site
    }

    pub(super) fn union(&mut self, a: Site, b: Site) -> bool {
        if a == b {
            return false;
        }
        let a_rep = self.find_mut(a);
        let b_rep = self.find_mut(b);

        if a_rep == b_rep {
            return false;
        }

        let a_rank = self.rank[<[usize; 2]>::from(a_rep)];
        let b_rank = self.rank[<[usize; 2]>::from(b_rep)];

        match a_rank.cmp(&b_rank) {
            Ordering::Less => self.parent[<[usize; 2]>::from(a_rep)] = b_rep,
            Ordering::Greater => self.parent[<[usize; 2]>::from(b_rep)] = a_rep,
            Ordering::Equal => {
                self.parent[<[usize; 2]>::from(b_rep)] = a_rep;
                self.rank[<[usize; 2]>::from(a_rep)] += 1;
            }
        }
        true
    }
}

impl From<Site> for [usize; 2] {
    fn from(s: Site) -> Self {
        [s.seq, s.pos]
    }
}

////! `UnionFind<K>` is a disjoint-set data structure.
//
// use super::graph::IndexType;
// use std::cmp::Ordering;
//
// /// `UnionFind<K>` is a disjoint-set data structure. It tracks set membership of *n* elements
// /// indexed from *0* to *n - 1*. The scalar type is `K` which must be an unsigned integer type.
// ///
// /// <http://en.wikipedia.org/wiki/Disjoint-set_data_structure>
// ///
// /// Too awesome not to quote:
// ///
// /// “The amortized time per operation is **O(α(n))** where **α(n)** is the
// /// inverse of **f(x) = A(x, x)** with **A** being the extremely fast-growing Ackermann function.”
// #[derive(Debug, Clone)]
// pub struct UnionFind<K> {
//     // For element at index *i*, store the index of its parent; the representative itself
//     // stores its own index. This forms equivalence classes which are the disjoint sets, each
//     // with a unique representative.
//     parent: Vec<K>,
//     // It is a balancing tree structure,
//     // so the ranks are logarithmic in the size of the container -- a byte is more than enough.
//     //
//     // Rank is separated out both to save space and to save cache in when searching in the parent
//     // vector.
//     rank: Vec<u8>,
// }
//
// #[inline]
// unsafe fn get_unchecked<K>(xs: &[K], index: usize) -> &K {
//     debug_assert!(index < xs.len());
//     xs.get_unchecked(index)
// }
//
// #[inline]
// unsafe fn get_unchecked_mut<K>(xs: &mut [K], index: usize) -> &mut K {
//     debug_assert!(index < xs.len());
//     xs.get_unchecked_mut(index)
// }
//
// impl<K> UnionFind<K>
// where
//     K: IndexType,
// {
//     /// Create a new `UnionFind` of `n` disjoint sets.
//     pub fn new(n: usize) -> Self {
//         let rank = vec![0; n];
//         let parent = (0..n).map(K::new).collect::<Vec<K>>();
//
//         UnionFind { parent, rank }
//     }
//
//     /// Return the representative for `x`.
//     ///
//     /// **Panics** if `x` is out of bounds.
//     pub fn find(&self, x: K) -> K {
//         assert!(x.index() < self.parent.len());
//         unsafe {
//             let mut x = x;
//             loop {
//                 // Use unchecked indexing because we can trust the internal set ids.
//                 let xparent = *get_unchecked(&self.parent, x.index());
//                 if xparent == x {
//                     break;
//                 }
//                 x = xparent;
//             }
//             x
//         }
//     }
//
//     /// Return the representative for `x`.
//     ///
//     /// Write back the found representative, flattening the internal
//     /// datastructure in the process and quicken future lookups.
//     ///
//     /// **Panics** if `x` is out of bounds.
//     pub fn find_mut(&mut self, x: K) -> K {
//         assert!(x.index() < self.parent.len());
//         unsafe { self.find_mut_recursive(x) }
//     }
//
//     unsafe fn find_mut_recursive(&mut self, mut x: K) -> K {
//         let mut parent = *get_unchecked(&self.parent, x.index());
//         while parent != x {
//             let grandparent = *get_unchecked(&self.parent, parent.index());
//             *get_unchecked_mut(&mut self.parent, x.index()) = grandparent;
//             x = parent;
//             parent = grandparent;
//         }
//         x
//     }
//
//     /// Returns `true` if the given elements belong to the same set, and returns
//     /// `false` otherwise.
//     pub fn equiv(&self, x: K, y: K) -> bool {
//         self.find(x) == self.find(y)
//     }
//
//     /// Unify the two sets containing `x` and `y`.
//     ///
//     /// Return `false` if the sets were already the same, `true` if they were unified.
//     ///
//     /// **Panics** if `x` or `y` is out of bounds.
//     pub fn union(&mut self, x: K, y: K) -> bool {
//         if x == y {
//             return false;
//         }
//         let xrep = self.find_mut(x);
//         let yrep = self.find_mut(y);
//
//         if xrep == yrep {
//             return false;
//         }
//
//         let xrepu = xrep.index();
//         let yrepu = yrep.index();
//         let xrank = self.rank[xrepu];
//         let yrank = self.rank[yrepu];
//
//         // The rank corresponds roughly to the depth of the treeset, so put the
//         // smaller set below the larger
//         match xrank.cmp(&yrank) {
//             Ordering::Less => self.parent[xrepu] = yrep,
//             Ordering::Greater => self.parent[yrepu] = xrep,
//             Ordering::Equal => {
//                 self.parent[yrepu] = xrep;
//                 self.rank[xrepu] += 1;
//             }
//         }
//         true
//     }
//
//     /// Return a vector mapping each element to its representative.
//     pub fn into_labeling(mut self) -> Vec<K> {
//         // write in the labeling of each element
//         unsafe {
//             for ix in 0..self.parent.len() {
//                 let k = *get_unchecked(&self.parent, ix);
//                 let xrep = self.find_mut_recursive(k);
//                 *self.parent.get_unchecked_mut(ix) = xrep;
//             }
//         }
//         self.parent
//     }
// }
