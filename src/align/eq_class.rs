use crate::align::gabios::TransitiveClosure;
use crate::align::micro_alignment::{ScoredMicroAlignment, Site};
use crate::data_loaders::Sequence;
use fxhash::{FxHashMap, FxHashSet};
use itertools::{repeat_n, Itertools};
use petgraph::algo::toposort;
use petgraph::Graph;
use std::ops::{Deref, Not};
use std::vec::IntoIter;

#[derive(Debug, Clone, Default)]
pub struct EqClasses {
    classes: Vec<FxHashSet<Site>>,
}

impl EqClasses {
    pub fn new(diagonals: &[ScoredMicroAlignment], closure: &TransitiveClosure) -> Self {
        let mut classes: Vec<FxHashSet<Site>> = vec![];
        for diag in diagonals {
            'outer: for (site_a, site_b) in diag.site_pair_iter() {
                for class in &mut classes {
                    let all_aligned = class
                        .iter()
                        .all(|site| closure.sites_are_aligned(*site, site_a));
                    let any_aligned = class
                        .iter()
                        .any(|site| closure.sites_are_aligned(*site, site_a));

                    if any_aligned && all_aligned.not() {
                        println!("Broke transitivity again!")
                    }

                    match class.iter().take(1).next() {
                        Some(&site) => {
                            if closure.sites_are_aligned(site, site_a) {
                                class.insert(site_a);
                                class.insert(site_b);

                                if closure.sites_are_aligned(site, site_b).not() {
                                    println!(
                                        "Broke transitivity of alignment: {:?} {:?} {:?}\n {:#?}",
                                        site, site_a, site_b, diag
                                    );
                                }
                                continue 'outer;
                            }
                        }
                        None => continue,
                    }
                }
                classes.push([site_a, site_b].iter().cloned().collect())
            }
        }

        // check no multiple sites of same sequence in same class invariant
        for class in &classes {
            let unique_seq_cnt = class.iter().unique_by(|site| site.seq).count();
            assert_eq!(
                unique_seq_cnt,
                class.len(),
                "Class has sites with duplicate seq: {:?}",
                &class
            );
        }

        let unsorted_self = Self { classes };
        unsorted_self.sort_self(closure)
    }

    fn sort_self(mut self, closure: &TransitiveClosure) -> Self {
        let mut graph = Graph::<&FxHashSet<Site>, ()>::new();
        let mut node_indices = vec![];
        for eq_class in &self.classes {
            node_indices.push(graph.add_node(eq_class));
        }

        for idx1 in &node_indices {
            for idx2 in &node_indices {
                let class1 = graph.node_weight(*idx1).unwrap();
                let class2 = graph.node_weight(*idx2).unwrap();
                let repr1 = get_el_from_hashset(class1);
                let repr2 = get_el_from_hashset(class2);
                if closure.less(*repr1, *repr2) {
                    graph.add_edge(*idx1, *idx2, ());
                }
            }
        }
        let sorted_indices = toposort(&graph, None).unwrap();
        let sorted_classes = sorted_indices
            .into_iter()
            .map(|idx| graph.node_weight(idx).unwrap().deref().clone())
            .collect_vec();
        self.classes = sorted_classes;
        self
    }

    pub fn align_sequences(&self, seqs: &mut [Sequence]) {
        let mut classes = self.classes.clone();
        classes.reverse();
        while let Some(class) = classes.pop() {
            let max_pos = class
                .iter()
                .max_by_key(|site| site.pos)
                .expect("Unexpected empty eq class")
                .pos;
            let mut shifted_by: FxHashMap<usize, usize> = FxHashMap::default();
            shifted_by.reserve(class.len());
            //            println!("Handling class {:?}", class);
            for site in &class {
                seqs[site.seq].data[site.pos].make_ascii_uppercase();
                let shift = max_pos - site.pos;
                if shift == 0 {
                    continue;
                }
                //                println!("Adding shift of {} at: {:?}", shift, &site);
                seqs[site.seq]
                    .data
                    .splice(site.pos..site.pos, repeat_n(b'-', shift));
                shifted_by.insert(site.seq, shift);
            }
            classes = classes
                .into_iter()
                .map(|shifted_class| {
                    shifted_class
                        .into_iter()
                        .map(|mut site| match shifted_by.get(&site.seq) {
                            Some(shift) => {
                                site.pos += *shift;
                                site
                            }
                            None => site,
                        })
                        .collect()
                })
                .collect();
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = &FxHashSet<Site>> {
        self.classes.iter()
    }
}

fn get_el_from_hashset<E>(hashset: &FxHashSet<E>) -> &E {
    hashset.iter().next().unwrap()
}

impl IntoIterator for EqClasses {
    type Item = FxHashSet<Site>;
    type IntoIter = IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.classes.into_iter()
    }
}
