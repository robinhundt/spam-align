use crate::align::gabios::Closure as TransitiveClosure;
use crate::align::micro_alignment::{ScoredMicroAlignment, Site};
use crate::data_loaders::Sequence;
use fxhash::{FxHashMap, FxHashSet};
use itertools::{repeat_n, Itertools};
use petgraph::algo::toposort;
use petgraph::Graph;
use std::ops::{Deref, Not};
use std::time::Instant;
use std::vec::IntoIter;

#[derive(Debug, Clone, Default)]
pub struct EqClasses {
    classes: Vec<Vec<Site>>,
}

impl EqClasses {
    pub fn new(diagonals: &[ScoredMicroAlignment], closure: &TransitiveClosure) -> Self {
        let classes = closure.eq_classes();
        let unsorted_self = Self { classes };
        unsorted_self.sort_self(closure)
    }

    fn sort_self(mut self, closure: &TransitiveClosure) -> Self {
        let mut graph = Graph::<&Vec<Site>, ()>::with_capacity(self.classes.len(), 0);
        let mut node_indices = vec![];
        for eq_class in &self.classes {
            node_indices.push(graph.add_node(eq_class));
        }
        for idx1 in &node_indices {
            for idx2 in &node_indices {
                let class1 = graph.node_weight(*idx1).unwrap();
                let class2 = graph.node_weight(*idx2).unwrap();
                let repr1 = class1.first().expect("Empty eq class");
                let repr2 = class2.first().expect("Empty eq class");
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
                //                println!("Adding shift of {} at: {:?}", shift, &site);
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

    pub fn iter(&self) -> impl Iterator<Item = &Vec<Site>> {
        self.classes.iter()
    }
}

impl IntoIterator for EqClasses {
    type Item = Vec<Site>;
    type IntoIter = IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.classes.into_iter()
    }
}
