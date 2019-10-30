use crate::align::gabios::PartialAlignment;
use crate::align::{Diagonal, EnumeratedSequence, Site};
use fxhash::{FxHashMap, FxHashSet};
use itertools::{repeat_n, Itertools};
use std::cmp::Ordering;
use std::iter::FromIterator;

#[derive(Debug, Clone, Default)]
pub struct EqClasses {
    classes: Vec<FxHashSet<Site>>,
}

impl EqClasses {
    pub fn new(diagonals: &[Diagonal], part_alig: &PartialAlignment) -> Self {
        let mut classes: Vec<FxHashSet<Site>> = vec![];
        for diag in diagonals {
            'outer: for (site_a, site_b) in diag.site_iter() {
                for class in &mut classes {
                    match class.iter().take(1).next() {
                        Some(site) => {
                            if part_alig.sites_are_aligned(*site, site_a) {
                                let s = site.clone();
                                class.insert(site_a);
                                class.insert(site_b);

                                if class.len() > 3 {
                                    dbg!(s, &site_a, &site_b);
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
            )
        }

        let mut unsorted_self = Self { classes };
        unsorted_self.sort_self(part_alig)
    }

    fn sort_self(mut self, part_alig: &PartialAlignment) -> Self {
        self.classes.sort_by(|a, b| {
            for site_a in a {
                for site_b in b {
                    return if part_alig.le(*site_a, *site_b) {
                        Ordering::Less
                    } else if part_alig.ge(*site_a, *site_b) {
                        Ordering::Greater
                    } else {
                        Ordering::Equal
                    };
                }
            }
            unreachable!("Called sort_self with empty eq classes!")
        });
        self
    }

    pub fn align_sequences(&self, seqs: &mut [EnumeratedSequence]) {
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
            println!("Handling class {:?}", class);
            for site in &class {
                let shift = max_pos - site.pos;
                if shift == 0 {
                    continue;
                }
                println!("Adding shift of {} at: {:?}", shift, &site);
                seqs[site.seq]
                    .seq
                    .data
                    .splice(site.pos..site.pos, repeat_n(b'-', max_pos - site.pos));
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
}
