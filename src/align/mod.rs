use itertools::{max, repeat_n, Itertools};
use ndarray::Array2;

use crate::align::gabios::PartialAlignment;
use crate::data_loaders::Sequence;
use crate::score::score_prot_pairwise;
use crate::spaced_word::{
    find_word_matches, generate_random_patterns, read_patterns_from_file, Pattern,
};
use fxhash::hash;
use std::iter::{FromIterator, Map, Zip};
use std::ops::RangeInclusive;

pub mod diagonal;
pub mod eq_class;
pub mod gabios;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Diagonal {
    /// end point (inclusive) in first sequence
    i: Site,
    /// end point (inclusive) in second sequence
    j: Site,
    /// diagonal length
    k: usize,
}

impl Diagonal {
    fn overlaps_with(&self, other: &Diagonal) -> bool {
        let overlap = |end_site1: Site, end_site2: Site| {
            end_site1.seq == end_site2.seq
                && end_site1.pos - self.k <= end_site2.pos
                && end_site2.pos - other.k <= end_site1.pos
        };

        overlap(self.i, other.i)
            || overlap(self.i, other.j)
            || overlap(self.j, other.i)
            || overlap(self.j, other.j)
    }

    fn site_iter<'a>(&'a self) -> impl Iterator<Item = (Site, Site)> + 'a {
        (self.i.pos - self.k..=self.i.pos)
            .zip((self.j.pos - self.k..=self.j.pos))
            .map(move |(pos_i, pos_j)| {
                let site_i = Site {
                    pos: pos_i,
                    seq: self.i.seq,
                };
                let site_j = Site {
                    pos: pos_j,
                    seq: self.j.seq,
                };
                (site_i, site_j)
            })
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub struct Site {
    seq: usize,
    pos: usize,
}

#[derive(Debug, Clone)]
pub struct EnumeratedSequence {
    idx: usize,
    seq: Sequence,
}

impl EnumeratedSequence {
    pub fn new(idx: usize, seq: Sequence) -> Self {
        Self { idx, seq }
    }
}

pub fn align(sequences: &[EnumeratedSequence]) -> (Vec<Diagonal>, PartialAlignment) {
    //    let patterns = generate_random_patterns(5, 3, 10, 0, 3);
    let patterns = read_patterns_from_file("sample.pat").unwrap();
    dbg!(&patterns);
    let mut diagonals = find_diagonals(sequences, &patterns).collect_vec();
    let max_seq_len = sequences
        .iter()
        .max_by_key(|seq| seq.seq.data.len())
        .unwrap()
        .seq
        .data
        .len();

    println!("# of diagonals: {}", diagonals.len());

    let mut partial_alignment = PartialAlignment::new(max_seq_len, sequences.len());

    diagonals.sort_by_cached_key(|diag| {
        let sequences = [&sequences[diag.i.seq], &sequences[diag.j.seq]];
        -weight(diag, sequences)
    });
    dbg!(hash(&diagonals));

    let diagonals = diagonals
        .into_iter()
        .enumerate()
        .filter(|(idx, diag)| partial_alignment.add_diagonal(diag))
        .map(|(idx, diag)| {
            println!("Added: {:?}", diag);
            diag
        })
        .collect();
    (diagonals, partial_alignment)
    // sort diagonals by weight
    // if diagonal consistent -> add to alignment
    // return added diagonals
}

pub fn align_seq_pair(sequences: [&EnumeratedSequence; 2]) -> Vec<Diagonal> {
    let patterns = generate_random_patterns(30, 3, 10, 0, 3);
    let diagonals: Vec<Diagonal> = dbg!(find_diagonals_in_pair(sequences, &patterns).collect());
    let len1 = sequences[0].seq.data.len();
    let len2 = sequences[1].seq.data.len();

    let mut prec: Array2<i32> = Array2::from_elem([len1, len2], -1);

    let mut sigma_d: Vec<i32> = Vec::with_capacity(diagonals.len());
    let mut pi_d: Vec<i32> = Vec::with_capacity(diagonals.len());
    let mut score: Array2<i32> = Array2::zeros([len1, len2]);
    let mut sigma_d_ij: Array2<usize> = Array2::zeros([len1, len2]);

    for i in 1..len1 {
        for j in 1..len2 {
            score[[i, j]] = calc_score_ij(i, j, sequences, &diagonals, &score);
            prec[[i, j]] = calc_prec_ij(i, j, sequences, &diagonals, &score, &prec);
        }
    }

    // TODO error handling when no diagonals with positive score exist
    let mut aligned_diagonals = match prec[[len1 - 1, len2 - 1]] {
        -1 => return Vec::new(),
        diag => vec![diag],
    };

    while let Some(diag_idx) = calc_pi(
        &diagonals[dbg!(*aligned_diagonals.last().unwrap()) as usize],
        &prec,
    ) {
        let diag = &diagonals[diag_idx as usize];
        aligned_diagonals.push(prec[[diag.i.pos, diag.j.pos]]);
    }

    aligned_diagonals
        .into_iter()
        .map(|diag_idx| diagonals[diag_idx as usize].clone())
        .rev()
        .collect()
}

fn find_diagonals_in_pair<'a>(
    sequences: [&'a EnumeratedSequence; 2],
    patterns: &'a [Pattern],
) -> impl Iterator<Item = Diagonal> + 'a {
    patterns.into_iter().flat_map(move |pattern| {
        find_word_matches(
            &pattern,
            &[sequences[0].seq.clone(), sequences[1].seq.clone()],
        )
        .map(move |word_match| {
            let i = Site {
                seq: sequences[0].idx,
                pos: word_match.fst_word.pos_in_seq + pattern.len() - 1,
            };
            let j = Site {
                seq: sequences[1].idx,
                pos: word_match.snd_word.pos_in_seq + pattern.len() - 1,
            };
            Diagonal {
                i,
                j,
                k: pattern.len() - 1,
            }
        })
    })
}

fn find_diagonals<'a>(
    sequences: &'a [EnumeratedSequence],
    patterns: &'a [Pattern],
) -> impl Iterator<Item = Diagonal> + 'a {
    sequences
        .iter()
        .tuple_combinations::<(_, _)>()
        .flat_map(move |(s1, s2)| find_diagonals_in_pair([s1, s2], &patterns))
}

fn calc_sigma_d(sequences: [&EnumeratedSequence; 2], diag: &Diagonal, score: &Array2<i32>) -> i32 {
    if (diag.i.pos - diag.k) as i32 - 1 < 0 || (diag.j.pos - diag.k) as i32 - 1 < 0 {
        return weight(diag, sequences);
    }
    score[[diag.i.pos - diag.k - 1, diag.j.pos - diag.k - 1]] + weight(diag, sequences)
}

fn calc_sigma_d_ij(
    i: usize,
    j: usize,
    sequences: [&EnumeratedSequence; 2],
    diagonals: &[Diagonal],
    score: &Array2<i32>,
) -> Option<(usize, i32)> {
    diagonals
        .iter()
        .enumerate()
        .filter(|(_, diag)| diag.i.pos == i && diag.j.pos == j)
        .map(|(diag_cnt, diag)| (diag_cnt, calc_sigma_d(sequences, diag, score)))
        .fold(
            vec![],
            |mut max_set: Vec<(usize, i32)>, (diag_cnt, sigma_diag)| {
                match max_set.get(0) {
                    Some(&(max_diag_cnt, max_sigma_diag)) => {
                        if sigma_diag > max_sigma_diag {
                            max_set.clear();
                            max_set.push((diag_cnt, sigma_diag));
                        }
                    }
                    None => max_set.push((diag_cnt, sigma_diag)),
                };
                max_set
            },
        )
        .into_iter()
        .max_by_key(|(diag_cnt, sigma_diag)| diagonals[*diag_cnt].k)
}

fn calc_score_ij(
    i: usize,
    j: usize,
    sequences: [&EnumeratedSequence; 2],
    diagonals: &[Diagonal],
    score: &Array2<i32>,
) -> i32 {
    *max(&[
        score[[i, j - 1]],
        score[[i - 1, j]],
        calc_sigma_d_ij(i, j, sequences, diagonals, score)
            .unwrap_or((0, 0))
            .1,
    ])
    .expect("No maximum in calc_score_ij")
}

fn calc_prec_ij(
    i: usize,
    j: usize,
    sequences: [&EnumeratedSequence; 2],
    diagonals: &[Diagonal],
    score: &Array2<i32>,
    prec: &Array2<i32>,
) -> i32 {
    if score[[i, j]] == score[[i, j - 1]] {
        prec[[i, j - 1]]
    } else if score[[i, j - 1]] < score[[i, j]] && score[[i, j]] == score[[i - 1, j]] {
        prec[[i - 1, j]]
    } else {
        calc_sigma_d_ij(i, j, sequences, diagonals, score)
            .unwrap_or_else(|| panic!("No diagonals ending in i:{}, j:{}, this is a bug", i, j))
            .0 as i32
    }
}

fn weight(diagonal: &Diagonal, sequences: [&EnumeratedSequence; 2]) -> i32 {
    let diag2 = &sequences[1].seq.data[diagonal.j.pos - diagonal.k..=diagonal.j.pos];
    let diag1 = &sequences[0].seq.data[diagonal.i.pos - diagonal.k..=diagonal.i.pos];
    score_prot_pairwise(diag1, diag2)
}

fn calc_pi(diag: &Diagonal, prec: &Array2<i32>) -> Option<usize> {
    if (diag.i.pos - diag.k) as i64 - 1 < 0 || (diag.j.pos - diag.k) as i64 - 1 < 0 {
        return None;
    }
    match prec[[diag.i.pos - diag.k - 1, diag.j.pos - diag.k - 1]] {
        -1 => None,
        prec_idx => Some(prec_idx as usize),
    }
}

fn construct_pair_alignment_from_diagonals(
    diagonals: &[Diagonal],
    sequences: [&Sequence; 2],
) -> [Sequence; 2] {
    let s1 = &sequences[0].data;
    let s2 = &sequences[1].data;
    let mut s1_align = Vec::with_capacity(sequences[0].data.len());
    let mut s2_align = Vec::with_capacity(sequences[1].data.len());

    let mut i = 0;
    let mut j = 0;
    for diag in diagonals {
        loop {
            if i == diag.i.pos - diag.k && j < diag.j.pos - diag.k {
                s1_align.push(b'-');
                s2_align.push(s2[j]);
                j += 1;
            } else if i < diag.i.pos - diag.k && j == diag.j.pos - diag.k {
                s1_align.push(s1[i]);
                s2_align.push(b'-');
                i += 1;
            } else {
                s1_align.push(s1[i]);
                s2_align.push(s2[j]);
                i += 1;
                j += 1;
            }
            if i == diag.i.pos + 1 && j == diag.j.pos + 1 {
                break;
            }
        }
    }
    loop {
        if i == s1.len() && j == s2.len() {
            break;
        }
        if i == s1.len() && j < s2.len() {
            s1_align.push(b'-');
            s2_align.push(s2[j]);
            j += 1;
        } else if i < s1.len() && j == s2.len() {
            s1_align.push(s1[i]);
            s2_align.push(b'-');
            i += 1;
        } else {
            s1_align.push(s1[i]);
            s2_align.push(s2[j]);
            i += 1;
            j += 1;
        }
    }

    [
        Sequence::new(sequences[0].name.clone(), s1_align),
        Sequence::new(sequences[1].name.clone(), s2_align),
    ]
}

//pub fn insert_gaps(sequences: &mut [EnumeratedSequence], diagonals: Vec<Diagonal>) {
//    let mut connected_diagonals = ConnectedDiagonals::new(diagonals);
//    connected_diagonals.sort();
//    connected_diagonals.reverse();
//
//    //    let mut added_gaps = Vec::from_iter(repeat_n(0, sequences.len()));
//
//    while let Some(mut component) = connected_diagonals.pop() {
//        //        added_gaps.iter_mut().for_each(|gap_count| *gap_count = 0);
//        while let Some(diagonal) = component.pop() {
//            let max_end_pos = component.iter().fold(
//                std::cmp::max(diagonal.i.pos, diagonal.j.pos),
//                |acc, diag| *max(&[acc, diag.i.pos, diag.j.pos]).unwrap(),
//            );
//            if max_end_pos == 0 {
//                panic!(
//                    "Bug: Empty component in insert_gaps; len: {}",
//                    component.len()
//                )
//            }
//            let mut splice_seq = |end_site: Site, k: usize| {
//                let gap_start_idx = end_site.pos - k;
//                let gap_count = max_end_pos - end_site.pos;
//                let gap_iter = repeat_n(b'-', gap_count);
//                dbg!(gap_count);
//                sequences[end_site.seq]
//                    .seq
//                    .data
//                    .splice(gap_start_idx..gap_start_idx, gap_iter);
//
//                let update_diagonal = |diagonal: &mut Diagonal| {
//                    if diagonal.i.seq == end_site.seq {
//                        diagonal.i.pos += gap_count;
//                    } else if diagonal.j.seq == end_site.seq {
//                        diagonal.j.pos += gap_count;
//                    }
//                };
//                for other_component in connected_diagonals.iter_mut() {
//                    for diagonal in other_component.iter_mut() {
//                        update_diagonal(diagonal);
//                    }
//                }
//                for diagonal in component.iter_mut() {
//                    update_diagonal(diagonal);
//                }
//            };
//            splice_seq(diagonal.i, diagonal.k);
//            splice_seq(diagonal.j, diagonal.k);
//        }
//    }
//}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_loaders::balibase;
    use bio::alignment::sparse::find_kmer_matches_seq1_hashed;
    use eq_class::EqClasses;

    #[test]
    fn test_align_seq_pair() {
        let seq1 = Sequence {
            name: String::from("s1"),
            data: "HHHHHHEEEEEEEEEEEMMMMMMMMMMMSSSSSSSLLLL".into(),
        };
        let seq2 = Sequence {
            name: String::from("s2"),
            data: "VVVVHHHHTTTTTTEEEEEEEEEPPPPPPMMMMMMMMMMLLLLLL".into(),
        };
        let seq1 = EnumeratedSequence { idx: 0, seq: seq1 };
        let seq2 = EnumeratedSequence { idx: 1, seq: seq2 };
        let diags = align_seq_pair([&seq1, &seq2]);

        dbg!(construct_pair_alignment_from_diagonals(
            &diags,
            [&seq1.seq, &seq2.seq]
        ));
        panic!();
    }

    #[test]
    fn test_align() {
        let alignment = balibase::parse_xml_file("data/bb3_release/RV30/BBS30004.xml").unwrap();
        let mut sequences = alignment
            .unaligned_data
            .into_iter()
            .enumerate()
            .map(|(idx, seq)| EnumeratedSequence { idx, seq })
            .collect_vec();
        //        let s1 = Sequence {
        //            name: String::from("s1"),
        //            data: "HHHHHHEEEEEEEEEEEMMMMMMMMMMMSSSSSSSLLLL".into(),
        //        };
        //        let s2 = Sequence {
        //            name: String::from("s2"),
        //            data: "VVVVHHHHTTTTTTEEEEEEEEEPPPPPPMMMMMMMMMMLLLLLL".into(),
        //        };
        //        let seq1 = EnumeratedSequence { idx: 0, seq: s1 };
        //        let seq2 = EnumeratedSequence {
        //            idx: 1,
        //            seq: s2.clone(),
        //        };
        //        let seq3 = EnumeratedSequence { idx: 2, seq: s2 };
        //
        //        let mut sequences = vec![seq1, seq2, seq3];

        let (diagonals, part_alig) = align(&sequences);
        let eq_classes = EqClasses::new(&diagonals, &part_alig);
        eq_classes.align_sequences(&mut sequences);
        dbg!(sequences);
        panic!()
    }
}
