use std::collections::HashMap;
use std::error::Error;

use bio::io::fasta::Reader;
use itertools::Itertools;
use rayon::prelude::*;

use fxhash::hash;
use spam_align::align::eq_class::EqClasses;
use spam_align::align::{align, EnumeratedSequence};
use spam_align::data_loaders::balibase::{parse_xml_file, BBAlignment, FilterXmlFile};
use spam_align::score::{score_prot_msa, score_prot_pairwise};
use spam_align::spaced_word::{find_word_match_buckets, Pattern};

fn main() -> Result<(), Box<dyn Error>> {
    let alignment = parse_xml_file("data/bb3_release/RV30/BBS30004.xml").unwrap();
    let mut sequences = alignment
        .unaligned_data
        .into_iter()
        .enumerate()
        .map(|(idx, seq)| EnumeratedSequence::new(idx, seq))
        .collect_vec();
    println!("Running align...");
    let (diagonals, part_alig) = align(&sequences);
    let eq_classes = EqClasses::new(&diagonals, &part_alig);
    eq_classes.align_sequences(&mut sequences);
    dbg!(sequences);
    Ok(())
}
