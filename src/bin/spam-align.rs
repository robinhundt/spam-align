use itertools::Itertools;
use spam_align::align::align;
use spam_align::align::eq_class::EqClasses;
use spam_align::data_loaders::{balibase, format_as_fasta};
use spam_align::spaced_word::read_patterns_from_file;
use spam_align::Sequences;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let alignment = balibase::parse_xml_file("data/bb3_release/RV30/BBS30004.xml").unwrap();
    let sequences = Sequences::new(&alignment);
    let patterns = read_patterns_from_file("./sample.pat").unwrap();
    let (scored_diagonals, closure) = align(&sequences, &patterns);
    println!("Added {} diagonals", scored_diagonals.len());
    let diagonals = scored_diagonals
        .into_iter()
        .map(|scored| scored.micro_alignment)
        .collect_vec();
    let mut orig_sequences = alignment.unaligned_data;
    let eq_classes = EqClasses::new(&diagonals, &closure);
    eq_classes.align_sequences(&mut orig_sequences);
    println!("{:#?}", diagonals);
    println!("{}", format_as_fasta(&orig_sequences));
    Ok(())
}
