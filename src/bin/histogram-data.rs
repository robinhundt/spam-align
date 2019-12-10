use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use bio::io::fasta::Reader;
use itertools::Itertools;
use rayon::prelude::*;
use structopt::StructOpt;

use spam_align::data_loaders::balibase;
use spam_align::data_loaders::balibase::{BBAlignment, FilterXmlFile};
use spam_align::score::{score_prot_msa, score_prot_pairwise};
use spam_align::spaced_word::{find_word_match_buckets, find_word_matches, Pattern};

#[derive(StructOpt, Debug)]
struct Opt {
    #[structopt(name = "FILE", parse(from_os_str))]
    in_file: PathBuf,
    #[structopt(name = "OUT", parse(from_os_str))]
    out_file: PathBuf,
    #[structopt(name = "PATTERN")]
    pattern: Pattern,
}

fn main() -> Result<(), Box<dyn Error>> {
    let opt = Opt::from_args();

    let alignment = balibase::parse_xml_file(opt.in_file).unwrap();
    let pattern = opt.pattern;

    let mut buffered_file = BufWriter::new(File::create(opt.out_file).unwrap());

    //    for spaced_word_match in find_word_matches(&pattern, &alignment.unaligned_data) {
    //        let consistent = spaced_word_match.is_consistent(pattern.len(), &alignment);
    //        let consistent = match consistent {
    //            MatchConsistency::Inconsistent => -1,
    //            MatchConsistency::PartiallyConsistent => 0,
    //            MatchConsistency::Consistent => 1,
    //            MatchConsistency::Unknown => 2,
    //        };
    //        buffered_file.write(format!("{};{}\n", spaced_word_match.score, consistent).as_bytes());
    //    }

    Ok(())
}
