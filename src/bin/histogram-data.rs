use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use structopt::StructOpt;

use spam_align::data_loaders::balibase;
//use spam_align::data_loaders::balibase::{BBAlignment, FilterXmlFile};
use spam_align::spaced_word::Pattern;

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

    let _alignment = balibase::parse_xml_file(opt.in_file).unwrap();
    let _pattern = opt.pattern;

    let mut _buffered_file = BufWriter::new(File::create(opt.out_file).unwrap());

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
