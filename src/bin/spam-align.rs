#[macro_use]
extern crate anyhow;

use std::ffi::OsStr;
use std::path::PathBuf;

use anyhow::Result;
use spam_align::align::{align, AlignProgress};
use spam_align::spaced_word::read_patterns_from_file;
use spam_align::{read_fasta, write_as_fasta, Sequence};
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
struct Opt {
    #[structopt(short = "I", long, name = "FILE", parse(from_os_str))]
    in_file: PathBuf,
    #[structopt(short = "P", long, name = "PATTERN_SET", parse(from_os_str))]
    pattern_set_path: PathBuf,
    #[structopt(short = "o", long, name = "OUT", parse(from_os_str))]
    out_file: PathBuf,
    /// Show progress information
    #[structopt(short = "p", parse(from_flag))]
    show_progress: AlignProgress,
}

fn main() -> Result<()> {
    let opt: Opt = Opt::from_args();
    let mut sequences = load_sequences(opt.in_file)?;
    let patterns = read_patterns_from_file(opt.pattern_set_path)?;
    align(&mut sequences, &patterns, opt.show_progress);
    write_as_fasta(opt.out_file, &sequences)?;
    Ok(())
}

fn load_sequences(path: PathBuf) -> Result<Vec<Sequence>> {
    match path.extension().and_then(OsStr::to_str) {
        Some("fasta") | Some("tfa") | Some("fa") => read_fasta(path),
        Some(unknown) => Err(anyhow!("Unsupported filename extension: {}", unknown)),
        None => Err(anyhow!("Invalid filename extension")),
    }
}
