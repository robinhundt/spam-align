#[macro_use]
extern crate anyhow;

use std::ffi::OsStr;
use std::path::PathBuf;

use anyhow::Result;
use spam_align::align::{align, AlignProgress, Strategy};
use spam_align::spaced_word::read_patterns_from_file;
use spam_align::{format_as_fasta, read_fasta, write_as_fasta, Sequence};
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
/// SpaM-Align is an experimental multiple alignment
/// tool written by Robin Hundt as part of his
/// bachelor thesis @gobics Göttingen.
struct Opt {
    /// Sequences that should be aligned in fasta format
    #[structopt(short = "I", long, name = "FILE", parse(from_os_str))]
    in_file: PathBuf,
    /// Path to a file containing one pattern per line, consisting of
    /// 1 (match position) and 0 (don't care position). Lines starting
    /// with a '#' are ignored
    #[structopt(short = "P", long, name = "PATTERN_SET", parse(from_os_str))]
    pattern_set_path: PathBuf,
    /// Where to write the resulting alignment in fasta format. If not
    /// provided the result will be printed to stdout
    #[structopt(short = "o", long, name = "OUT", parse(from_os_str))]
    out_file: Option<PathBuf>,
    /// Show progress information
    #[structopt(short = "p", long, parse(from_flag))]
    show_progress: AlignProgress,
    /// Use micro alignments with a dynamic dimension to construct
    /// alignment
    #[structopt(short = "d", long, parse(from_flag))]
    dyn_dim: Strategy,
}

fn main() -> Result<()> {
    env_logger::init();
    let opt: Opt = Opt::from_args();
    let mut sequences = load_sequences(opt.in_file)?;
    let patterns = read_patterns_from_file(opt.pattern_set_path)?;
    align(&mut sequences, &patterns, opt.dyn_dim, opt.show_progress);
    match opt.out_file {
        Some(path) => write_as_fasta(path, &sequences)?,
        None => println!("{}", format_as_fasta(&sequences)),
    }
    Ok(())
}

fn load_sequences(path: PathBuf) -> Result<Vec<Sequence>> {
    match path.extension().and_then(OsStr::to_str) {
        Some("fasta") | Some("tfa") | Some("fa") => read_fasta(path),
        Some(unknown) => Err(anyhow!("Unsupported filename extension: {}", unknown)),
        None => Err(anyhow!("Invalid filename extension")),
    }
}
