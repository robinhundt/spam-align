#[macro_use]
extern crate anyhow;
use anyhow::Result;

use spam_align::align::eq_class::EqClasses;

use spam_align::align::{align, AlignProgress};
use spam_align::spaced_word::read_patterns_from_file;
use spam_align::{write_as_fasta, Alignment};

use std::ffi::OsStr;

use std::path::PathBuf;

use structopt::StructOpt;

#[derive(StructOpt, Debug)]
struct Opt {
    #[structopt(short = "I", long, name = "FILE", parse(from_os_str))]
    in_file: PathBuf,
    #[structopt(short = "P", long, name = "PATTERN_SET", parse(from_os_str))]
    pattern_set_path: PathBuf,
    #[structopt(short = "o", long, name = "OUT", parse(from_os_str))]
    out_file: PathBuf,
    #[structopt(short = "p", parse(from_flag))]
    show_progrss: AlignProgress,
}

fn main() -> Result<()> {
    let opt: Opt = Opt::from_args();
    let alignment = load_alignment(opt.in_file)?;
    let patterns = read_patterns_from_file(opt.pattern_set_path)?;
    let (_, closure) = align(&alignment.unaligned_data, &patterns, AlignProgress::Show);
    let eq_classes = EqClasses::new(&closure);
    let mut alignment = alignment;
    eq_classes.align_sequences(&mut alignment.unaligned_data);
    write_as_fasta(opt.out_file, &alignment.unaligned_data)?;
    Ok(())
}

fn load_alignment(path: PathBuf) -> Result<Alignment> {
    let alignment = match path.extension().and_then(OsStr::to_str) {
        Some("fasta") | Some("tfa") | Some("fa") => Alignment::read_fasta(path)?,
        Some(unknown) => return Err(anyhow!("Unsupported filename extension: {}", unknown)),
        None => return Err(anyhow!("Invalid filename extension")),
    };
    Ok(alignment)
}
