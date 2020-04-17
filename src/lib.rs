#[macro_use]
extern crate anyhow;
use std::fmt::Formatter;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::str;
use std::{fmt, io};

use anyhow::{Context, Result};
use bio::io::fasta;

use crate::align::micro_alignment::Site;
use std::time::Instant;

pub mod align;
pub mod score;
pub mod spaced_word;

#[derive(Clone, PartialEq)]
pub struct Sequence {
    pub name: String,
    pub data: Vec<u8>,
}

impl Sequence {
    pub fn len(&self) -> usize {
        self.data.len()
    }
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}

impl fmt::Debug for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Sequence {{ name: {}, data: {} }}",
            self.name,
            str::from_utf8(&self.data).unwrap()
        )
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, ">{}\n{}", self.name, str::from_utf8(&self.data).unwrap())
    }
}

pub fn format_as_fasta(seqs: &[Sequence]) -> String {
    use std::fmt::Write;
    let mut buf = String::new();
    for seq in seqs {
        writeln!(buf, "{}", seq).expect("Write fasta buffer")
    }
    buf
}

pub fn write_as_fasta(path: impl AsRef<Path>, seqs: &[Sequence]) -> io::Result<()> {
    let mut file = File::create(path)?;
    let formatted = format_as_fasta(seqs);
    file.write_all(formatted.as_bytes())?;
    Ok(())
}

pub fn read_fasta(path: impl AsRef<Path>) -> Result<Vec<Sequence>> {
    let file = File::open(&path)?;
    let reader = fasta::Reader::new(file);
    let sequences: Result<Vec<_>> = reader
        .records()
        .map(|record| match record {
            Ok(record) => Ok(Sequence {
                name: record.id().to_string(),
                data: record.seq().to_vec(),
            }),
            Err(err) => Err(err.into()),
        })
        .collect();
    sequences.context("Reading record")
}

fn log_elapsed_time<T, F: FnOnce() -> T>(name: &str, f: F) -> T {
    let now = Instant::now();
    let ret = f();
    log::info!("<{}>: {}ms", name, now.elapsed().as_millis());
    ret
}
