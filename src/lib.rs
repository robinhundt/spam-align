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

pub mod align;
pub mod score;
pub mod spaced_word;

#[derive(Clone, Debug)]
pub struct Alignment {
    pub name: String,
    pub aligned_data: Vec<Sequence>,
    pub unaligned_data: Vec<Sequence>,
    unaligned_to_aligned_pos: Vec<Vec<usize>>,
    pub core_blocks: Vec<bool>,
}

#[derive(Clone, PartialEq)]
pub struct Sequence {
    pub name: String,
    pub data: Vec<u8>,
}

#[derive(Eq, PartialEq, Copy, Clone, Debug, Hash)]
pub enum PositionAlignment {
    Correct,
    Incorrect,
    Unknown,
}

#[derive(Eq, PartialEq, Hash, Debug)]
pub enum MicroAlignmentCheck {
    Correct,
    PartiallyCorrectIncorrect,
    PartiallyCorrectUnknown,
    Incorrect,
    Unknown,
}

impl Alignment {
    pub const GAP_CHARACTER: u8 = b'-';
    pub fn new(name: String, aligned_data: Vec<Sequence>, core_blocks: Vec<bool>) -> Self {
        let (unaligned_data, pos_mapping) = aligned_data
            .iter()
            .map(|seq| seq.clone().into_unaligned())
            .unzip();

        Self {
            name,
            aligned_data,
            unaligned_data,
            unaligned_to_aligned_pos: pos_mapping,
            core_blocks,
        }
    }

    pub fn read_fasta(path: impl AsRef<Path>) -> Result<Self> {
        let file = File::open(&path)?;
        let reader = fasta::Reader::new(file);
        let aligned_data: Result<Vec<_>> = reader
            .records()
            .map(|record| {
                match record {
                    Ok(record) => Ok(Sequence::new(
                        record.id().to_string(),
                        record.seq().to_vec(),
                    )),
                    Err(err) => Err(err.into()),
                }
                // record.map(|record| Sequence::new(record.id().to_string(), record.seq().to_vec()))
            })
            .collect();
        let aligned_data = aligned_data.context("Reading record")?;
        let name = path
            .as_ref()
            .file_stem()
            .context("No filestem for fasta file")?
            .to_string_lossy()
            .into_owned();
        let alig_len = aligned_data.first().context("Empty alignment")?.data.len();
        let core_blocks = vec![true; alig_len];

        Ok(Self::new(name, aligned_data, core_blocks))
    }

    pub fn pos_aligned(&self, pos1: Site, pos2: Site) -> PositionAlignment {
        let pos_mapping = &self.unaligned_to_aligned_pos;
        let mapped_1 = pos_mapping[pos1.seq][pos1.pos];
        let mapped_2 = pos_mapping[pos2.seq][pos2.pos];
        if mapped_1 == mapped_2 && self.core_blocks[mapped_1] {
            PositionAlignment::Correct
        } else if mapped_1 != mapped_2 && (self.core_blocks[mapped_1] || self.core_blocks[mapped_2])
        {
            PositionAlignment::Incorrect
        } else {
            // mapped to position is not part of core block region
            PositionAlignment::Unknown
        }
    }

    pub fn core_block_data(&self) -> Vec<Sequence> {
        self.aligned_data
            .iter()
            .map(|seq| seq.select_core_blocks(&self.core_blocks))
            .collect()
    }

    pub fn core_block_columns(&self) -> usize {
        self.core_blocks.iter().copied().map(usize::from).sum()
    }
}

impl Sequence {
    pub fn new(name: String, data: Vec<u8>) -> Self {
        Self { name, data }
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Returns new Sequence with stripped gaps und Vec that matches each unaligned position to
    /// it's original aligned position
    fn into_unaligned(self) -> (Self, Vec<usize>) {
        let (positions, no_gap_data) = self
            .data
            .into_iter()
            .enumerate()
            .filter(|(_, el)| *el != Alignment::GAP_CHARACTER)
            .unzip();
        let unaligned_sequence = Self {
            name: self.name,
            data: no_gap_data,
        };
        (unaligned_sequence, positions)
    }

    fn select_core_blocks(&self, core_blocks: &[bool]) -> Self {
        if self.data.len() != core_blocks.len() {
            panic!("Sequence len must be equal to core blocks len")
        }
        let data = self
            .data
            .iter()
            .zip(core_blocks)
            .filter_map(|(el, core_block)| if *core_block { Some(*el) } else { None })
            .collect();
        Self {
            name: self.name.clone(),
            data,
        }
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
        .map(|record| {
            match record {
                Ok(record) => Ok(Sequence::new(
                    record.id().to_string(),
                    record.seq().to_vec(),
                )),
                Err(err) => Err(err.into()),
            }
            // record.map(|record| Sequence::new(record.id().to_string(), record.seq().to_vec()))
        })
        .collect();
    sequences.context("Reading record")
}
