use crate::align::micro_alignment::Site;
use serde::export::Formatter;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::str;
use std::{fmt, io};

pub mod balibase;

#[derive(Clone, Debug)]
pub struct Alignment {
    pub name: String,
    pub aligned_data: Vec<Sequence>,
    pub unaligned_data: Vec<Sequence>,
    unaligned_to_aligned_pos: Vec<Vec<usize>>,
    pub core_blocks: Vec<bool>,
}

#[derive(Clone)]
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

    pub fn pos_aligned(&self, pos1: Site, pos2: Site) -> PositionAlignment {
        let pos_mapping = &self.unaligned_to_aligned_pos;
        let mapped_1 = pos_mapping[pos1.seq][pos1.pos];
        let mapped_2 = pos_mapping[pos2.seq][pos2.pos];
        if mapped_1 == mapped_2 && self.core_blocks[mapped_1] {
            PositionAlignment::Correct
        } else if mapped_1 != mapped_2 {
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
}

impl Sequence {
    pub fn new(name: String, data: Vec<u8>) -> Self {
        Self { name, data }
    }

    /// Returns new Sequence with stripped gaps und Vec that matches each unaligned position to
    /// it's original aligned position
    fn into_unaligned(self) -> (Self, Vec<usize>) {
        let (positions, no_gap_data) = self
            .data
            .into_iter()
            .enumerate()
            .filter(|(_, el)| *el != b'-')
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
    file.write(formatted.as_bytes())?;
    Ok(())
}
