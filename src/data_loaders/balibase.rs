use std::error::Error;
use std::ffi::OsStr;
use std::fs::{self, File};
use std::path::Path;

use anyhow::Result;
use quick_xml::de::from_reader;
use serde::Deserialize;

use FilterXmlFile::*;

use crate::data_loaders::{Alignment, Sequence};
use std::io::BufReader;

#[derive(Debug, Deserialize)]
struct XMLRoot {
    alignment: BBAlignment,
}

#[derive(Debug, Deserialize)]
pub struct BBAlignment {
    #[serde(rename = "aln-name")]
    pub name: String,
    #[serde(rename = "sequence")]
    pub sequences: Vec<BBSequence>,
    #[serde(rename = "column-score")]
    pub core_block: BBCoreBlock,
}

#[derive(Debug, Deserialize)]
pub struct BBSequence {
    #[serde(rename = "seq-name")]
    pub name: String,
    #[serde(rename = "seq-data")]
    pub data: String,
}

#[derive(Debug, Deserialize)]
pub struct BBCoreBlock {
    #[serde(rename = "colsco-data")]
    pub data: String,
}

pub enum FilterXmlFile {
    CoreBlockDataOnly,
    AllData,
}

impl From<BBAlignment> for Alignment {
    fn from(bb_alignment: BBAlignment) -> Self {
        let core_blocks = bb_alignment
            .core_block
            .data
            .split(" ")
            .map(|block| {
                let block_num: i32 = block.trim().parse().unwrap_or_else(|_| {
                    panic!("Encountered non i32 in core block data: {}", block)
                });
                block_num == 1
            })
            .collect();
        let seqs = bb_alignment
            .sequences
            .into_iter()
            .map(|seq| seq.into())
            .collect();

        Self::new(bb_alignment.name, seqs, core_blocks)
    }
}

impl From<BBSequence> for Sequence {
    fn from(bb_seq: BBSequence) -> Self {
        let data = bb_seq
            .data
            .bytes()
            .filter(|el| !el.is_ascii_whitespace())
            .collect();

        Self {
            name: bb_seq.name,
            data,
        }
    }
}

pub fn parse_xml_file(path: impl AsRef<Path>) -> Result<Alignment> {
    let bb_alignment = BBAlignment::from_xml_file(path)?;
    Ok(bb_alignment.into())
}

pub fn parse_xml_files_in_dir(
    path: impl AsRef<Path>,
    file_filter: FilterXmlFile,
) -> Result<Vec<Alignment>> {
    let bb_alignments = BBAlignment::from_xml_files_in_dir(path, file_filter)?;
    Ok(bb_alignments.into_iter().map(|a| a.into()).collect())
}

impl BBAlignment {
    pub fn from_xml_file(path: impl AsRef<Path>) -> Result<BBAlignment> {
        let file = BufReader::new(File::open(&path)?);
        let xml_root: XMLRoot = from_reader(file)?;
        Ok(xml_root.alignment)
    }

    pub fn from_xml_files_in_dir(
        path: impl AsRef<Path>,
        file_filter: FilterXmlFile,
    ) -> Result<Vec<BBAlignment>> {
        fs::read_dir(path)?
            .filter(|dir_entry| match dir_entry {
                Ok(entry) => entry.path().extension() == Some(OsStr::new("xml")),
                Err(_) => false,
            })
            .filter(|dir_entry| {
                let dir_entry = dir_entry.as_ref().unwrap();
                let starts_with_bbs = dir_entry
                    .path()
                    .file_name()
                    .unwrap()
                    .to_string_lossy()
                    .starts_with("BBS");
                match file_filter {
                    CoreBlockDataOnly => starts_with_bbs,
                    AllData => !starts_with_bbs,
                }
            })
            .map(|dir_entry| Self::from_xml_file(dir_entry.unwrap().path()))
            .collect()
    }
}
