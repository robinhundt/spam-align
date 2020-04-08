use anyhow::Result;
use itertools::Itertools;
use rayon::prelude::*;
use std::ffi::OsStr;
use std::fs;
use std::fs::{DirEntry, File};
use std::ops::Not;
use std::path::PathBuf;
use std::process::{Command, Stdio};

fn main() -> Result<()> {
    compute_results_for_balibase("data/bb3_release".into(), "mafft_output".into()).unwrap();
    Ok(())
}

fn compute_results_for_balibase(balibase_path: PathBuf, out_path: PathBuf) -> Result<()> {
    let balibase_folders = ["RV11", "RV12", "RV20", "RV30", "RV50"];
    let balibase_folders = balibase_folders.iter().map(|folder| {
        let mut path = balibase_path.clone();
        path.push(folder);
        path
    });
    for folder in balibase_folders {
        let mut out_folder_path = out_path.clone();
        out_folder_path.push(folder.file_name().unwrap());

        fs::create_dir_all(&out_folder_path)?;
        let bb_files = fs::read_dir(folder)?
            .map(Result::unwrap)
            .map(|dir_entry: DirEntry| dir_entry.path())
            .filter(|path| {
                let is_fasta = path.extension() == Some(OsStr::new("tfa"));
                let is_not_bbs = path
                    .file_name()
                    .unwrap()
                    .to_string_lossy()
                    .starts_with("BBS")
                    .not();
                is_fasta && is_not_bbs
            })
            .collect_vec();

        bb_files.par_iter().for_each(|file: &PathBuf| {
            let mut out_path = out_folder_path.clone();
            out_path.push(file.file_name().unwrap());
            let out_file = File::create(out_path).unwrap();
            Command::new("mafft")
                .arg("--quiet")
                .arg("--auto")
                .arg(file)
                .stdout(out_file)
                .spawn()
                .unwrap()
                .wait_with_output();
        })
    }
    Ok(())
}
