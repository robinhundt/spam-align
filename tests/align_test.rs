use std::fs;

use anyhow::{Context, Result};
use spam_align::align::{align, AlignProgress};
use spam_align::spaced_word::read_patterns_from_file;
use spam_align::Alignment;

#[test]
fn compare_againts_previous_alignments() -> Result<()> {
    // test against precomputed alignemts in "tests/fixtures/alignments" in order
    // to ensure that refactoring/optimizations do not change the output

    let alignments: Result<Vec<_>> = fs::read_dir("tests/fixtures/alignments")
        .context("Reading alignments")?
        .map(|dir_entry| {
            let dir_entry = dir_entry.unwrap();
            Alignment::read_fasta(dir_entry.path())
        })
        .collect();
    let mut alignments = alignments?;

    let pattern_set =
        read_patterns_from_file("tests/fixtures/w-3_d-7.out").context("Reading patterns")?;

    for alignment in alignments.iter_mut() {
        align(
            &mut alignment.unaligned_data,
            &pattern_set,
            AlignProgress::Hide,
        );
        assert_eq!(alignment.aligned_data, alignment.unaligned_data)
    }

    Ok(())
}
