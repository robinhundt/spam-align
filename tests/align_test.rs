use std::fs;

use anyhow::{Context, Result};
use spam_align::align::{align, AlignProgress};
use spam_align::read_fasta;
use spam_align::spaced_word::read_patterns_from_file;

#[test]
fn compare_againts_previous_alignments() -> Result<()> {
    // test against precomputed alignemts in "tests/fixtures/alignments" in order
    // to ensure that refactoring/optimizations do not change the output

    let expected_alignments: Result<Vec<_>> = fs::read_dir("tests/fixtures/alignments")
        .context("Reading alignments")?
        .map(|dir_entry| {
            let dir_entry = dbg!(dir_entry.unwrap());
            read_fasta(dir_entry.path())
        })
        .collect();
    let mut expected_alignments = expected_alignments?;

    let input_data: Result<Vec<_>> = fs::read_dir("tests/fixtures/unaligned")
        .context("Reading unaligned data")?
        .map(|dir_entry| {
            let dir_entry = dbg!(dir_entry.unwrap());
            read_fasta(dir_entry.path())
        })
        .collect();
    let mut input_data = input_data?;

    let pattern_set =
        read_patterns_from_file("tests/fixtures/w-3_d-7.out").context("Reading patterns")?;

    for (sequences, alignment) in input_data.iter_mut().zip(expected_alignments) {
        align(sequences, &pattern_set, AlignProgress::Hide);
        assert_eq!(&alignment, sequences)
    }

    Ok(())
}
