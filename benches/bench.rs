use anyhow::Result;
use spam_align::align::{align, AlignProgress, Strategy};
use spam_align::read_fasta;
use spam_align::spaced_word::read_patterns_from_file;

use criterion::{criterion_group, criterion_main, Criterion};
use itertools::Itertools;
use spam_align::align::micro_alignment::construct_2dim_micro_alignments;
use spam_align::score::score_prot_msa;

fn bench_align(c: &mut Criterion) {
    let input_12012 = read_fasta("benches/data/BB12012.tfa").expect("fasta file");
    let input_20006 = read_fasta("benches/data/BB20006.tfa").expect("fasta file");
    let input_50001 = read_fasta("benches/data/BB50001.tfa").expect("fasta file");
    let pattern_set = read_patterns_from_file("benches/data/w-3_d-7.out").expect("pattern file");

    let mut group = c.benchmark_group("align");
    group.bench_function("do_align BB12012", |b| {
        b.iter(|| {
            let mut data = input_12012.clone();
            align(
                &mut data,
                &pattern_set,
                Strategy::TwoDim,
                AlignProgress::Hide,
            );
        })
    });
    group.sample_size(10);
    group.bench_function("do_align BB20006", |b| {
        b.iter(|| {
            let mut data = input_20006.clone();
            align(
                &mut data,
                &pattern_set,
                Strategy::TwoDim,
                AlignProgress::Hide,
            );
        })
    });
    group.bench_function("do_align BB50001", |b| {
        b.iter(|| {
            let mut data = input_50001.clone();
            align(
                &mut data,
                &pattern_set,
                Strategy::TwoDim,
                AlignProgress::Hide,
            );
        })
    });

    group.finish()
}

fn bench_find_micro_alignments(c: &mut Criterion) {
    let input = read_fasta("benches/data/BB30022.tfa").expect("fasta file");
    let pattern_set = read_patterns_from_file("benches/data/w-3_d-7.out").expect("pattern file");

    let mut group = c.benchmark_group("find micro_alignments");
    group.sample_size(20);

    group.bench_function("find micro_alignments BB30022", |b| {
        b.iter(|| {
            let _data =
                construct_2dim_micro_alignments(&pattern_set, &input, score_prot_msa).collect_vec();
        })
    });
    group.finish();
}

criterion_group!(benches, bench_align, bench_find_micro_alignments);
criterion_main!(benches);
