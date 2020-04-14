use anyhow::Result;
use criterion::{criterion_group, criterion_main, Criterion};
use spam_align::align::{align, AlignProgress};
use spam_align::read_fasta;
use spam_align::spaced_word::read_patterns_from_file;

fn bench_align(c: &mut Criterion) {
    let input_12012 = read_fasta("benches/data/BB12012.tfa").expect("fasta file");
    let input_20006 = read_fasta("benches/data/BB20006.tfa").expect("fasta file");
    let input_50001 = read_fasta("benches/data/BB50001.tfa").expect("fasta file");
    let pattern_set = read_patterns_from_file("benches/data/w-3_d-7.out").expect("pattern file");

    let mut group = c.benchmark_group("align");
    group.bench_function("align BB12012", |b| {
        b.iter(|| {
            let mut data = input_12012.clone();
            align(&mut data, &pattern_set, AlignProgress::Hide);
        })
    });
    group.sample_size(10);
    group.bench_function("align BB20006", |b| {
        b.iter(|| {
            let mut data = input_20006.clone();
            align(&mut data, &pattern_set, AlignProgress::Hide);
        })
    });
    group.bench_function("align BB50001", |b| {
        b.iter(|| {
            let mut data = input_50001.clone();
            align(&mut data, &pattern_set, AlignProgress::Hide);
        })
    });

    group.finish()
}

// fn criterion_benchmark(c: &mut Criterion) {
//     c.bench_function("fib 20", |b| b.iter(|| fibonacci(black_box(20))));
// }

criterion_group!(benches, bench_align);
criterion_main!(benches);
