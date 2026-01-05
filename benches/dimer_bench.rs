use criterion::{black_box, criterion_group, criterion_main, Criterion};
use primalschemers::primaldimer;

fn bench_calc_offset(c: &mut Criterion) {
    // Representative primer sequences from test_ensure_consistant_result
    let seq1 = b"ACACCTGTGCCTGTTAAACCAT";
    let seq2 = b"CAATTTGGTAATTGAACACCCATAAAGGT";
    let offset = -12;

    c.bench_function("calc_at_offset", |b| {
        b.iter(|| {
            // Test a specific offset that causes overlap
            primaldimer::calc_at_offset(black_box(seq1), black_box(seq2), black_box(offset))
        })
    });
}

fn bench_do_pool_interact(c: &mut Criterion) {
    let seq1 = b"ACACCTGTGCCTGTTAAACCAT".as_slice();
    let seq2 = b"CAATTTGGTAATTGAACACCCATAAAGGT".as_slice();

    let pool1 = vec![seq1];
    let pool2 = vec![seq2];

    c.bench_function("do_pool_interact_u8_slice", |b| {
        b.iter(|| {
            primaldimer::do_pool_interact_u8_slice(
                black_box(&pool1),
                black_box(&pool2),
                black_box(-9.0),
            )
        })
    });
}

criterion_group!(benches, bench_calc_offset, bench_do_pool_interact);
criterion_main!(benches);
