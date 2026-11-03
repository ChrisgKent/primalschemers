//! Profile `MisprimeDB::screen`: every primer_size_min (19) window of a
//! reference (both strands) as primers, screened against the reference as an MSA and a background.
//!
//! cargo run --release --example mispriming_profile -- REF.fasta BACKGROUND.fasta [fuzzy]
use primalschemers::mispriming::{for_each_record, MisprimeDB, QueryIndex};
use primalschemers::seqfuncs::reverse_complement;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let (reference, background) = (&args[1], &args[2]);
    let fuzzy = args.get(3).is_none_or(|f| f != "exact");
    // primer_size_min, which is also the kmer size
    let k = 19;

    let mut genome = Vec::new();
    for_each_record(reference, |_, seq| genome = seq).unwrap();
    let mut primers: Vec<Vec<u8>> = genome.windows(k).map(|w| w.to_vec()).collect();
    primers.extend(
        genome
            .windows(k)
            .map(reverse_complement)
            .collect::<Vec<_>>(),
    );
    println!("{} primers, fuzzy={fuzzy}", primers.len());

    let terminals: Vec<Vec<u8>> = primers.iter().map(|p| p[p.len() - k..].to_vec()).collect();
    let t = Instant::now();
    let index = QueryIndex::new(&terminals, k, fuzzy);
    println!("index build: {:.2}s", t.elapsed().as_secs_f64());
    drop(index);

    let t = Instant::now();
    let mut bases = 0usize;
    for_each_record(background, |_, seq| bases += seq.len()).unwrap();
    println!(
        "read only: {:.2}s ({:.2} Gb)",
        t.elapsed().as_secs_f64(),
        bases as f64 / 1e9
    );

    let mut db = MisprimeDB::new(
        vec![reference.clone()],
        vec![background.clone()],
        k,
        fuzzy,
        true,
    )
    .unwrap();
    let t = Instant::now();
    let n = db.screen(&primers).unwrap();
    println!(
        "screen: {:.2}s ({n} terminal kmers)",
        t.elapsed().as_secs_f64()
    );
}
