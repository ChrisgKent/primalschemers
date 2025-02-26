use pyo3::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

pub mod digest;
pub mod kmer;
pub mod primaldimer;
pub mod seqfuncs;
pub mod seqio;
pub mod tm;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

#[pyfunction]
fn digest_seq(msa_path: &str, ncores: usize) -> PyResult<(Vec<kmer::FKmer>, Vec<kmer::RKmer>)> {
    // Read in the MSA
    let (headers, seqs) = seqio::fasta_reader(msa_path);

    let seq_array = seqs.iter().map(|s| s.as_bytes()).collect::<Vec<&[u8]>>();
    // let seqs = seqio::remove_end_insertions(seq_array);

    // Create the threadpool

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(ncores)
        .build()
        .unwrap();

    pool.install(|| {
        // Create the digest
        let digested_f = digest::digest_f_primer(&seq_array);
        let fkmers: Vec<kmer::FKmer> = digested_f.into_iter().filter_map(Result::ok).collect();

        // Create the reverse digest
        let digested_r = digest::digest_r_primer(&seq_array);
        let rkmers: Vec<kmer::RKmer> = digested_r.into_par_iter().filter_map(Result::ok).collect();
        Ok((fkmers, rkmers))
    })
}

/// A Python module implemented in Rust.

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_class::<kmer::FKmer>()?;
    m.add_class::<kmer::RKmer>()?;
    m.add_function(wrap_pyfunction!(digest_seq, m)?)?;
    Ok(())
}
