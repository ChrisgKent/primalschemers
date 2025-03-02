use digest::{DigestConfig, IndexResult};
use indicatif::ProgressBar;
use pyo3::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{collections::HashMap, time::Duration};

pub mod digest;
pub mod kmer;
pub mod mapping;
pub mod msa;
pub mod primaldimer;
pub mod seqfuncs;
pub mod seqio;
pub mod tm;

#[pyfunction]
fn digest_seq(
    msa_path: &str,
    ncores: usize,
    remap: bool,
    findexes: Option<Vec<usize>>,
    rindexes: Option<Vec<usize>>,
    primer_len_min: Option<usize>,
    primer_len_max: Option<usize>,
    primer_gc_max: Option<f64>,
    primer_gc_min: Option<f64>,
    primer_tm_max: Option<f64>,
    primer_tm_min: Option<f64>,
    max_walk: Option<usize>,
    max_homopolymers: Option<usize>,
    min_freq: Option<f64>,
) -> PyResult<(Vec<kmer::FKmer>, Vec<kmer::RKmer>, Vec<String>)> {
    // Start the spinner
    let spinner = ProgressBar::new_spinner();
    spinner.set_message("Parsing MSA");
    spinner.enable_steady_tick(Duration::from_millis(100));

    // Create config
    let dconf = DigestConfig::new(
        primer_len_min,
        primer_len_max,
        primer_gc_max,
        primer_gc_min,
        primer_tm_max,
        primer_tm_min,
        max_walk,
        max_homopolymers,
        min_freq,
    );

    // Read in the MSA
    let (_headers, seqs) = seqio::fasta_reader(msa_path);
    let mut log_strs: Vec<String> = Vec::new();

    // Create the sequence array
    let seq_array = seqs
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect::<Vec<Vec<u8>>>();

    // Create the threadpool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(ncores)
        .build()
        .unwrap();

    pool.install(|| {
        // Parse input files
        let seq_array = seqio::remove_end_insertions(seq_array);

        // Create the mapping array
        let mapping_array = mapping::create_mapping_array(&seqs[0].as_bytes());

        // Create slices
        let seq_slice = seq_array
            .iter()
            .map(|s| s.as_slice())
            .collect::<Vec<&[u8]>>();

        spinner.finish_and_clear();
        // Create the digest
        let digested_f = digest::digest_f_primer(&seq_slice, &dconf, findexes);
        let digested_r = digest::digest_r_primer(&seq_slice, &dconf, rindexes);

        // Create the reverse digest
        // Start the spinner
        let spinner = ProgressBar::new_spinner();
        spinner.set_message("Processing Kmers");
        spinner.enable_steady_tick(Duration::from_millis(100));

        // Count the errors stats
        let mut fp_count: HashMap<&IndexResult, usize> = HashMap::new();
        for res in digested_f.iter() {
            match res {
                Ok(_) => {
                    let count = fp_count.entry(&IndexResult::Pass()).or_insert(0);
                    *count += 1;
                }
                Err(e) => {
                    let count = fp_count.entry(e).or_insert(0);
                    *count += 1;
                }
            }
        }
        let mut rp_count: HashMap<&IndexResult, usize> = HashMap::new();
        for res in digested_r.iter() {
            match res {
                Ok(_) => {
                    let count = rp_count.entry(&IndexResult::Pass()).or_insert(0);
                    *count += 1;
                }
                Err(e) => {
                    let count = rp_count.entry(e).or_insert(0);
                    *count += 1;
                }
            }
        }

        // Sort values and push to log string vec
        let mut values = fp_count.into_iter().collect::<Vec<(&IndexResult, usize)>>();
        values.sort_by(|a, b| b.1.cmp(&a.1));
        log_strs.push(format!("fprimer status:{:?}", values));
        let mut values = rp_count.into_iter().collect::<Vec<(&IndexResult, usize)>>();
        values.sort_by(|a, b| b.1.cmp(&a.1));
        log_strs.push(format!("fprimer status:{:?}", values));

        let fkmers: Vec<kmer::FKmer> = digested_f.into_par_iter().filter_map(Result::ok).collect();
        let rkmers: Vec<kmer::RKmer> = digested_r.into_par_iter().filter_map(Result::ok).collect();

        // Remap the kmers
        if remap {
            let mut rm_fk: Vec<kmer::FKmer> = Vec::with_capacity(fkmers.len());
            for mut fk in fkmers.into_iter() {
                match mapping_array[fk.end()] {
                    Some(i) => {
                        fk.remap(i);
                        rm_fk.push(fk);
                    }
                    None => {}
                }
            }

            let mut rm_rk: Vec<kmer::RKmer> = Vec::with_capacity(rkmers.len());
            for mut rk in rkmers.into_iter() {
                match mapping_array[rk.start()] {
                    Some(i) => {
                        rk.remap(i);
                        rm_rk.push(rk);
                    }
                    None => {}
                }
            }

            spinner.finish_and_clear();
            return Ok((rm_fk, rm_rk, log_strs));
        } else {
            spinner.finish_and_clear();
            return Ok((fkmers, rkmers, log_strs));
        }
    })
}

// Create mapping array
#[pyfunction]
fn do_seqs_interact(seq1: &[u8], seq2: &[u8], t: f64) -> bool {
    // Create the reverse complement of the sequences
    let mut seq1_rev: Vec<u8> = seq1.to_vec();
    seq1_rev.reverse();
    let mut seq2_rev: Vec<u8> = seq2.to_vec();
    seq2_rev.reverse();
    // Check for interactions
    if primaldimer::does_seq1_extend_no_alloc(&seq1, &seq2_rev, t)
        || primaldimer::does_seq1_extend_no_alloc(&seq2, &seq1_rev, t)
    {
        return true;
    }
    false
}

#[pyfunction]
fn do_pool_interact(seqs1: Vec<Vec<u8>>, seqs2: Vec<Vec<u8>>, t: f64) -> bool {
    // Create the reverse complement of the sequences
    let mut seqs1_rev: Vec<Vec<u8>> = seqs1.iter().map(|s| s.to_vec()).collect();
    for s in seqs1_rev.iter_mut() {
        s.reverse();
    }
    let mut seqs2_rev: Vec<Vec<u8>> = seqs2.iter().map(|s| s.to_vec()).collect();
    for s in seqs2_rev.iter_mut() {
        s.reverse();
    }

    // Check for interactions
    for seq1i in 0..seqs1.len() {
        for seq2i in 0..seqs2.len() {
            if primaldimer::does_seq1_extend_no_alloc(&seqs1[seq1i], &seqs2_rev[seq2i], t)
                || primaldimer::does_seq1_extend_no_alloc(&seqs2[seq2i], &seqs1_rev[seq1i], t)
            {
                return true;
            }
        }
    }
    false
}

/// A Python module implemented in Rust.
#[pyfunction]
fn create_mapping_array(seq: &[u8]) -> PyResult<Vec<Option<usize>>> {
    Ok(mapping::create_mapping_array(seq))
}

#[pymodule]
fn _core(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<kmer::FKmer>()?;
    m.add_class::<kmer::RKmer>()?;
    m.add_function(wrap_pyfunction!(digest_seq, m)?)?;
    m.add_function(wrap_pyfunction!(do_seqs_interact, m)?)?;
    m.add_function(wrap_pyfunction!(do_pool_interact, m)?)?;

    Ok(())
}
