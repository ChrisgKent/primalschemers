use config::{DigestConfig, ThermoType};
use digest::{thermo_check, DigestError, IndexResult, ThermoResult};
use indicatif::ProgressBar;
use pyo3::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{collections::HashMap, time::Duration};

pub mod config;
pub mod digest;
pub mod kmer;
pub mod mapping;
pub mod primaldimer;
pub mod seqfuncs;
pub mod seqio;
pub mod tm;

#[pyclass]
struct Digester {
    remap: bool,
    _seq_array: Vec<Vec<u8>>,
    _mapping_array: Vec<Option<usize>>,
    _ref_to_msa_array: Vec<usize>,
    _thread_pool: rayon::ThreadPool,
}
#[pymethods]
impl Digester {
    #[new]
    #[pyo3(signature = (msa_path, ncores, remap))]
    fn new(msa_path: &str, ncores: usize, remap: bool) -> Self {
        // Build the thread pool
        let _thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(ncores)
            .build()
            .unwrap();

        // Read in the MSA
        let (_headers, seqs) = seqio::fasta_reader(&msa_path);

        // parse the fasta to a seq array
        let mut _seq_array = seqs
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect::<Vec<Vec<u8>>>();
        _seq_array = seqio::remove_end_insertions(_seq_array);

        // Create the mapping array
        let _mapping_array: Vec<Option<usize>> =
            mapping::create_mapping_array(&seqs[0].as_bytes(), remap);
        let _ref_to_msa_array = mapping::create_ref_to_msa(&_mapping_array);

        // Return the Digester
        Digester {
            remap,
            _seq_array,
            _thread_pool,
            _mapping_array,
            _ref_to_msa_array,
        }
    }
    fn create_seq_slice(&self) -> Vec<&[u8]> {
        self._seq_array
            .iter()
            .map(|s| s.as_slice())
            .collect::<Vec<&[u8]>>()
    }

    #[pyo3(signature = (findexes=None, rindexes=None, primer_len_min=None, primer_len_max=None, primer_gc_max=None, primer_gc_min=None, primer_tm_max=None, primer_tm_min=None, primer_annealing_prop=None, annealing_temp_c=None, max_walk=None, max_homopolymers=None, min_freq=None, ignore_n=None, dimerscore=None, return_counts=false ))]
    pub fn digest(
        &self,
        findexes: Option<Vec<usize>>,
        rindexes: Option<Vec<usize>>,
        primer_len_min: Option<usize>,
        primer_len_max: Option<usize>,
        primer_gc_max: Option<f64>,
        primer_gc_min: Option<f64>,
        primer_tm_max: Option<f64>,
        primer_tm_min: Option<f64>,
        primer_annealing_prop: Option<f64>,
        annealing_temp_c: Option<f64>,
        max_walk: Option<usize>,
        max_homopolymers: Option<usize>,
        min_freq: Option<f64>,
        ignore_n: Option<bool>,
        dimerscore: Option<f64>,
        return_counts: bool,
    ) -> PyResult<(Vec<kmer::FKmer>, Vec<kmer::RKmer>, Vec<String>)> {
        // If both annealing are set use annealing
        let thermo_type = match (primer_annealing_prop, annealing_temp_c) {
            (Some(_), Some(_)) => ThermoType::ANNEALING,
            _ => ThermoType::TM,
        };

        let dconf = DigestConfig::new(
            primer_len_min,
            primer_len_max,
            primer_gc_max,
            primer_gc_min,
            primer_tm_max,
            primer_tm_min,
            primer_annealing_prop,
            annealing_temp_c,
            Some(thermo_type),
            max_walk,
            max_homopolymers,
            min_freq,
            ignore_n,
            dimerscore,
        );

        let seq_slice = self.create_seq_slice();

        let mut log_strs = Vec::new();

        self._thread_pool.install(|| {
            // Create the digest
            let digested_f = digest::digest_f_primer(&seq_slice, &dconf, findexes);
            let digested_r = digest::digest_r_primer(&seq_slice, &dconf, rindexes);

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
            log_strs.push(format!("rprimer status:{:?}", values));

            let fkmers: Vec<kmer::FKmer> =
                digested_f.into_par_iter().filter_map(Result::ok).collect();
            let rkmers: Vec<kmer::RKmer> =
                digested_r.into_par_iter().filter_map(Result::ok).collect();

            // Remap the kmers
            if self.remap {
                let mut rm_fk: Vec<kmer::FKmer> = Vec::with_capacity(fkmers.len());
                for mut fk in fkmers.into_iter() {
                    // Check for being on last base
                    if fk.end() == self._mapping_array.len() {
                        match self._mapping_array[fk.end() - 1] {
                            Some(i) => {
                                fk.remap(i + 1);
                                rm_fk.push(fk);
                            }
                            None => {}
                        }
                    } else {
                        match self._mapping_array[fk.end()] {
                            Some(i) => {
                                fk.remap(i);
                                rm_fk.push(fk);
                            }
                            None => {}
                        }
                    }
                }

                let mut rm_rk: Vec<kmer::RKmer> = Vec::with_capacity(rkmers.len());
                for mut rk in rkmers.into_iter() {
                    match self._mapping_array[rk.start()] {
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
    #[getter]
    fn _seq_array(&self) -> &Vec<Vec<u8>> {
        &self._seq_array
    }
    #[getter]
    fn _mapping_array(&self) -> &Vec<Option<usize>> {
        &self._mapping_array
    }
    #[getter]
    fn _ref_to_msa_array(&self) -> &Vec<usize> {
        &self._ref_to_msa_array
    }
}

#[pyfunction]
#[pyo3(signature = (msa_path, ncores, remap, findexes=None, rindexes=None, primer_len_min=None, primer_len_max=None, primer_gc_max=None, primer_gc_min=None, primer_tm_max=None, primer_tm_min=None, primer_annealing_prop=None, annealing_temp_c=None, max_walk=None, max_homopolymers=None, min_freq=None, ignore_n=None, dimerscore=None))]
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
    // tm
    primer_tm_max: Option<f64>,
    primer_tm_min: Option<f64>,
    // annealing
    primer_annealing_prop: Option<f64>,
    annealing_temp_c: Option<f64>,

    max_walk: Option<usize>,
    max_homopolymers: Option<usize>,
    min_freq: Option<f64>,
    ignore_n: Option<bool>,
    dimerscore: Option<f64>,
) -> PyResult<(Vec<kmer::FKmer>, Vec<kmer::RKmer>, Vec<String>)> {
    // Start the spinner
    let spinner = ProgressBar::new_spinner();
    spinner.set_message("Parsing MSA");
    spinner.enable_steady_tick(Duration::from_millis(100));

    // If both annealing are set use annealing
    let thermo_type = match (primer_annealing_prop, annealing_temp_c) {
        (Some(_), Some(_)) => ThermoType::ANNEALING,
        _ => ThermoType::TM,
    };

    // Create config
    let dconf = DigestConfig::new(
        primer_len_min,
        primer_len_max,
        primer_gc_max,
        primer_gc_min,
        primer_tm_max,
        primer_tm_min,
        primer_annealing_prop,
        annealing_temp_c,
        Some(thermo_type),
        max_walk,
        max_homopolymers,
        min_freq,
        ignore_n,
        dimerscore,
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
        let mapping_array = mapping::create_mapping_array(&seqs[0].as_bytes(), remap);

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
        log_strs.push(format!("rprimer status:{:?}", values));

        let fkmers: Vec<kmer::FKmer> = digested_f.into_par_iter().filter_map(Result::ok).collect();
        let rkmers: Vec<kmer::RKmer> = digested_r.into_par_iter().filter_map(Result::ok).collect();

        // Remap the kmers
        if remap {
            let mut rm_fk: Vec<kmer::FKmer> = Vec::with_capacity(fkmers.len());
            for mut fk in fkmers.into_iter() {
                // Check for being on last base
                if fk.end() == mapping_array.len() {
                    match mapping_array[fk.end() - 1] {
                        Some(i) => {
                            fk.remap(i + 1);
                            rm_fk.push(fk);
                        }
                        None => {}
                    }
                } else {
                    match mapping_array[fk.end()] {
                        Some(i) => {
                            fk.remap(i);
                            rm_fk.push(fk);
                        }
                        None => {}
                    }
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

#[pyfunction]
fn calc_at_offset_py(seq1: &str, seq2: &str, offset: i32) -> f64 {
    //Provide strings in 5'-3'
    // This will return the score for this offset
    let seq1: Vec<u8> = seq1.as_bytes().to_vec();
    let mut seq2 = seq2.as_bytes().to_vec();
    seq2.reverse();

    match primaldimer::calc_at_offset(&seq1, &seq2, offset) {
        Some(score) => return score,
        None => return 100.,
    };
}

#[pymodule]
fn _core(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Object classes
    m.add_class::<kmer::FKmer>()?;
    m.add_class::<kmer::RKmer>()?;
    m.add_class::<Digester>()?;

    // Functions
    m.add_function(wrap_pyfunction!(digest_seq, m)?)?;
    m.add_function(wrap_pyfunction!(calc_at_offset_py, m)?)?;
    m.add_function(wrap_pyfunction!(do_seqs_interact, m)?)?;
    m.add_function(wrap_pyfunction!(do_pool_interact, m)?)?;
    m.add_function(wrap_pyfunction!(kmer::generate_primerpairs_py, m)?)?;

    Ok(())
}
