use config::{DigestConfig, ThermoType};
use digest::IndexResult;
use indicatif::ProgressBar;
use pyo3::exceptions::PyValueError;
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
pub mod thal;
pub mod tm;

/// An MSA's sequences and its MSA-to-reference mapping.
type Msa = (Vec<Vec<u8>>, Vec<Option<usize>>);

/// Read an MSA for digestion: its sequences (end insertions removed) and the
/// MSA-to-reference mapping. ValueError for an unreadable or empty file, or
/// sequences of different lengths.
fn load_msa(msa_path: &str, remap: bool) -> PyResult<Msa> {
    let (_headers, seqs) = seqio::try_fasta_reader(msa_path).map_err(PyValueError::new_err)?;
    let Some(first) = seqs.first() else {
        return Err(PyValueError::new_err(format!(
            "{msa_path}: no sequences in the MSA"
        )));
    };
    if let Some(other) = seqs.iter().find(|s| s.len() != first.len()) {
        return Err(PyValueError::new_err(format!(
            "{msa_path}: sequences are not the same length ({} and {})",
            first.len(),
            other.len()
        )));
    }
    let mapping_array = mapping::create_mapping_array(first.as_bytes(), remap);
    let seq_array = seqs.into_iter().map(String::into_bytes).collect();
    Ok((seqio::remove_end_insertions(seq_array), mapping_array))
}

/// Reaction conditions for primer Tm and annealing: the digest's defaults
/// ([`config::conditions`]) with any given value replaced.
fn digest_conditions(
    mv_conc: Option<f64>,
    dv_conc: Option<f64>,
    dntp_conc: Option<f64>,
    dna_conc: Option<f64>,
    formamide_conc: Option<f64>,
) -> primer3_thal_rs::ThalArgs {
    let d = config::conditions();
    primer3_thal_rs::ThalArgs {
        mv: mv_conc.unwrap_or(d.mv),
        dv: dv_conc.unwrap_or(d.dv),
        dntp: dntp_conc.unwrap_or(d.dntp),
        dna_conc: dna_conc.unwrap_or(d.dna_conc),
        formamide_conc: formamide_conc.unwrap_or(d.formamide_conc),
        ..d
    }
}

/// Count each digest outcome, most common first, for the log.
fn status_counts<K: KmerPassed>(results: &[Result<K, IndexResult>]) -> Vec<(&IndexResult, usize)> {
    let mut counts: HashMap<&IndexResult, usize> = HashMap::new();
    for res in results {
        let status = match res {
            Ok(kmer) if kmer.passed() => &digest::PASSED,
            Ok(_) => &digest::RETURNED_FAILING,
            Err(e) => e,
        };
        *counts.entry(status).or_insert(0) += 1;
    }
    let mut values: Vec<(&IndexResult, usize)> = counts.into_iter().collect();
    values.sort_by(|a, b| b.1.cmp(&a.1));
    values
}

trait KmerPassed {
    fn passed(&self) -> bool;
}
impl KmerPassed for kmer::FKmer {
    fn passed(&self) -> bool {
        kmer::FKmer::passed(self)
    }
}
impl KmerPassed for kmer::RKmer {
    fn passed(&self) -> bool {
        kmer::RKmer::passed(self)
    }
}

/// Digest every requested index, then remap the kmers to reference
/// coordinates if asked, dropping those on positions absent from it.
fn run_digest(
    seqs: &[Vec<u8>],
    mapping_array: &[Option<usize>],
    remap: bool,
    dconf: &DigestConfig,
    findexes: Option<Vec<usize>>,
    rindexes: Option<Vec<usize>>,
) -> (Vec<kmer::FKmer>, Vec<kmer::RKmer>, Vec<String>) {
    let seq_slice: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
    let digested_f = digest::digest_f_primer(&seq_slice, dconf, findexes);
    let digested_r = digest::digest_r_primer(&seq_slice, dconf, rindexes);

    let spinner = ProgressBar::new_spinner();
    spinner.set_message("Processing Kmers");
    spinner.enable_steady_tick(Duration::from_millis(100));

    let log_strs = vec![
        format!("fprimer status:{:?}", status_counts(&digested_f)),
        format!("rprimer status:{:?}", status_counts(&digested_r)),
    ];

    let fkmers = digested_f.into_par_iter().filter_map(Result::ok);
    let rkmers = digested_r.into_par_iter().filter_map(Result::ok);
    let result = if remap {
        let fkmers = fkmers
            .filter_map(|mut fk| {
                // On the last base the end maps one past the last position
                let mapped = if fk.end() == mapping_array.len() {
                    mapping_array[fk.end() - 1].map(|i| i + 1)
                } else {
                    mapping_array[fk.end()]
                };
                mapped.map(|i| {
                    fk.remap(i);
                    fk
                })
            })
            .collect();
        let rkmers = rkmers
            .filter_map(|mut rk| {
                mapping_array[rk.start()].map(|i| {
                    rk.remap(i);
                    rk
                })
            })
            .collect();
        (fkmers, rkmers, log_strs)
    } else {
        (fkmers.collect(), rkmers.collect(), log_strs)
    };
    spinner.finish_and_clear();
    result
}

/// Annealing mode when a target annealing proportion is given (at
/// `annealing_temp_c`, or the config's default temperature), else Tm mode.
fn thermo_type(primer_annealing_prop: Option<f64>) -> ThermoType {
    if primer_annealing_prop.is_some() {
        ThermoType::ANNEALING
    } else {
        ThermoType::TM
    }
}

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
    /// ValueError for an unreadable or empty MSA, or sequences of different
    /// lengths.
    #[new]
    #[pyo3(signature = (msa_path, ncores, remap))]
    fn new(msa_path: &str, ncores: usize, remap: bool) -> PyResult<Self> {
        let _thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(ncores)
            .build()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let (_seq_array, _mapping_array) = load_msa(msa_path, remap)?;
        let _ref_to_msa_array: Vec<usize> = mapping::create_ref_to_msa(&_mapping_array);
        Ok(Digester {
            remap,
            _seq_array,
            _thread_pool,
            _mapping_array,
            _ref_to_msa_array,
        })
    }
    fn create_seq_slice(&self) -> Vec<&[u8]> {
        self._seq_array
            .iter()
            .map(|s| s.as_slice())
            .collect::<Vec<&[u8]>>()
    }

    /// Digest the MSA. Runs without the GIL, so other Python threads keep
    /// going. Concentrations (mM, DNA nM) and formamide (mol/l) set the
    /// Tm/annealing conditions, defaulting to primalscheme3's; the hairpin
    /// check uses the same salts without formamide, as primer3.
    #[pyo3(signature = (findexes=None, rindexes=None, primer_len_min=None, primer_len_max=None, primer_gc_max=None, primer_gc_min=None, primer_tm_max=None, primer_tm_min=None, primer_annealing_prop=None, primer_annealing_diff=None, annealing_temp_c=None, max_walk=None, max_homopolymers=None, min_freq=None, ignore_n=None, dimerscore=None, thermo_check=false, primer_hairpin_th_max=None, return_all=false, mv_conc=None, dv_conc=None, dntp_conc=None, dna_conc=None, formamide_conc=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn digest(
        &self,
        py: Python<'_>,
        findexes: Option<Vec<usize>>,
        rindexes: Option<Vec<usize>>,
        primer_len_min: Option<usize>,
        primer_len_max: Option<usize>,
        primer_gc_max: Option<f64>,
        primer_gc_min: Option<f64>,
        primer_tm_max: Option<f64>,
        primer_tm_min: Option<f64>,
        primer_annealing_prop: Option<f64>,
        primer_annealing_diff: Option<f64>,
        annealing_temp_c: Option<f64>,
        max_walk: Option<usize>,
        max_homopolymers: Option<usize>,
        min_freq: Option<f64>,
        ignore_n: Option<bool>,
        dimerscore: Option<f64>,
        thermo_check: bool,
        primer_hairpin_th_max: Option<f64>,
        return_all: bool,
        mv_conc: Option<f64>,
        dv_conc: Option<f64>,
        dntp_conc: Option<f64>,
        dna_conc: Option<f64>,
        formamide_conc: Option<f64>,
    ) -> PyResult<(Vec<kmer::FKmer>, Vec<kmer::RKmer>, Vec<String>)> {
        let dconf = DigestConfig::new(
            primer_len_min,
            primer_len_max,
            primer_gc_max,
            primer_gc_min,
            primer_tm_max,
            primer_tm_min,
            primer_annealing_prop,
            annealing_temp_c,
            primer_annealing_diff,
            Some(thermo_type(primer_annealing_prop)),
            max_walk,
            max_homopolymers,
            min_freq,
            ignore_n,
            dimerscore,
            Some(thermo_check),
            primer_hairpin_th_max,
            Some(return_all),
            Some(digest_conditions(
                mv_conc,
                dv_conc,
                dntp_conc,
                dna_conc,
                formamide_conc,
            )),
        );
        Ok(py.detach(|| {
            self._thread_pool.install(|| {
                run_digest(
                    &self._seq_array,
                    &self._mapping_array,
                    self.remap,
                    &dconf,
                    findexes,
                    rindexes,
                )
            })
        }))
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

/// Read and digest an MSA in one call: `Digester(msa_path, ncores,
/// remap).digest(...)`, except that `thermo_check` defaults to on.
#[pyfunction]
#[pyo3(signature = (msa_path,ncores, remap, findexes=None, rindexes=None, primer_len_min=None, primer_len_max=None, primer_gc_max=None, primer_gc_min=None, primer_tm_max=None, primer_tm_min=None, primer_annealing_prop=None, primer_annealing_diff=None,annealing_temp_c=None, max_walk=None, max_homopolymers=None, min_freq=None, ignore_n=None, dimerscore=None, thermo_check=true, primer_hairpin_th_max=None, return_all=false, mv_conc=None, dv_conc=None, dntp_conc=None, dna_conc=None, formamide_conc=None))]
#[allow(clippy::too_many_arguments)]
fn digest_seq(
    py: Python<'_>,
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
    primer_annealing_prop: Option<f64>,
    primer_annealing_diff: Option<f64>,
    annealing_temp_c: Option<f64>,
    max_walk: Option<usize>,
    max_homopolymers: Option<usize>,
    min_freq: Option<f64>,
    ignore_n: Option<bool>,
    dimerscore: Option<f64>,
    thermo_check: bool,
    primer_hairpin_th_max: Option<f64>,
    return_all: bool,
    mv_conc: Option<f64>,
    dv_conc: Option<f64>,
    dntp_conc: Option<f64>,
    dna_conc: Option<f64>,
    formamide_conc: Option<f64>,
) -> PyResult<(Vec<kmer::FKmer>, Vec<kmer::RKmer>, Vec<String>)> {
    let digester = Digester::new(msa_path, ncores, remap)?;
    digester.digest(
        py,
        findexes,
        rindexes,
        primer_len_min,
        primer_len_max,
        primer_gc_max,
        primer_gc_min,
        primer_tm_max,
        primer_tm_min,
        primer_annealing_prop,
        primer_annealing_diff,
        annealing_temp_c,
        max_walk,
        max_homopolymers,
        min_freq,
        ignore_n,
        dimerscore,
        thermo_check,
        primer_hairpin_th_max,
        return_all,
        mv_conc,
        dv_conc,
        dntp_conc,
        dna_conc,
        formamide_conc,
    )
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
    if primaldimer::does_seq1_extend_no_alloc(seq1, &seq2_rev, t)
        || primaldimer::does_seq1_extend_no_alloc(seq2, &seq1_rev, t)
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

    primaldimer::calc_at_offset(&seq1, &seq2, offset).unwrap_or(100.)
}

/// Percentage of `sequence` bound to its complement at `annealing_temp_c`,
/// on the same model as the digest and `calc_tm`. DMSO and formamide don't
/// affect binding and are accepted for compatibility. Raises ValueError for
/// an empty sequence or a base other than A/C/G/T.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn calc_annealing(
    sequence: &str,
    dna_nm: f64,
    k_mm: f64,
    divalent_conc: f64,
    dntp_conc: f64,
    dmso_conc: f64,
    dmso_fact: f64,
    formamide_conc: f64,
    annealing_temp_c: f64,
) -> PyResult<f64> {
    let args = primer3_thal_rs::ThalArgs {
        mv: k_mm,
        dv: divalent_conc,
        dntp: dntp_conc,
        dna_conc: dna_nm,
        dmso_conc,
        dmso_fact,
        formamide_conc,
        ..Default::default()
    };
    tm::oligo_annealing(sequence.as_bytes(), &args, annealing_temp_c)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
}

#[pymodule]
fn primalschemers(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Object classes
    m.add_class::<kmer::FKmer>()?;
    m.add_class::<kmer::RKmer>()?;
    m.add_class::<Digester>()?;
    m.add_class::<thal::PyThermoResult>()?;

    // Functions
    m.add_function(wrap_pyfunction!(digest_seq, m)?)?;
    m.add_function(wrap_pyfunction!(calc_at_offset_py, m)?)?;
    m.add_function(wrap_pyfunction!(do_seqs_interact, m)?)?;
    m.add_function(wrap_pyfunction!(do_pool_interact, m)?)?;
    m.add_function(wrap_pyfunction!(kmer::generate_primerpairs_py, m)?)?;

    // Thermo
    m.add_function(wrap_pyfunction!(calc_annealing, m)?)?;
    m.add_function(wrap_pyfunction!(thal::calc_hairpin, m)?)?;
    m.add_function(wrap_pyfunction!(thal::calc_homodimer, m)?)?;
    m.add_function(wrap_pyfunction!(thal::calc_heterodimer, m)?)?;
    m.add_function(wrap_pyfunction!(thal::calc_hairpin_batch, m)?)?;
    m.add_function(wrap_pyfunction!(thal::calc_heterodimer_batch, m)?)?;
    m.add_function(wrap_pyfunction!(thal::calc_tm, m)?)?;
    m.add_function(wrap_pyfunction!(thal::calc_annealing_profile, m)?)?;
    m.add_function(wrap_pyfunction!(thal::calc_annealing_hetero, m)?)?;

    Ok(())
}
