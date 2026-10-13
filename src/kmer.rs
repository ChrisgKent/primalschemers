use indicatif::{ParallelProgressIterator, ProgressBar, ProgressIterator, ProgressStyle};
use pyo3::prelude::*;
use rayon::prelude::*;
use std::str::from_utf8;

use crate::primaldimer;

/// Status tag for sequences below `min_freq` (dropped by a normal digest).
pub const LOW_FREQ: &str = "LowFreq";
/// Status tag for errors a normal digest skips, e.g. `ContainsN` with
/// `ignore_n`. Reported for visibility but doesn't fail the kmer.
pub const IGNORED: &str = "Ignored";

/// Everything a `return_all` digest reports for one kmer: each sequence
/// with its count and failed checks, kmer-level failures, and alignment
/// positions that produced no sequence.
pub struct KmerReport {
    pub seqs: Vec<Vec<u8>>,
    pub counts: Vec<f64>,
    /// Failed checks per sequence (empty = passed), aligned with `seqs`.
    pub statuses: Vec<Vec<String>>,
    /// Failures of the kmer as a whole, e.g. `PrimerDimer`.
    pub kmer_statuses: Vec<String>,
    /// Positions with no sequence: their errors and count.
    pub errors: Vec<(Vec<String>, f64)>,
    /// Interacting sequence pairs `(i, j)`, `i <= j`, indexing `seqs`.
    pub dimers: Vec<(usize, usize)>,
}

/// Per-sequence status fields of a kmer built by a `return_all` digest.
#[derive(Debug, Default, PartialEq, PartialOrd)]
struct Statuses {
    per_seq: Vec<Vec<String>>,
    kmer: Vec<String>,
    errors: Vec<(Vec<String>, f64)>,
    /// Interacting pairs, indexing the kmer's sorted `seqs`.
    dimers: Vec<(usize, usize)>,
}

impl Statuses {
    /// Would a normal digest keep this kmer? Low-frequency and ignored
    /// entries don't count, as the normal digest drops or skips them.
    fn passed(&self) -> bool {
        let kept = |tags: &Vec<String>| !tags.iter().any(|t| t == LOW_FREQ || t == IGNORED);
        self.kmer.is_empty()
            && self.per_seq.iter().any(kept)
            && self.per_seq.iter().filter(|t| kept(t)).all(|t| t.is_empty())
            && self.errors.iter().filter(|(t, _)| kept(t)).count() == 0
    }
}

/// Sort a report's sequences, keeping counts, statuses and dimer indices
/// aligned. Sequences from a digest are unique, so no merging is needed.
fn sorted_report(report: KmerReport) -> (Vec<Vec<u8>>, Vec<f64>, Statuses) {
    let mut rows: Vec<(Vec<u8>, f64, Vec<String>, usize)> = report
        .seqs
        .into_iter()
        .zip(report.counts)
        .zip(report.statuses)
        .enumerate()
        .map(|(orig, ((seq, count), status))| (seq, count, status, orig))
        .collect();
    rows.sort_by(|a, b| a.0.cmp(&b.0));
    let mut seqs = Vec::with_capacity(rows.len());
    let mut counts = Vec::with_capacity(rows.len());
    let mut per_seq = Vec::with_capacity(rows.len());
    // new_index[original index] = position after sorting
    let mut new_index = vec![0; rows.len()];
    for (pos, (seq, count, status, orig)) in rows.into_iter().enumerate() {
        seqs.push(seq);
        counts.push(count);
        per_seq.push(status);
        new_index[orig] = pos;
    }
    let mut dimers: Vec<(usize, usize)> = report
        .dimers
        .into_iter()
        .map(|(a, b)| {
            let (a, b) = (new_index[a], new_index[b]);
            (a.min(b), a.max(b))
        })
        .collect();
    dimers.sort();
    let statuses = Statuses {
        per_seq,
        kmer: report.kmer_statuses,
        errors: report.errors,
        dimers,
    };
    (seqs, counts, statuses)
}

#[derive(Debug)]
#[pyclass]
pub struct FKmer {
    seqs: Vec<Vec<u8>>,
    counts: Option<Vec<f64>>,
    end: usize,
    /// Set only by a `return_all` digest.
    statuses: Option<Statuses>,
}

// Equality and order compare only the sequences and position, ignoring
// counts and statuses, so they agree with each other.
impl PartialEq for FKmer {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}
impl Eq for FKmer {}
impl PartialOrd for FKmer {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for FKmer {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Compare only seqs and end, ignoring counts
        match self.seqs.cmp(&other.seqs) {
            std::cmp::Ordering::Equal => self.end.cmp(&other.end),
            ordering => ordering,
        }
    }
}
impl FKmer {
    /// Sorted, deduplicated sequences (summing the counts of duplicates).
    /// Panics if `counts` has a different length to `seqs`.
    pub fn new(seqs: Vec<Vec<u8>>, end: usize, counts: Option<Vec<f64>>) -> FKmer {
        let (sorted_seqs, sorted_counts) = match counts {
            Some(counts_vec) => {
                if counts_vec.iter().len() != seqs.iter().len() {
                    panic!("Different number of seqs and counts")
                }
                // Create pairs of (seq, count), sort by seq, then separate
                let mut pairs: Vec<(Vec<u8>, f64)> =
                    seqs.into_iter().zip(counts_vec.into_iter()).collect();
                pairs.sort_by(|a, b| a.0.cmp(&b.0));

                // Remove duplicates while preserving counts (sum counts for duplicates)
                let mut unique_pairs: Vec<(Vec<u8>, f64)> = Vec::with_capacity(pairs.len());
                for (seq, count) in pairs {
                    if let Some((last_seq, last_count)) = unique_pairs.last_mut() {
                        if *last_seq == seq {
                            *last_count += count;
                            continue;
                        }
                    }
                    unique_pairs.push((seq, count));
                }

                let (seqs, counts): (Vec<Vec<u8>>, Vec<f64>) = unique_pairs.into_iter().unzip();
                (seqs, Some(counts))
            }
            None => {
                let mut sorted_seqs = seqs;
                sorted_seqs.sort();
                sorted_seqs.dedup();
                (sorted_seqs, None)
            }
        };
        FKmer {
            seqs: sorted_seqs,
            end: end,
            counts: sorted_counts,
            statuses: None,
        }
    }

    /// A kmer from a `return_all` digest, carrying every sequence's status.
    pub fn from_report(report: KmerReport, end: usize) -> FKmer {
        let (seqs, counts, statuses) = sorted_report(report);
        FKmer {
            seqs,
            counts: Some(counts),
            end,
            statuses: Some(statuses),
        }
    }
}

#[pymethods]
impl FKmer {
    /// ValueError if `counts` is given with a different length to `seqs`.
    #[new]
    #[pyo3(signature = (seqs, end, counts=None,))]
    fn py_new(seqs: Vec<Vec<u8>>, end: usize, counts: Option<Vec<f64>>) -> PyResult<FKmer> {
        if counts.as_ref().is_some_and(|c| c.len() != seqs.len()) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "Different number of seqs and counts",
            ));
        }
        Ok(FKmer::new(seqs, end, counts))
    }
    pub fn starts(&self) -> Vec<usize> {
        // Returns the start positions of the sequences.
        self.seqs
            .iter()
            .map(|s| match self.end.checked_sub(s.len()) {
                Some(s) => s,
                None => 0,
            })
            .collect()
    }
    #[getter]
    pub fn end(&self) -> usize {
        self.end
    }
    pub fn counts(&self) -> &Option<Vec<f64>> {
        &self.counts
    }
    pub fn len(&self) -> Vec<usize> {
        self.seqs.iter().map(|s| s.len()).collect()
    }
    pub fn num_seqs(&self) -> usize {
        self.seqs.len()
    }
    pub fn seqs(&self) -> Vec<String> {
        // Return the sequences as strings
        self.seqs
            .iter()
            .map(|s| from_utf8(s).unwrap().to_string())
            .collect()
    }
    pub fn seqs_bytes(&self) -> Vec<&[u8]> {
        // Return the sequences as utf8 bytes
        self.seqs.iter().map(|s| s.as_slice()).collect()
    }
    pub fn region(&self) -> (usize, usize) {
        (*self.starts().iter().min().unwrap(), self.end)
    }
    pub fn to_bed(&self, chrom: String, amplicon_prefix: String, pool: usize) -> String {
        let mut string = String::new();
        for (index, seq) in self.seqs().iter().enumerate() {
            let start_pos = match self.end.checked_sub(seq.len()) {
                Some(pos) => pos,
                None => 0,
            };
            let seq_count = self.counts.as_ref().and_then(|c| c.get(index));
            let attr_str: String;

            match seq_count {
                Some(c) => {
                    attr_str = format!("\tpc={}", c);
                }
                None => attr_str = "".to_string(),
            }
            string.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}{}\n",
                chrom,
                start_pos,
                self.end,
                format!("{}_LEFT_{}", amplicon_prefix, index + 1), // +1 for 1-based indexing of primer_number
                pool,
                "+",
                seq,
                attr_str
            ));
        }
        string
    }
    pub fn remap(&mut self, end: usize) {
        self.end = end;
    }
    pub fn gc(&self) -> Vec<f64> {
        // Calculate the GC content of each sequence
        self.seqs
            .iter()
            .map(|s| {
                let gc_count = s.iter().filter(|&&b| b == b'G' || b == b'C').count();
                gc_count as f64 / s.len() as f64
            })
            .collect()
    }

    pub fn counts_sum(&self) -> Option<f64> {
        self.counts.as_ref().map(|c| c.iter().sum())
    }

    /// Failed checks per sequence, aligned with `seqs()` (empty = passed);
    /// None unless digested with `return_all`.
    pub fn statuses(&self) -> Option<Vec<Vec<String>>> {
        self.statuses.as_ref().map(|s| s.per_seq.clone())
    }

    /// Failures of the kmer as a whole (e.g. `PrimerDimer`).
    pub fn kmer_statuses(&self) -> Vec<String> {
        self.statuses.as_ref().map(|s| s.kmer.clone()).unwrap_or_default()
    }

    /// Alignment positions that produced no sequence: (errors, count).
    /// Errors a normal digest skips are tagged `Ignored` (e.g. `ContainsN`
    /// with `ignore_n`) and don't affect `passed()`.
    pub fn errors(&self) -> Vec<(Vec<String>, f64)> {
        self.statuses.as_ref().map(|s| s.errors.clone()).unwrap_or_default()
    }

    /// Sequence pairs that form a primer dimer, as indices into `seqs()`
    /// with `i <= j` (`(i, i)` is a self-dimer). Empty unless digested with
    /// `return_all`; non-empty exactly when `kmer_statuses()` has `PrimerDimer`.
    pub fn dimers(&self) -> Vec<(usize, usize)> {
        self.statuses.as_ref().map(|s| s.dimers.clone()).unwrap_or_default()
    }

    /// Whether a normal digest would keep this kmer. Always true for kmers
    /// not digested with `return_all`.
    pub fn passed(&self) -> bool {
        self.statuses.as_ref().is_none_or(Statuses::passed)
    }
}

#[pyclass]
#[derive(Debug)]
pub struct RKmer {
    seqs: Vec<Vec<u8>>,
    counts: Option<Vec<f64>>,
    start: usize,
    /// Set only by a `return_all` digest.
    statuses: Option<Statuses>,
}

// Equality and order compare only the sequences and position, ignoring
// counts and statuses, so they agree with each other.
impl PartialEq for RKmer {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}
impl Eq for RKmer {}
impl PartialOrd for RKmer {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for RKmer {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Compare only seqs and start, ignoring counts
        match self.seqs.cmp(&other.seqs) {
            std::cmp::Ordering::Equal => self.start.cmp(&other.start),
            ordering => ordering,
        }
    }
}

impl RKmer {
    /// Sorted, deduplicated sequences (summing the counts of duplicates).
    /// Panics if `counts` has a different length to `seqs`.
    pub fn new(seqs: Vec<Vec<u8>>, start: usize, counts: Option<Vec<f64>>) -> RKmer {
        let (sorted_seqs, sorted_counts) = match counts {
            Some(counts_vec) => {
                if counts_vec.iter().len() != seqs.iter().len() {
                    panic!("Different number of seqs and counts")
                }

                // Create pairs of (seq, count), sort by seq, then separate
                let mut pairs: Vec<(Vec<u8>, f64)> =
                    seqs.into_iter().zip(counts_vec.into_iter()).collect();
                pairs.sort_by(|a, b| a.0.cmp(&b.0));

                // Remove duplicates while preserving counts (sum counts for duplicates)
                let mut unique_pairs: Vec<(Vec<u8>, f64)> = Vec::with_capacity(pairs.len());
                for (seq, count) in pairs {
                    if let Some((last_seq, last_count)) = unique_pairs.last_mut() {
                        if *last_seq == seq {
                            *last_count += count;
                            continue;
                        }
                    }
                    unique_pairs.push((seq, count));
                }

                let (seqs, counts): (Vec<Vec<u8>>, Vec<f64>) = unique_pairs.into_iter().unzip();
                (seqs, Some(counts))
            }
            None => {
                let mut sorted_seqs = seqs;
                sorted_seqs.sort();
                sorted_seqs.dedup();
                (sorted_seqs, None)
            }
        };
        RKmer {
            seqs: sorted_seqs,
            counts: sorted_counts,
            start: start,
            statuses: None,
        }
    }

    /// A kmer from a `return_all` digest, carrying every sequence's status.
    pub fn from_report(report: KmerReport, start: usize) -> RKmer {
        let (seqs, counts, statuses) = sorted_report(report);
        RKmer {
            seqs,
            counts: Some(counts),
            start,
            statuses: Some(statuses),
        }
    }
}

#[pymethods]
impl RKmer {
    /// ValueError if `counts` is given with a different length to `seqs`.
    #[new]
    #[pyo3(signature = (seqs, start, counts=None,))]
    fn py_new(seqs: Vec<Vec<u8>>, start: usize, counts: Option<Vec<f64>>) -> PyResult<RKmer> {
        if counts.as_ref().is_some_and(|c| c.len() != seqs.len()) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "Different number of seqs and counts",
            ));
        }
        Ok(RKmer::new(seqs, start, counts))
    }
    pub fn seqs(&self) -> Vec<String> {
        // Return the sequences as strings
        self.seqs
            .iter()
            .map(|s| from_utf8(s).unwrap().to_string())
            .collect()
    }
    #[getter]
    pub fn start(&self) -> usize {
        self.start
    }
    pub fn counts(&self) -> &Option<Vec<f64>> {
        &self.counts
    }
    pub fn ends(&self) -> Vec<usize> {
        self.seqs.iter().map(|s| self.start + s.len()).collect()
    }
    pub fn lens(&self) -> Vec<usize> {
        self.seqs.iter().map(|s| s.len()).collect()
    }
    pub fn num_seqs(&self) -> usize {
        self.seqs.len()
    }
    pub fn seqs_bytes(&self) -> Vec<&[u8]> {
        // Return the sequences as utf8 bytes
        self.seqs.iter().map(|s| s.as_slice()).collect()
    }
    pub fn region(&self) -> (usize, usize) {
        (self.start, *self.ends().iter().max().unwrap())
    }
    pub fn to_bed(&self, chrom: String, amplicon_prefix: String, pool: usize) -> String {
        let mut string = String::new();
        for (index, seq) in self.seqs().iter().enumerate() {
            let seq_count = self.counts.as_ref().and_then(|c| c.get(index));
            let attr_str: String;

            match seq_count {
                Some(c) => {
                    attr_str = format!("\tpc={}", c);
                }
                None => attr_str = "".to_string(),
            }

            string.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}{}\n",
                chrom,
                self.start,
                self.start + seq.len(),
                format!("{}_RIGHT_{}", amplicon_prefix, index + 1), // +1 for 1-based indexing
                pool,
                "-",
                seq,
                attr_str
            ));
        }
        string
    }
    pub fn remap(&mut self, start: usize) {
        self.start = start;
    }

    pub fn gc(&self) -> Vec<f64> {
        // Calculate the GC content of each sequence
        self.seqs
            .iter()
            .map(|s| {
                let gc_count = s.iter().filter(|&&b| b == b'G' || b == b'C').count();
                gc_count as f64 / s.len() as f64
            })
            .collect()
    }

    pub fn counts_sum(&self) -> Option<f64> {
        self.counts.as_ref().map(|c| c.iter().sum())
    }

    /// Failed checks per sequence, aligned with `seqs()` (empty = passed);
    /// None unless digested with `return_all`.
    pub fn statuses(&self) -> Option<Vec<Vec<String>>> {
        self.statuses.as_ref().map(|s| s.per_seq.clone())
    }

    /// Failures of the kmer as a whole (e.g. `PrimerDimer`).
    pub fn kmer_statuses(&self) -> Vec<String> {
        self.statuses.as_ref().map(|s| s.kmer.clone()).unwrap_or_default()
    }

    /// Alignment positions that produced no sequence: (errors, count).
    /// Errors a normal digest skips are tagged `Ignored` (e.g. `ContainsN`
    /// with `ignore_n`) and don't affect `passed()`.
    pub fn errors(&self) -> Vec<(Vec<String>, f64)> {
        self.statuses.as_ref().map(|s| s.errors.clone()).unwrap_or_default()
    }

    /// Sequence pairs that form a primer dimer, as indices into `seqs()`
    /// with `i <= j` (`(i, i)` is a self-dimer). Empty unless digested with
    /// `return_all`; non-empty exactly when `kmer_statuses()` has `PrimerDimer`.
    pub fn dimers(&self) -> Vec<(usize, usize)> {
        self.statuses.as_ref().map(|s| s.dimers.clone()).unwrap_or_default()
    }

    /// Whether a normal digest would keep this kmer. Always true for kmers
    /// not digested with `return_all`.
    pub fn passed(&self) -> bool {
        self.statuses.as_ref().is_none_or(Statuses::passed)
    }
}

#[pyfunction]
pub fn generate_primerpairs_py(
    py: Python<'_>,
    fkmers: Vec<Py<FKmer>>,
    rkmers: Vec<Py<RKmer>>,
    t: f64,
    amplicon_size_min: usize,
    amplicon_size_max: usize,
) -> PyResult<Vec<(Py<FKmer>, Py<RKmer>)>> {
    // Set up pb
    let progress_bar = ProgressBar::new(fkmers.len() as u64);
    progress_bar.set_message("primerpair generation");
    progress_bar.set_style(
        ProgressStyle::default_bar()
            .template("{msg} [{elapsed}] {wide_bar:.cyan/blue} {pos:>7}/{len:7} {eta}")
            .unwrap(),
    );

    // Create two arrays on set bases for rapid GIL free lookup
    let rkmer_ends = rkmers
        .iter()
        .map(|r| r.borrow(py).start())
        .collect::<Vec<usize>>();

    // ensure rkmers are sorted
    if !rkmer_ends.is_sorted() {
        panic!("RKmer list is not sorted")
    }

    // Generate the primer pairs
    let nested_pp: Vec<Vec<(Py<FKmer>, Py<RKmer>)>> = fkmers
        .iter()
        .progress_with(progress_bar)
        .map(|fkmer| {
            let rkmer_window_start = fkmer.borrow(py).end() + amplicon_size_min;
            let rkmer_window_end = &fkmer.borrow(py).end() + amplicon_size_max;

            // Get the start position of the rkmer window
            let pos_rkmer_start = match rkmer_ends.binary_search(&rkmer_window_start) {
                Ok(mut pos) => {
                    while rkmer_ends[pos] >= rkmer_window_start && pos > 0 {
                        pos -= 1;
                    }
                    pos
                }
                Err(pos) => pos,
            };

            let mut primer_pairs: Vec<(Py<FKmer>, Py<RKmer>)> = Vec::new();
            for i in pos_rkmer_start..rkmers.len() {
                let rkmer = &rkmers[i];
                if rkmer.borrow(py).start() > rkmer_window_end {
                    break;
                }
                if primaldimer::do_pool_interact_u8_slice(
                    &fkmer.borrow(py).seqs_bytes(),
                    &rkmer.borrow(py).seqs_bytes(),
                    t,
                ) {
                    primer_pairs.push((fkmer.clone_ref(py), rkmer.clone_ref(py)));
                }
            }
            primer_pairs
        })
        .collect();

    let pp: Vec<(Py<FKmer>, Py<RKmer>)> = nested_pp.into_iter().flatten().collect();
    Ok(pp)
}

pub fn generate_primerpairs<'a>(
    fkmers: &Vec<&'a FKmer>,
    rkmers: &Vec<&'a RKmer>,
    dimerscore: f64,
    amplicon_size_min: usize,
    amplicon_size_max: usize,
) -> Vec<(&'a FKmer, &'a RKmer)> {
    // Set up pb
    let progress_bar = ProgressBar::new(fkmers.len() as u64);
    progress_bar.set_message("primerpair generation");
    progress_bar.set_style(
        ProgressStyle::default_bar()
            .template("{msg} [{elapsed}] {wide_bar:.cyan/blue} {pos:>7}/{len:7} {eta}")
            .unwrap(),
    );

    // ensure rkmers are sorted
    if !rkmers.is_sorted_by(|a, b| a.start() < b.start()) {
        panic!("RKmer list is not sorted")
    }

    // Generate the primer pairs
    let nested_pp: Vec<Vec<(&'a FKmer, &'a RKmer)>> = fkmers
        .par_iter()
        .progress_with(progress_bar)
        .map(|fkmer| {
            let rkmer_window_start = fkmer.end() + amplicon_size_min;
            let rkmer_window_end = fkmer.end() + amplicon_size_max;

            // Get the start position of the rkmer window
            let pos_rkmer_start =
                match rkmers.binary_search_by(|rk| rk.start().cmp(&rkmer_window_start)) {
                    Ok(mut pos) => {
                        while rkmers[pos].start() >= rkmer_window_start && pos > 0 {
                            pos -= 1;
                        }
                        pos
                    }
                    Err(pos) => pos,
                };

            let mut primer_pairs: Vec<(&'a FKmer, &'a RKmer)> = Vec::new();
            for i in pos_rkmer_start..rkmers.len() {
                let rkmer = &rkmers[i];
                if rkmer.start() > rkmer_window_end {
                    break;
                }
                if primaldimer::do_pool_interact_u8(&fkmer.seqs, &rkmer.seqs, dimerscore) {
                    primer_pairs.push((fkmer, rkmer));
                }
            }
            primer_pairs
        })
        .collect();

    let pp: Vec<(&'a FKmer, &'a RKmer)> = nested_pp.into_iter().flatten().collect();
    pp
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_equality_ignores_counts() {
        use std::cmp::Ordering;
        let f = |end, count| FKmer::new(vec![b"ACGT".to_vec()], end, Some(vec![count]));
        assert_eq!(f(10, 1.0), f(10, 5.0));
        assert_eq!(f(10, 1.0).cmp(&f(10, 5.0)), Ordering::Equal);
        assert_ne!(f(10, 1.0), f(11, 1.0));
        let r = |start, count| RKmer::new(vec![b"ACGT".to_vec()], start, Some(vec![count]));
        assert_eq!(r(10, 1.0), r(10, 5.0));
        assert_ne!(r(10, 1.0), r(11, 1.0));
        assert!(r(10, 1.0) < r(11, 1.0));
    }

    fn tags(t: &[&str]) -> Vec<String> {
        t.iter().map(|s| s.to_string()).collect()
    }

    fn report(statuses: Vec<Vec<String>>, kmer: Vec<String>, errors: Vec<(Vec<String>, f64)>) -> KmerReport {
        let n = statuses.len();
        KmerReport {
            seqs: (0..n).map(|i| format!("ACGT{i}").into_bytes()).collect(),
            counts: vec![1.0; n],
            statuses,
            kmer_statuses: kmer,
            errors,
            dimers: vec![],
        }
    }

    #[test]
    fn test_passed_rules() {
        let cases: Vec<(&str, KmerReport, bool)> = vec![
            ("all clean", report(vec![tags(&[]), tags(&[])], vec![], vec![]), true),
            ("one failing sequence", report(vec![tags(&[]), tags(&["Hairpin"])], vec![], vec![]), false),
            ("kmer status", report(vec![tags(&[])], tags(&["PrimerDimer"]), vec![]), false),
            ("error", report(vec![tags(&[])], vec![], vec![(tags(&["WalkedOutLeft"]), 1.0)]), false),
            // A normal digest drops low-frequency entries before checking them
            ("low-freq failing sequence", report(vec![tags(&[]), tags(&["HighGC", "LowFreq"])], vec![], vec![]), true),
            ("low-freq error", report(vec![tags(&[])], vec![], vec![(tags(&["ContainsN", "LowFreq"]), 1.0)]), true),
            ("ignored error", report(vec![tags(&[])], vec![], vec![(tags(&["ContainsN", "Ignored"]), 1.0)]), true),
            ("ignored error alone", report(vec![tags(&["LowFreq"])], vec![], vec![(tags(&["ContainsN", "Ignored"]), 1.0)]), false),
            ("only low-freq sequences", report(vec![tags(&["LowFreq"])], vec![], vec![]), false),
            ("no sequences", report(vec![], vec![], vec![]), false),
        ];
        for (name, r, expected) in cases {
            assert_eq!(FKmer::from_report(clone_report(&r), 10).passed(), expected, "FKmer {name}");
            assert_eq!(RKmer::from_report(r, 10).passed(), expected, "RKmer {name}");
        }
    }

    fn clone_report(r: &KmerReport) -> KmerReport {
        KmerReport {
            seqs: r.seqs.clone(),
            counts: r.counts.clone(),
            statuses: r.statuses.clone(),
            kmer_statuses: r.kmer_statuses.clone(),
            errors: r.errors.clone(),
            dimers: r.dimers.clone(),
        }
    }

    #[test]
    fn test_normal_kmers_have_no_statuses() {
        let f = FKmer::new(vec![b"ACGT".to_vec()], 4, None);
        let r = RKmer::new(vec![b"ACGT".to_vec()], 0, Some(vec![1.0]));
        assert!(f.passed() && r.passed());
        assert_eq!(f.statuses(), None);
        assert_eq!(r.statuses(), None);
        assert!(f.kmer_statuses().is_empty() && f.errors().is_empty() && f.dimers().is_empty());
        assert!(r.kmer_statuses().is_empty() && r.errors().is_empty() && r.dimers().is_empty());
    }

    #[test]
    fn test_rkmer_from_report_aligned() {
        let r = RKmer::from_report(
            KmerReport {
                seqs: vec![b"TTTT".to_vec(), b"AAAA".to_vec()],
                counts: vec![1.0, 2.0],
                statuses: vec![tags(&["LowTm"]), tags(&[])],
                kmer_statuses: tags(&["PrimerDimer"]),
                errors: vec![(tags(&["MaxWalk"]), 3.0)],
                dimers: vec![(0, 1), (0, 0)],
            },
            7,
        );
        assert_eq!(r.start(), 7);
        assert_eq!(r.seqs(), vec!["AAAA", "TTTT"]);
        assert_eq!(r.counts(), &Some(vec![2.0, 1.0]));
        assert_eq!(r.statuses(), Some(vec![tags(&[]), tags(&["LowTm"])]));
        assert_eq!(r.kmer_statuses(), tags(&["PrimerDimer"]));
        assert_eq!(r.errors(), vec![(tags(&["MaxWalk"]), 3.0)]);
        assert_eq!(r.dimers(), vec![(0, 1), (1, 1)]);
        assert!(!r.passed());
    }

    #[test]
    fn test_from_report_keeps_dimer_indices_aligned() {
        // Unsorted input: indices must follow their sequences through the sort.
        let report = KmerReport {
            seqs: vec![b"TTTT".to_vec(), b"AAAA".to_vec(), b"GGGG".to_vec()],
            counts: vec![1.0, 2.0, 3.0],
            statuses: vec![vec![], vec!["HighGC".to_string()], vec![]],
            kmer_statuses: vec!["PrimerDimer".to_string()],
            errors: vec![],
            dimers: vec![(0, 2), (1, 1)], // TTTT-GGGG, AAAA-AAAA
        };
        let k = FKmer::from_report(report, 10);
        assert_eq!(k.seqs(), vec!["AAAA", "GGGG", "TTTT"]);
        assert_eq!(k.counts(), &Some(vec![2.0, 3.0, 1.0]));
        assert_eq!(k.statuses().unwrap()[0], vec!["HighGC".to_string()]);
        assert_eq!(k.dimers(), vec![(0, 0), (1, 2)]);
    }

    #[test]
    fn test_fkmer_start() {
        let seqs = vec![b"ATCG".to_vec()];
        let fkmer = FKmer::new(seqs, 100, None);
        assert_eq!(fkmer.starts(), vec![96]);
    }
    #[test]
    fn test_fkmer_start_lt_zero() {
        let seqs = vec![b"ATCG".to_vec()];
        let fkmer = FKmer::new(seqs, 1, None);
        assert_eq!(fkmer.starts(), vec![0]);
    }
    #[test]
    fn test_fkmer_dedup() {
        let seqs = vec![b"ATCG".to_vec(), b"ATCG".to_vec()];
        let fkmer = FKmer::new(seqs, 100, None);
        assert_eq!(fkmer.seqs().len(), 1);
    }
    #[test]
    fn test_fkmer_counts_order() {
        let seqs = vec![b"GGGG".to_vec(), b"AAAA".to_vec()];
        let counts = vec![10.0, 51.0];

        let fkmer = FKmer::new(seqs.clone(), 100, Some(counts.clone()));

        // Check the counts maintained paired when seqs are reordered
        assert_ne!(seqs, fkmer.seqs);
        assert_ne!(Some(counts.clone()), fkmer.counts);
    }

    #[test]
    fn test_fkmer_counts_dupe() {
        let seqs = vec![b"AAAA".to_vec(), b"AAAA".to_vec()];
        let counts = vec![10.0, 51.0];

        let fkmer = FKmer::new(seqs.clone(), 100, Some(counts.clone()));

        // Check the duplicate counts are summed
        assert_eq!(Some(vec![counts.iter().sum()]), fkmer.counts);
    }

    #[test]
    fn test_rkmer_end() {
        let seqs = vec![b"ATCG".to_vec()];
        let rkmer = RKmer::new(seqs, 100, None);
        assert_eq!(rkmer.ends(), vec![104]);
    }
    #[test]
    fn test_rkmer_end_lt_zero() {
        let seqs = vec![b"ATCG".to_vec()];
        let rkmer = RKmer::new(seqs, 1, None);
        assert_eq!(rkmer.ends(), vec![5]);
    }

    #[test]
    fn test_rkmer_lens() {
        let seqs = vec![b"ATCG".to_vec(), b"ATCG".to_vec()];
        let rkmer = RKmer::new(seqs, 100, None);
        assert_eq!(rkmer.lens(), vec![4]);
    }
    #[test]
    fn test_rkmer_counts_order() {
        let seqs = vec![b"GGGG".to_vec(), b"AAAA".to_vec()];
        let counts = vec![10.0, 51.0];

        let rkmer = RKmer::new(seqs.clone(), 100, Some(counts.clone()));

        // Check the counts maintained paired when seqs are reordered
        assert_ne!(seqs, rkmer.seqs);
        assert_ne!(Some(counts.clone()), rkmer.counts);
    }

    #[test]
    fn test_rkmer_counts_dupe() {
        let seqs = vec![b"AAAA".to_vec(), b"AAAA".to_vec()];
        let counts = vec![10.0, 51.0];

        let rkmer = RKmer::new(seqs.clone(), 100, Some(counts.clone()));

        // Check the duplicate counts are summed
        assert_eq!(Some(vec![counts.iter().sum()]), rkmer.counts);
    }

    #[test]
    fn test_counts_dupe_non_adjacent() {
        // Duplicates separated in the input are merged after sorting
        let seqs = vec![
            b"GGGG".to_vec(),
            b"AAAA".to_vec(),
            b"CCCC".to_vec(),
            b"AAAA".to_vec(),
            b"GGGG".to_vec(),
            b"AAAA".to_vec(),
        ];
        let counts = vec![1.0, 2.0, 4.0, 8.0, 16.0, 32.0];
        let expected_seqs = vec![b"AAAA".to_vec(), b"CCCC".to_vec(), b"GGGG".to_vec()];
        let expected_counts = Some(vec![42.0, 4.0, 17.0]);

        let fkmer = FKmer::new(seqs.clone(), 100, Some(counts.clone()));
        assert_eq!(fkmer.seqs, expected_seqs);
        assert_eq!(fkmer.counts, expected_counts);

        let rkmer = RKmer::new(seqs, 100, Some(counts));
        assert_eq!(rkmer.seqs, expected_seqs);
        assert_eq!(rkmer.counts, expected_counts);
    }
}
