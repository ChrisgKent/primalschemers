use indicatif::{ParallelProgressIterator, ProgressBar, ProgressIterator, ProgressStyle};
use pyo3::prelude::*;
use rayon::prelude::*;
use std::str::from_utf8;

use crate::primaldimer;

#[derive(Debug, PartialEq, PartialOrd)]
#[pyclass]
pub struct FKmer {
    seqs: Vec<Vec<u8>>,
    counts: Option<Vec<f64>>,
    end: usize,
}

// Implement Eq manually by ignoring the counts field
impl Eq for FKmer {}

// Implement Ord manually by ignoring the counts field for comparison
impl Ord for FKmer {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Compare only seqs and end, ignoring counts
        match self.seqs.cmp(&other.seqs) {
            std::cmp::Ordering::Equal => self.end.cmp(&other.end),
            ordering => ordering,
        }
    }
}
#[pymethods]
impl FKmer {
    #[pyo3(signature = (seqs, end, counts=None,))]
    #[new]
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
                let mut unique_pairs: Vec<(Vec<u8>, f64)> = Vec::new();
                for (seq, count) in pairs {
                    if let Some(existing) = unique_pairs
                        .iter_mut()
                        .find(|(existing_seq, _)| existing_seq == &seq)
                    {
                        existing.1 += count; // Sum the counts for duplicate sequences
                    } else {
                        unique_pairs.push((seq, count));
                    }
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
        }
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
            string.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                chrom,
                start_pos,
                self.end,
                format!("{}_LEFT_{}", amplicon_prefix, index + 1), // +1 for 1-based indexing of primer_number
                pool,
                "+",
                seq
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
}

#[pyclass]
#[derive(Debug, PartialEq, PartialOrd)]
pub struct RKmer {
    seqs: Vec<Vec<u8>>,
    counts: Option<Vec<f64>>,
    start: usize,
}

// Implement Eq manually by ignoring the counts field
impl Eq for RKmer {}

// Implement Ord manually by ignoring the counts field for comparison
impl Ord for RKmer {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Compare only seqs and end, ignoring counts
        match self.seqs.cmp(&other.seqs) {
            std::cmp::Ordering::Equal => self.start.cmp(&other.start),
            ordering => ordering,
        }
    }
}

#[pymethods]
impl RKmer {
    #[new]
    #[pyo3(signature = (seqs, start, counts=None,))]
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
                let mut unique_pairs: Vec<(Vec<u8>, f64)> = Vec::new();
                for (seq, count) in pairs {
                    if let Some(existing) = unique_pairs
                        .iter_mut()
                        .find(|(existing_seq, _)| existing_seq == &seq)
                    {
                        existing.1 += count; // Sum the counts for duplicate sequences
                    } else {
                        unique_pairs.push((seq, count));
                    }
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
        }
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
            string.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                chrom,
                self.start,
                self.start + seq.len(),
                format!("{}_RIGHT_{}", amplicon_prefix, index + 1), // +1 for 1-based indexing
                pool,
                "-",
                seq
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
}
