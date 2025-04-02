use indicatif::{ParallelProgressIterator, ProgressBar, ProgressIterator, ProgressStyle};
use pyo3::prelude::*;
use rayon::prelude::*;
use std::str::from_utf8;

use crate::primaldimer;

#[derive(Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
#[pyclass]
pub struct FKmer {
    seqs: Vec<Vec<u8>>,
    end: usize,
}
#[pymethods]
impl FKmer {
    #[new]
    pub fn new(mut seqs: Vec<Vec<u8>>, end: usize) -> FKmer {
        seqs.sort(); // Sort the sequences by base
        seqs.dedup();
        FKmer {
            seqs: seqs,
            end: end,
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
            string.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                chrom,
                self.end - seq.len(),
                self.end,
                format!("{}_LEFT_{}", amplicon_prefix, index),
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
}

#[pyclass]
#[derive(Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct RKmer {
    seqs: Vec<Vec<u8>>,
    start: usize,
}
#[pymethods]
impl RKmer {
    #[new]
    pub fn new(mut seqs: Vec<Vec<u8>>, start: usize) -> RKmer {
        seqs.sort(); // Sort the sequences by base
        seqs.dedup();
        RKmer {
            seqs: seqs,
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
                format!("{}_RIGHT_{}", amplicon_prefix, index),
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
        let fkmer = FKmer::new(seqs, 100);
        assert_eq!(fkmer.starts(), vec![96]);
    }
    #[test]
    fn test_fkmer_start_lt_zero() {
        let seqs = vec![b"ATCG".to_vec()];
        let fkmer = FKmer::new(seqs, 1);
        assert_eq!(fkmer.starts(), vec![0]);
    }
    #[test]
    fn test_fkmer_dedup() {
        let seqs = vec![b"ATCG".to_vec(), b"ATCG".to_vec()];
        let fkmer = FKmer::new(seqs, 100);
        assert_eq!(fkmer.seqs().len(), 1);
    }
    #[test]
    fn test_rkmer_end() {
        let seqs = vec![b"ATCG".to_vec()];
        let rkmer = RKmer::new(seqs, 100);
        assert_eq!(rkmer.ends(), vec![104]);
    }
    #[test]
    fn test_rkmer_end_lt_zero() {
        let seqs = vec![b"ATCG".to_vec()];
        let rkmer = RKmer::new(seqs, 1);
        assert_eq!(rkmer.ends(), vec![5]);
    }

    #[test]
    fn test_rkmer_lens() {
        let seqs = vec![b"ATCG".to_vec(), b"ATCG".to_vec()];
        let rkmer = RKmer::new(seqs, 100);
        assert_eq!(rkmer.lens(), vec![4]);
    }
}
