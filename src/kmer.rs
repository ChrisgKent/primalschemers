use pyo3::prelude::*;
use std::str::from_utf8;

trait Kmer {
    fn seqs(&self) -> Vec<String>;
    fn seqs_bytes(&self) -> Vec<&[u8]>;
}

impl Kmer for FKmer {
    fn seqs(&self) -> Vec<String> {
        // Return the sequences as strings
        self.seqs
            .iter()
            .map(|s| from_utf8(s).unwrap().to_string())
            .collect()
    }
    fn seqs_bytes(&self) -> Vec<&[u8]> {
        // Return the sequences as utf8 bytes
        self.seqs.iter().map(|s| s.as_slice()).collect()
    }
}
impl Kmer for RKmer {
    fn seqs(&self) -> Vec<String> {
        // Return the sequences as strings
        self.seqs
            .iter()
            .map(|s| from_utf8(s).unwrap().to_string())
            .collect()
    }
    fn seqs_bytes(&self) -> Vec<&[u8]> {
        // Return the sequences as utf8 bytes
        self.seqs.iter().map(|s| s.as_slice()).collect()
    }
}

#[derive(Debug)]
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
    #[getter]
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
    #[getter]
    pub fn lens(&self) -> Vec<usize> {
        self.seqs.iter().map(|s| s.len()).collect()
    }
    #[getter]
    pub fn seqs(&self) -> Vec<String> {
        // Return the sequences as strings
        self.seqs
            .iter()
            .map(|s| from_utf8(s).unwrap().to_string())
            .collect()
    }
    #[getter]
    pub fn seqs_bytes(&self) -> Vec<&[u8]> {
        // Return the sequences as utf8 bytes
        self.seqs.iter().map(|s| s.as_slice()).collect()
    }
}

#[pyclass]
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
    #[getter]
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
    #[getter]
    pub fn ends(&self) -> Vec<usize> {
        self.seqs.iter().map(|s| self.start + s.len()).collect()
    }
    #[getter]
    pub fn lens(&self) -> Vec<usize> {
        self.seqs.iter().map(|s| s.len()).collect()
    }
    #[getter]
    pub fn seqs_bytes(&self) -> Vec<&[u8]> {
        // Return the sequences as utf8 bytes
        self.seqs.iter().map(|s| s.as_slice()).collect()
    }
}

fn do_kmers_interact(kmer1: &impl Kmer, kmer2: &impl Kmer) -> bool {
    // Check if the kmers interact
    let seqs1 = kmer1.seqs_bytes();
    let seqs2 = kmer2.seqs_bytes();
    for s1 in seqs1.iter() {
        for s2 in seqs2.iter() {
            if s1 == s2 {
                return true;
            }
        }
    }
    false
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
