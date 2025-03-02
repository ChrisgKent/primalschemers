use indicatif::{ParallelProgressIterator, ProgressStyle};
use pyo3::prelude::*;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use std::str::from_utf8;

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

pub fn generate_primerpairs<'a>(
    fkmers: &Vec<&'a FKmer>,
    rkmers: &Vec<&'a RKmer>,
    t: f64,
    amplicion_size_min: usize,
    amplicion_size_max: usize,

    min_overlap: usize,
) -> Vec<(&'a FKmer, &'a RKmer)> {
    //Set up pb
    let progress_bar =
        ProgressStyle::with_template("[{elapsed}] {wide_bar:40.red/white} {pos:>7}/{len:7} {eta}")
            .unwrap();

    // Generate the primer pairs in parallel
    let nested_pp: Vec<Vec<(&FKmer, &RKmer)>> = fkmers
        .par_iter()
        .progress_with_style(progress_bar)
        .map(|fkmer| {
            let rkmer_window = &fkmer.end() + amplicion_size_min;
            // Get the start position of the rkmer window
            let pos_rkmer_start = match rkmers.binary_search_by(|r| r.start().cmp(&rkmer_window)) {
                Ok(mut pos) => {
                    while rkmers[pos].start() == rkmer_window && pos > 0 {
                        pos -= 1;
                    }
                    pos
                }
                Err(pos) => pos,
            };

            let max_index = &fkmer.end() + amplicion_size_max;

            let mut primer_pairs: Vec<(&FKmer, &RKmer)> = Vec::new();
            for i in pos_rkmer_start..rkmers.len() {
                let rkmer = &rkmers[i];
                if rkmer.start() > max_index {
                    break;
                }
                // if do_kmers_interact(*fkmer, *rkmer, t) {
                //     primer_pairs.push((fkmer, rkmer));
                // }
            }
            primer_pairs
        })
        .collect();

    nested_pp.into_par_iter().flatten().collect()
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

    #[test]
    fn test_do_kmers_interact() {
        let seqs1 = vec![b"ACACCTGTGCCTGTTAAACCAT".to_vec()];
        let seqs2 = vec![b"TGGAAATACCCACAAGTTAATGGTTTAAC".to_vec()];
        let fkmer = FKmer::new(seqs1, 100);
        let rkmer = RKmer::new(seqs2, 100);

        // assert_eq!(do_kmers_interact(&fkmer, &rkmer, -26.0), true);
    }
    #[test]
    fn test_do_kmers_not_interact() {
        let seqs1 = vec![b"CCAAACAAAGTTGGGTAAGGATCGA".to_vec()];
        let seqs2 = vec![b"ACTGGTCACTCTGAACAACCTCC".to_vec()];
        let fkmer = FKmer::new(seqs1, 100);
        let rkmer = RKmer::new(seqs2, 100);

        // assert_eq!(do_kmers_interact(&fkmer, &rkmer, -26.0), false);
    }
}
