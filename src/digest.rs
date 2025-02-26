use crate::kmer::{FKmer, RKmer};
use crate::seqfuncs::{
    atcg_only, complement_base, contains_ambs, expand_amb_base, expand_amb_sequence, gc_content,
    max_homopolymer,
};
use crate::tm;

use std::collections::HashMap;

use indicatif::{ParallelProgressIterator, ProgressStyle};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

static MIN_PRIMER_LEN: usize = 19;
static MAX_PRIMER_LEN: usize = 50;
static MAX_WALK: usize = 100;
static PRIMER_TM_MAX: f64 = 62.5;
static PRIMER_TM_MIN: f64 = 59.5;
static MAX_HOMOPOLYMERS: usize = 4;

#[derive(Eq, PartialEq, Hash, Debug, Clone)]
pub enum DigestError {
    InvalidBase,
    WalkedOutLeft,
    WalkedOutRight,
    GapOnSetBase,
    NoValidPrimer,
    MaxWalk,
}

#[derive(Clone, Debug)]
pub enum ThermoResult {
    Pass,
    HighGC,
    LowGC,
    Hairpin,
    Homopolymer,
    HighTm,
    LowTm,
}

#[derive(Clone, Debug)]

pub enum IndexResult {
    ThermoResult(ThermoResult),
    DigestError(DigestError),
    Pass,
}

pub struct DigestionKmers {
    pub seq: Option<Vec<u8>>,
    pub status: Option<IndexResult>,
    pub count: f64,
}
impl DigestionKmers {
    pub fn new(seq: Option<Vec<u8>>, status: Option<IndexResult>, count: f64) -> DigestionKmers {
        DigestionKmers {
            seq: seq,
            status: status,
            count: count,
        }
    }
    pub fn calc_freq(&self, total: f64) -> f64 {
        self.count / total
    }
    pub fn thermo_check(&mut self) -> &Option<IndexResult> {
        // Returns the None if valid primer
        let result = match &self.seq {
            Some(seq) => match thermo_check(&seq) {
                ThermoResult::Pass => IndexResult::Pass,
                _ => IndexResult::ThermoResult(thermo_check(&seq)),
            },
            // If
            None => IndexResult::DigestError(DigestError::NoValidPrimer),
        };
        self.status = Some(result);
        &self.status
    }
    pub fn status(&self) -> &Option<IndexResult> {
        &self.status
    }
    pub fn sequence(&self) -> &Option<Vec<u8>> {
        &self.seq
    }
    fn rc(&mut self) -> &Option<Vec<u8>> {
        self.seq = match &self.seq {
            Some(seq) => Some(seq.iter().map(|b| complement_base(*b)).collect()),
            None => None,
        };
        &self.seq
    }
}

fn thermo_check(kmer: &[u8]) -> ThermoResult {
    // Check GC content
    let gc = gc_content(kmer);
    if gc < 0.4 {
        return ThermoResult::LowGC;
    }
    if gc > 0.6 {
        return ThermoResult::HighGC;
    }

    // Check homopolymer
    if max_homopolymer(kmer) > MAX_HOMOPOLYMERS {
        return ThermoResult::Homopolymer;
    }

    // Check hairpin
    // Check Tm
    let (tm, _an) = tm::oligotm_utf8(
        kmer,
        15.0,
        100.0,
        2.0,
        0.8,
        0.0,
        0.0,
        0.8,
        0.0,
        tm::TmMethod::SantaLucia2004,
    );
    if tm >= PRIMER_TM_MAX {
        return ThermoResult::HighTm;
    }
    if tm < PRIMER_TM_MIN {
        return ThermoResult::LowTm;
    }

    //TODO calc hairpin
    ThermoResult::Pass
}

fn process_seqs(seq_counts: HashMap<Result<Vec<u8>, DigestError>, f64>) -> Vec<DigestionKmers> {
    let mut digested: Vec<DigestionKmers> = Vec::new();
    for (k, v) in seq_counts.into_iter() {
        let dk = match k {
            Ok(s) => {
                // If sequence provided check thermo
                let mut dk = DigestionKmers::new(Some(s), None, v as f64);
                dk.thermo_check();
                dk
            }
            Err(e) => {
                let dk = DigestionKmers::new(None, Some(IndexResult::DigestError(e)), v as f64);
                dk
            }
        };
        digested.push(dk);
    }

    // Apply a freq filter
    let total: f64 = digested.iter().map(|d| d.count).sum();
    digested.retain(|d| d.calc_freq(total) >= 0.01);

    digested
}

pub fn walk_right(
    seq: &[u8],
    l_index: usize,
    r_index: usize,
    kmer: tm::Oligo,
) -> Vec<Result<Vec<u8>, DigestError>> {
    // Check tm
    let tm = kmer.calc_tm(15.0, 100.0, 2.0, 0.8, 0.0, 0.0, 0.8);

    if tm >= PRIMER_TM_MIN {
        return vec![Ok(kmer.seq)];
    }

    // Check if we've reached the end of the sequence
    if r_index < l_index || r_index - l_index >= MAX_WALK {
        return vec![Err(DigestError::MaxWalk)];
    }
    // Check bounds
    if r_index == seq.len() {
        return vec![Err(DigestError::WalkedOutRight)];
    }

    let new_base = seq[r_index];

    // If base is gap keep walking
    if new_base == b'-' {
        let new_results = walk_right(seq, l_index, r_index + 1, kmer);
        return new_results;
    }

    // if base is ambiguous, expand it
    let new_bases = match expand_amb_base(new_base) {
        Some(bases) => bases,
        None => return vec![Err(DigestError::InvalidBase)],
    };

    let mut results: Vec<Result<Vec<u8>, DigestError>> = Vec::new();

    // Clone the kmer
    let mut kmer_clones = Vec::new();
    for _ in 0..new_bases.len() - 1 {
        kmer_clones.push(kmer.clone());
    }
    kmer_clones.push(kmer);

    for (base, mut kmer_c) in new_bases.iter().zip(kmer_clones) {
        kmer_c.add_base(*base);
        let new_results = walk_right(seq, l_index, r_index + 1, kmer_c);
        results.extend(new_results);
    }
    results
}

pub fn digest_r_to_count(
    seqs: &Vec<&[u8]>,
    index: usize,
) -> HashMap<Result<Vec<u8>, DigestError>, f64> {
    // Returns the kmer as read (left -> right) in the genomes

    let mut kmer_count: HashMap<Result<Vec<u8>, DigestError>, f64> = HashMap::new();

    // Check bounds
    if index > seqs[0].len() - MIN_PRIMER_LEN {
        kmer_count.insert(Err(DigestError::WalkedOutRight), seqs.len() as f64);
        return kmer_count;
    }
    let rhs = index + MIN_PRIMER_LEN;

    // For each sequence, digest at the index
    for seq in seqs.iter() {
        // Check for gap on set base
        if seq[index] == b'-' {
            let c = kmer_count
                .entry(Err(DigestError::GapOnSetBase))
                .or_insert(0.0);
            *c += 1.0;
            continue;
        }

        // Create the kmer slice
        let mut kmer: Vec<u8> = Vec::with_capacity(MAX_PRIMER_LEN);
        kmer.extend_from_slice(&seq[index..rhs]);

        // Remove gaps from the kmer
        kmer = kmer
            .into_iter()
            .filter(|b| *b != b'-' && *b != b' ')
            .collect();

        if atcg_only(&kmer) {
            // No ambiguous bases
            let results = walk_right(seq, index, rhs, tm::Oligo::new(kmer));
            for r in results {
                let count = kmer_count.entry(r).or_insert(0.0);
                *count += 1.0;
            }
        } else if contains_ambs(&kmer) {
            // Handle ambiguous bases
            let expanded_kmer = expand_amb_sequence(&kmer);

            match expanded_kmer {
                None => {
                    let c = kmer_count
                        .entry(Err(DigestError::InvalidBase))
                        .or_insert(0.0);
                    *c += 1.0;
                }
                Some(expanded_kmer) => {
                    let num_kmers = expanded_kmer.len();
                    for ek in expanded_kmer.into_iter() {
                        let results = walk_right(seq, index, rhs, tm::Oligo::new(ek));
                        for r in results {
                            let count = kmer_count.entry(r).or_insert(0.0);
                            *count += 1.0 / num_kmers as f64;
                        }
                    }
                }
            }
        } else {
            let c = kmer_count
                .entry(Err(DigestError::InvalidBase))
                .or_insert(0.0);
            *c += 1.0;
        }
    }

    kmer_count
}

fn digest_r_at_index(seqs: &Vec<&[u8]>, index: usize) -> Result<RKmer, IndexResult> {
    let kmer_count: HashMap<Result<Vec<u8>, DigestError>, f64> = digest_r_to_count(seqs, index);

    // Process the results
    let mut dks = process_seqs(kmer_count);

    // Reverse complement the sequences
    for dk in dks.iter_mut() {
        dk.rc();
    }

    for dk in dks.iter() {
        match &dk.status {
            // Pass
            Some(IndexResult::ThermoResult(ThermoResult::Pass)) => {}
            Some(IndexResult::Pass) => {}
            // Fail
            Some(indexresult) => {
                return Err(indexresult.clone());
            }
            None => {
                return Err(IndexResult::DigestError(DigestError::NoValidPrimer));
            }
        }
    }

    Ok(RKmer::new(
        dks.into_iter().map(|dk| dk.seq.unwrap()).collect(),
        index,
    ))
}

pub fn digest_r_primer(seq_array: &Vec<&[u8]>) -> Vec<Result<RKmer, IndexResult>> {
    let indexes: Vec<usize> = (0..seq_array[0].len() + 1).collect();

    // Check that all sequences are the same length
    for seq in seq_array.iter() {
        if seq.len() != seq_array[0].len() {
            panic!("Sequences are not the same length");
        }
    }

    let progress_bar =
        ProgressStyle::with_template("[{elapsed}] {wide_bar:40.cyan/blue} {pos:>7}/{len:7} {eta}")
            .unwrap();
    let results: Vec<Result<RKmer, IndexResult>> = indexes
        .par_iter()
        .progress_with_style(progress_bar)
        .map(|i| digest_r_at_index(&seq_array, *i))
        .collect();

    results
}

pub fn walk_left(
    seq: &[u8],
    l_index: usize,
    r_index: usize,
    kmer: tm::Oligo,
) -> Vec<Result<Vec<u8>, DigestError>> {
    // Check tm
    let tm = kmer.calc_tm(15.0, 100.0, 2.0, 0.8, 0.0, 0.0, 0.8);

    if tm >= PRIMER_TM_MIN {
        return vec![Ok(kmer.seq)];
    }

    // Check if we've reached the end of the sequence
    if l_index > r_index || r_index - l_index >= MAX_WALK {
        return vec![Err(DigestError::MaxWalk)];
    }
    // Check bounds
    if l_index == 0 {
        return vec![Err(DigestError::WalkedOutLeft)];
    }

    let new_base = seq[l_index - 1];

    // TODO - Check if new_base is a valid base

    // If base is gap keep walking
    if new_base == b'-' {
        let new_results = walk_left(seq, l_index - 1, r_index, kmer);
        return new_results;
    }

    // if base is ambiguous, expand it
    let new_bases = match expand_amb_base(new_base) {
        Some(bases) => bases,
        None => return vec![Err(DigestError::InvalidBase)],
    };

    let mut results: Vec<Result<Vec<u8>, DigestError>> = Vec::new();

    // Clone the kmer
    let mut kmer_clones = Vec::new();
    for _ in 0..new_bases.len() - 1 {
        kmer_clones.push(kmer.clone());
    }
    kmer_clones.push(kmer);

    for (base, mut kmer_c) in new_bases.iter().zip(kmer_clones) {
        kmer_c.add_base(*base);
        let new_results = walk_left(seq, l_index - 1, r_index, kmer_c);
        results.extend(new_results);
    }
    results
}

pub fn digest_f_to_count(
    seqs: &Vec<&[u8]>,
    index: usize,
) -> HashMap<Result<Vec<u8>, DigestError>, f64> {
    let mut kmer_count: HashMap<Result<Vec<u8>, DigestError>, f64> = HashMap::new();

    // Check bounds
    if index < MIN_PRIMER_LEN {
        kmer_count.insert(Err(DigestError::WalkedOutLeft), seqs.len() as f64);
        return kmer_count;
    }
    let lhs = index - MIN_PRIMER_LEN;

    // For each sequence, digest at the index
    for seq in seqs.iter() {
        // Check for gap on set base
        if seq[index - 1] == b'-' {
            let c = kmer_count
                .entry(Err(DigestError::GapOnSetBase))
                .or_insert(0.0);
            *c += 1.0;
            continue;
        }
        // Create the kmer slice
        let mut kmer: Vec<u8> = Vec::with_capacity(MAX_PRIMER_LEN);
        kmer.extend_from_slice(&seq[lhs..index]);
        kmer.reverse(); // Reverse is used as push is .O(1) and .insert(0) is O(n)

        // Remove gaps from the kmer
        kmer = kmer
            .into_iter()
            .filter(|b| *b != b'-' && *b != b' ')
            .collect();

        if atcg_only(&kmer) {
            // No ambiguous bases
            let results = walk_left(seq, lhs, index, tm::Oligo::new(kmer));
            for r in results {
                let count = kmer_count.entry(r).or_insert(0.0);
                *count += 1.0;
            }
        } else if contains_ambs(&kmer) {
            // Handle ambiguous bases
            let expanded_kmer = expand_amb_sequence(&kmer);
            match expanded_kmer {
                None => {
                    let c = kmer_count
                        .entry(Err(DigestError::InvalidBase))
                        .or_insert(0.0);
                    *c += 1.0;
                }
                Some(expanded_kmer) => {
                    let num_kmers = expanded_kmer.len();
                    for ek in expanded_kmer.into_iter() {
                        let results = walk_left(seq, lhs, index, tm::Oligo::new(ek));
                        for r in results {
                            let count = kmer_count.entry(r).or_insert(0.0);
                            *count += 1.0 / num_kmers as f64;
                        }
                    }
                }
            }
        } else {
            let c = kmer_count
                .entry(Err(DigestError::InvalidBase))
                .or_insert(0.0);
            *c += 1.0;
        }
    }
    // Un-reverse the kmer
    let mut un_reversed: HashMap<Result<Vec<u8>, DigestError>, f64> = HashMap::new();
    for (k, v) in kmer_count.into_iter() {
        match k {
            Ok(mut kmer) => {
                kmer.reverse();
                un_reversed.insert(Ok(kmer), v);
            }
            Err(dr) => {
                un_reversed.insert(Err(dr), v);
            }
        }
    }
    un_reversed
}

fn digest_f_at_index(seqs: &Vec<&[u8]>, index: usize) -> Result<FKmer, IndexResult> {
    let kmer_count: HashMap<Result<Vec<u8>, DigestError>, f64> = digest_f_to_count(seqs, index);

    // Process the results
    let dks = process_seqs(kmer_count);

    // If any errors are found, return the error
    for dk in dks.iter() {
        match &dk.status {
            // Pass
            Some(IndexResult::ThermoResult(ThermoResult::Pass)) => {}
            Some(IndexResult::Pass) => {}
            // Fail
            Some(indexresult) => {
                return Err(indexresult.clone());
            }
            None => {
                return Err(IndexResult::DigestError(DigestError::NoValidPrimer));
            }
        }
    }

    // Create the FKmer
    Ok(FKmer::new(
        dks.into_iter().map(|dk| dk.seq.unwrap()).collect(),
        index,
    ))
}

pub fn digest_f_primer(seq_array: &Vec<&[u8]>) -> Vec<Result<FKmer, IndexResult>> {
    let indexes: Vec<usize> = (0..seq_array[0].len() + 1).collect();
    // Check that all sequences are the same length
    for seq in seq_array.iter() {
        if seq.len() != seq_array[0].len() {
            panic!("Sequences are not the same length");
        }
    }

    let progress_bar =
        ProgressStyle::with_template("[{elapsed}] {wide_bar:40.cyan/blue} {pos:>7}/{len:7} {eta}")
            .unwrap();

    let results: Vec<Result<FKmer, IndexResult>> = indexes
        .par_iter()
        .progress_with_style(progress_bar)
        .map(|i| digest_f_at_index(&seq_array, *i))
        .collect();

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_digest_r_to_count() {
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT".as_bytes()];

        let digested = digest_r_to_count(&seqs, 30);
        // Check num of seqs
        assert_eq!(digested.len(), 1);

        // Check sequence
        let exp_seq: Result<Vec<u8>, DigestError> =
            Ok("ACCAACCAACTTTCGATCTCTTGTAGA".as_bytes().to_vec());

        // Check sequence
        for key in digested.keys() {
            println!(
                "{:?}",
                match key {
                    Ok(k) => String::from_utf8(k.clone()).unwrap(),
                    Err(e) => format!("{:?}", e),
                }
            );
        }

        assert_eq!(digested.contains_key(&exp_seq), true);
        // Check count
        assert_eq!(digested.get(&exp_seq).unwrap(), &1.0);
    }

    #[test]
    fn test_digest_f_to_count() {
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT".as_bytes()];

        let digested = digest_f_to_count(&seqs, 40);
        // Check num of seqs
        assert_eq!(digested.len(), 1);

        // Check sequence
        let exp_seq: Result<Vec<u8>, DigestError> =
            Ok("TCCCAGGTAACAAACCAACCAAC".as_bytes().to_vec());

        assert_eq!(digested.contains_key(&exp_seq), true);
        // Check count
        assert_eq!(digested.get(&exp_seq).unwrap(), &1.0);
    }

    #[test]
    fn test_digest_f_to_count_ps3() {
        let seqs =
            vec!["CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT".as_bytes()];

        let digested = digest_f_to_count(&seqs, 60);
        // Check num of seqs
        assert_eq!(digested.len(), 1);

        // Check sequence
        let exp_seq = Ok("TGTCCAATGGTGCAAAAGGTATAATCATTAAT".as_bytes().to_vec());

        assert_eq!(digested.contains_key(&exp_seq), true);
        // Check count
        assert_eq!(digested.get(&exp_seq).unwrap(), &1.0);
    }

    #[test]
    fn test_digest_f_to_count_gap() {
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAA-TTTCGATCTCTTGTAGATCT".as_bytes()];

        let digested = digest_f_to_count(&seqs, 40);
        // Check num of seqs
        assert_eq!(digested.len(), 1);
        // Check sequence
        assert_eq!(digested.contains_key(&Err(DigestError::GapOnSetBase)), true);
    }

    #[test]
    fn test_digest_f_to_count_wl() {
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAA-TTTCGATCTCTTGTAGATCT".as_bytes()];

        let digested = digest_f_to_count(&seqs, 4);
        // Check num of seqs
        assert_eq!(digested.len(), 1);
        // Check sequence
        assert_eq!(
            digested.contains_key(&Err(DigestError::WalkedOutLeft)),
            true
        );
    }

    #[test]
    fn test_digest_f_to_count_invalid_base() {
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAANTTTCGATCTCTTGTAGATCT".as_bytes()];

        let digested = digest_f_to_count(&seqs, 40);
        // Check num of seqs
        assert_eq!(digested.len(), 1);
        // Check sequence
        assert_eq!(digested.contains_key(&Err(DigestError::InvalidBase)), true);
    }
}
