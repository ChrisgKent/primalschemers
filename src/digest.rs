use crate::config::{DigestConfig, ThermoType};
use crate::kmer::{FKmer, RKmer};
use crate::seqfuncs::{
    check_kmer, complement_base, expand_amb_base, expand_amb_sequence, gc_content, max_homopolymer,
    reverse_complement, KmerCheck,
};
use crate::{primaldimer, tm};

use std::collections::HashMap;

use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

static DNA_NM: f64 = 15.0;
static K_MM: f64 = 100.0;
static DIVALENT_CONC: f64 = 2.0;
static DNTP_CONC: f64 = 0.8;
static DMSO_CONC: f64 = 0.0;
static DMSO_FACT: f64 = 0.0;
static FORMAMIDE_CONC: f64 = 0.8;

static ANNEALING_DIFF: f64 = 5.0;

#[derive(Eq, PartialEq, Hash, Debug, Clone)]
pub enum DigestError {
    InvalidBase,
    WalkedOutLeft,
    WalkedOutRight,
    GapOnSetBase,
    NoValidPrimer,
    MaxWalk,
    EndOfSequence,
    ToLong,
    ContainsN,
}

#[derive(Eq, PartialEq, Hash, Debug, Clone)]
pub enum ThermoResult {
    Pass,
    HighGC,
    LowGC,
    Hairpin,
    Homopolymer,
    HighTm,
    LowTm,
    PrimerDimer,
    HighAnnealing,
    LowAnnealing,
    Fail,
}

#[derive(Eq, PartialEq, Hash, Debug, Clone)]
pub enum IndexResult {
    ThermoResult(ThermoResult),
    DigestError(DigestError),
    Pass(),
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
    pub fn _thermo_check(&mut self, dconf: &DigestConfig) -> &Option<IndexResult> {
        // Returns the None if valid primer
        let result = match &self.seq {
            Some(seq) => match thermo_check(&seq, dconf) {
                ThermoResult::Pass => IndexResult::Pass(),
                _ => IndexResult::ThermoResult(thermo_check(&seq, dconf)),
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
}

pub fn thermo_check(kmer: &[u8], dconf: &DigestConfig) -> ThermoResult {
    // Check GC content
    let gc = gc_content(kmer);
    if gc < dconf.primer_gc_min {
        return ThermoResult::LowGC;
    }
    if gc > dconf.primer_gc_max {
        return ThermoResult::HighGC;
    }

    // Check homopolymer
    if max_homopolymer(kmer) > dconf.max_homopolymers {
        return ThermoResult::Homopolymer;
    }

    // Check Tm if using TM Method

    match dconf.thermo_type {
        ThermoType::TM => {
            let tm = tm::oligo_tm_utf8(
                kmer,
                DNA_NM,
                K_MM,
                DIVALENT_CONC,
                DNTP_CONC,
                DMSO_CONC,
                DMSO_FACT,
                FORMAMIDE_CONC,
                tm::TmMethod::SantaLucia2004,
            );
            if tm >= dconf.primer_tm_max {
                return ThermoResult::HighTm;
            }
            if tm < dconf.primer_tm_min {
                return ThermoResult::LowTm;
            }
        }
        ThermoType::ANNEALING => {
            let prop = tm::oligo_annealing_utf8(
                kmer,
                DNA_NM,
                K_MM,
                DIVALENT_CONC,
                DNTP_CONC,
                DMSO_CONC,
                DMSO_FACT,
                FORMAMIDE_CONC,
                dconf.annealing_temp_c,
                tm::TmMethod::SantaLucia2004,
            );

            if prop > dconf.primer_annealing_prop.unwrap() + ANNEALING_DIFF {
                return ThermoResult::HighAnnealing;
            }
            if prop < dconf.primer_annealing_prop.unwrap() - ANNEALING_DIFF {
                return ThermoResult::LowAnnealing;
            }
        }
    }

    //TODO calc hairpin
    ThermoResult::Pass
}

fn process_seqs(
    seq_counts: HashMap<Result<Vec<u8>, DigestError>, f64>,
    dconf: &DigestConfig,
) -> Vec<DigestionKmers> {
    let mut digested: Vec<DigestionKmers> = Vec::with_capacity(seq_counts.len());
    for (k, v) in seq_counts.into_iter() {
        let dk = match k {
            Ok(s) => {
                // If sequence provided check thermo
                let mut dk = DigestionKmers::new(Some(s), None, v as f64);
                dk._thermo_check(dconf);
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
    digested.retain(|d| d.calc_freq(total) >= dconf.min_freq);

    digested
}

pub fn opto_binding(
    seq_counts: HashMap<Result<Vec<u8>, DigestError>, f64>,
    dconf: &DigestConfig,
) -> HashMap<Result<Vec<u8>, DigestError>, f64> {
    let mut return_map: HashMap<Result<Vec<u8>, DigestError>, f64> = HashMap::new();

    for (k, v) in seq_counts.into_iter() {
        let r = match k {
            Ok(seq) => {
                let truncated_seq = &seq[1..];

                match dconf.thermo_type {
                    ThermoType::TM => {
                        let target_tm = (dconf.primer_tm_max + dconf.primer_tm_min) / 2.0;

                        let full_tm_diff = (tm::oligo_tm_utf8(
                            &seq,
                            DNA_NM,
                            K_MM,
                            DIVALENT_CONC,
                            DNTP_CONC,
                            DMSO_CONC,
                            DMSO_FACT,
                            FORMAMIDE_CONC,
                            tm::TmMethod::SantaLucia2004,
                        ) - target_tm)
                            .abs();
                        let trunc_tm_diff = (tm::oligo_tm_utf8(
                            truncated_seq,
                            DNA_NM,
                            K_MM,
                            DIVALENT_CONC,
                            DNTP_CONC,
                            DMSO_CONC,
                            DMSO_FACT,
                            FORMAMIDE_CONC,
                            tm::TmMethod::SantaLucia2004,
                        ) - target_tm)
                            .abs();

                        // Return the seq with the closest tm
                        if full_tm_diff > trunc_tm_diff {
                            Ok(truncated_seq.to_vec())
                        } else {
                            Ok(seq)
                        }
                    }
                    ThermoType::ANNEALING => {
                        let full_annealing_diff = (tm::oligo_annealing_utf8(
                            &seq,
                            DNA_NM,
                            K_MM,
                            DIVALENT_CONC,
                            DNTP_CONC,
                            DMSO_CONC,
                            DMSO_FACT,
                            FORMAMIDE_CONC,
                            dconf.annealing_temp_c,
                            tm::TmMethod::SantaLucia2004,
                        ) - dconf.primer_annealing_prop.unwrap())
                        .abs();
                        let trunc_annealing_diff = (tm::oligo_annealing_utf8(
                            truncated_seq,
                            DNA_NM,
                            K_MM,
                            DIVALENT_CONC,
                            DNTP_CONC,
                            DMSO_CONC,
                            DMSO_FACT,
                            FORMAMIDE_CONC,
                            dconf.annealing_temp_c,
                            tm::TmMethod::SantaLucia2004,
                        ) - dconf.primer_annealing_prop.unwrap())
                        .abs();
                        // Return the seq with the closest tm
                        if full_annealing_diff > trunc_annealing_diff {
                            Ok(truncated_seq.to_vec())
                        } else {
                            Ok(seq)
                        }
                    }
                }
            }

            Err(e) => Err(e),
        };
        // Sum Sam seq
        let c = return_map.entry(r).or_insert(0.0);
        *c += v;
    }
    return return_map;
}

pub fn walk_right(
    seq: &[u8],
    l_index: usize,
    r_index: usize,
    kmer: tm::Oligo,
    dconf: &DigestConfig,
) -> Vec<Result<Vec<u8>, DigestError>> {
    // Check tm

    match dconf.thermo_type {
        ThermoType::TM => {
            if kmer.calc_tm(DNA_NM, K_MM, DIVALENT_CONC, DNTP_CONC, 0.0, 0.0, 0.8)
                >= dconf.primer_tm_min
            {
                return vec![Ok(kmer.seq)];
            }
        }
        ThermoType::ANNEALING => {
            if kmer.calc_annealing(
                dconf.annealing_temp_c,
                DNA_NM,
                K_MM,
                DIVALENT_CONC,
                DNTP_CONC,
            ) >= dconf.primer_annealing_prop.unwrap()
            {
                return vec![Ok(kmer.seq)];
            }
        }
    }

    // Check len
    if kmer.seq.len() > dconf.primer_len_max {
        return vec![Err(DigestError::ToLong)];
    }

    // Check if we've reached the end of the sequence
    if r_index < l_index || r_index - l_index >= dconf.max_walk {
        return vec![Err(DigestError::MaxWalk)];
    }
    // Check bounds
    if r_index == seq.len() {
        return vec![Err(DigestError::WalkedOutRight)];
    }

    let new_base = seq[r_index];

    match new_base {
        b' ' => return vec![Err(DigestError::EndOfSequence)],
        b'-' => return walk_right(seq, l_index, r_index + 1, kmer, dconf),
        b'N' => return vec![Err(DigestError::ContainsN)],
        _ => (),
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
        let new_results = walk_right(seq, l_index, r_index + 1, kmer_c, dconf);
        results.extend(new_results);
    }
    results
}

pub fn digest_r_to_count(
    seqs: &Vec<&[u8]>,
    index: usize,
    dconf: &DigestConfig,
) -> HashMap<Result<Vec<u8>, DigestError>, f64> {
    // Returns the kmer as read (left -> right) in the genomes

    let mut kmer_count: HashMap<Result<Vec<u8>, DigestError>, f64> = HashMap::new();

    // Check bounds
    if index > seqs[0].len() - dconf.primer_len_min {
        kmer_count.insert(Err(DigestError::WalkedOutRight), seqs.len() as f64);
        return kmer_count;
    }
    let rhs = index + dconf.primer_len_min;

    // For each sequence, digest at the index
    for seq in seqs.iter() {
        // Check for gap on set base
        match seq[index] {
            b'-' => {
                let c = kmer_count
                    .entry(Err(DigestError::GapOnSetBase))
                    .or_insert(0.0);
                *c += 1.0;
                continue;
            }
            b' ' => {
                let c = kmer_count
                    .entry(Err(DigestError::EndOfSequence))
                    .or_insert(0.0);
                *c += 1.0;
                continue;
            }
            _ => (),
        }

        // Create the kmer slice
        let mut kmer: Vec<u8> = Vec::with_capacity(dconf.primer_len_max);
        kmer.extend_from_slice(&seq[index..rhs]);

        // Remove gaps from the kmer
        kmer = kmer
            .into_iter()
            .filter(|b| *b != b'-' && *b != b' ')
            .collect();

        // Should not happen. As if set base is gap, it will be caught above
        if kmer.len() == 0 {
            let c = kmer_count
                .entry(Err(DigestError::NoValidPrimer))
                .or_insert(0.0);
            *c += 1.0;
            continue;
        }

        match check_kmer(&kmer) {
            KmerCheck::ATCGOnly => {
                // No ambiguous bases
                let results = walk_right(seq, index, rhs, tm::Oligo::new(kmer), dconf);
                let num_results = results.len() as f64;
                for r in results {
                    let count = kmer_count.entry(r).or_insert(0.0);
                    *count += 1.0 / num_results as f64;
                }
            }
            KmerCheck::ContainsAmbiguities => {
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
                            let results = walk_right(seq, index, rhs, tm::Oligo::new(ek), dconf);
                            let num_results = results.len() as f64;
                            for r in results {
                                let count = kmer_count.entry(r).or_insert(0.0);
                                *count += 1.0 / (num_results * num_kmers as f64);
                            }
                        }
                    }
                }
            }
            KmerCheck::ContainsInvalidBases => {
                let c = kmer_count
                    .entry(Err(DigestError::InvalidBase))
                    .or_insert(0.0);
                *c += 1.0;
            }
            KmerCheck::ContainsNs => {
                let c = kmer_count.entry(Err(DigestError::ContainsN)).or_insert(0.0);
                *c += 1.0;
            }
        }
    }
    // Rc the kmers
    let mut rc_kmer_count: HashMap<Result<Vec<u8>, DigestError>, f64> = HashMap::new();

    for (k, v) in kmer_count.into_iter() {
        let r = match k {
            Ok(seq) => Ok(reverse_complement(&seq)),
            Err(e) => Err(e),
        };
        rc_kmer_count.insert(r, v);
    }
    // Finally return opto kmers
    opto_binding(rc_kmer_count, dconf)
}

pub fn digest_r_at_index(
    seqs: &Vec<&[u8]>,
    index: usize,
    dconf: &DigestConfig,
) -> Result<RKmer, IndexResult> {
    let kmer_count: HashMap<Result<Vec<u8>, DigestError>, f64> =
        digest_r_to_count(seqs, index, dconf);

    // Process the results
    let dks = process_seqs(kmer_count, dconf);

    for dk in dks.iter() {
        match &dk.status {
            // Pass
            Some(IndexResult::ThermoResult(ThermoResult::Pass)) => {}
            Some(IndexResult::Pass()) => {}
            Some(IndexResult::DigestError(DigestError::EndOfSequence)) => {} // Ignore EOS errors
            // Maybe ignore N
            Some(IndexResult::DigestError(DigestError::ContainsN)) => {
                if dconf.ignore_n {
                    continue;
                } else {
                    return Err(IndexResult::DigestError(DigestError::ContainsN));
                }
            }
            // Fail
            Some(indexresult) => {
                return Err(indexresult.clone());
            }
            None => {
                return Err(IndexResult::DigestError(DigestError::NoValidPrimer));
            }
        }
    }
    // Filter EOS
    let dks: Vec<DigestionKmers> = dks.into_iter().filter(|dk| dk.seq.is_some()).collect();

    let seqs = dks.into_iter().map(|dk| dk.seq.unwrap()).collect();
    if primaldimer::do_pool_interact_u8(&seqs, &seqs, dconf.dimerscore) {
        return Err(IndexResult::ThermoResult(ThermoResult::PrimerDimer));
    }

    if seqs.len() == 0 {
        return Err(IndexResult::DigestError(DigestError::NoValidPrimer));
    }

    Ok(RKmer::new(seqs, index, None))
}

pub fn digest_r_dk(seq_array: &Vec<&[u8]>, dconf: &DigestConfig) -> Vec<Vec<DigestionKmers>> {
    let indexes: Vec<usize> = (0..seq_array[0].len() + 1).collect();

    // Check that all sequences are the same length
    for seq in seq_array.iter() {
        if seq.len() != seq_array[0].len() {
            panic!("Sequences are not the same length");
        }
    }

    let progress_bar = ProgressBar::new(indexes.len() as u64);
    progress_bar.set_message("rprimer digestion");
    progress_bar.set_style(
        ProgressStyle::default_bar()
            .template("{msg} [{elapsed}] {wide_bar:.cyan/blue} {pos:>7}/{len:7} {eta}")
            .unwrap(),
    );

    let results: Vec<Vec<DigestionKmers>> = indexes
        .par_iter()
        .progress_with(progress_bar)
        .map(|i| {
            let kmer_count: HashMap<Result<Vec<u8>, DigestError>, f64> =
                digest_r_to_count(seq_array, *i, dconf);

            // Process the results
            let dks = process_seqs(kmer_count, dconf);
            dks
        })
        .collect();

    results
}

pub fn digest_r_primer(
    seq_array: &Vec<&[u8]>,
    dconf: &DigestConfig,
    rindexes: Option<Vec<usize>>,
) -> Vec<Result<RKmer, IndexResult>> {
    // Allow custom indexes
    let indexes: Vec<usize> = match rindexes {
        // Could add checks
        Some(i) => i,
        None => (0..seq_array[0].len() + 1).collect(),
    };

    // Check that all sequences are the same length
    for seq in seq_array.iter() {
        if seq.len() != seq_array[0].len() {
            panic!("Sequences are not the same length");
        }
    }

    let progress_bar = ProgressBar::new(indexes.len() as u64);
    progress_bar.set_message("rprimer digestion");
    progress_bar.set_style(
        ProgressStyle::default_bar()
            .template("{msg} [{elapsed}] {wide_bar:.cyan/blue} {pos:>7}/{len:7} {eta}")
            .unwrap(),
    );

    let results: Vec<Result<RKmer, IndexResult>> = indexes
        .par_iter()
        .progress_with(progress_bar)
        .map(|i| digest_r_at_index(&seq_array, *i, dconf))
        .collect();

    results
}

pub fn walk_left(
    seq: &[u8],
    l_index: usize,
    r_index: usize,
    kmer: tm::Oligo,
    dconf: &DigestConfig,
) -> Vec<Result<Vec<u8>, DigestError>> {
    // Check tm
    match dconf.thermo_type {
        ThermoType::TM => {
            if kmer.calc_tm(DNA_NM, K_MM, DIVALENT_CONC, DNTP_CONC, 0.0, 0.0, 0.8)
                >= dconf.primer_tm_min
            {
                return vec![Ok(kmer.seq)];
            }
        }
        ThermoType::ANNEALING => {
            if kmer.calc_annealing(
                dconf.annealing_temp_c,
                DNA_NM,
                K_MM,
                DIVALENT_CONC,
                DNTP_CONC,
            ) >= dconf.primer_annealing_prop.unwrap()
            {
                return vec![Ok(kmer.seq)];
            }
        }
    }

    // Check len
    if kmer.seq.len() > dconf.primer_len_max {
        return vec![Err(DigestError::ToLong)];
    }

    // Check if we've reached the end of the sequence
    if l_index > r_index || r_index - l_index >= dconf.max_walk {
        return vec![Err(DigestError::MaxWalk)];
    }
    // Check bounds
    if l_index == 0 {
        return vec![Err(DigestError::WalkedOutLeft)];
    }

    let new_base = seq[l_index - 1];

    // If base is gap keep walking
    match new_base {
        b'-' => return walk_left(seq, l_index - 1, r_index, kmer, dconf),
        b' ' => return vec![Err(DigestError::EndOfSequence)],
        b'N' => return vec![Err(DigestError::ContainsN)],
        _ => (),
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
        kmer_c.add_base(complement_base(*base));
        let new_results = walk_left(seq, l_index - 1, r_index, kmer_c, dconf);
        results.extend(new_results);
    }
    results
}

pub fn digest_f_to_count(
    seqs: &Vec<&[u8]>,
    index: usize,
    dconf: &DigestConfig,
) -> HashMap<Result<Vec<u8>, DigestError>, f64> {
    let mut kmer_count: HashMap<Result<Vec<u8>, DigestError>, f64> = HashMap::new();

    // Check bounds
    if index < dconf.primer_len_min {
        kmer_count.insert(Err(DigestError::WalkedOutLeft), seqs.len() as f64);
        return kmer_count;
    }
    if index > seqs[0].len() {
        kmer_count.insert(Err(DigestError::WalkedOutRight), seqs.len() as f64);
        return kmer_count;
    }

    let lhs = index - dconf.primer_len_min;

    if index == 0 {
        kmer_count.insert(Err(DigestError::WalkedOutLeft), seqs.len() as f64);
        return kmer_count;
    }

    // For each sequence, digest at the index
    for seq in seqs.iter() {
        // Check for gap on set base
        match seq[index - 1] {
            b'-' => {
                let c = kmer_count
                    .entry(Err(DigestError::GapOnSetBase))
                    .or_insert(0.0);
                *c += 1.0;
                continue;
            }
            b' ' => {
                let c = kmer_count
                    .entry(Err(DigestError::EndOfSequence))
                    .or_insert(0.0);
                *c += 1.0;
                continue;
            }
            _ => (),
        }

        // Create the kmer slice
        let mut kmer: Vec<u8> = Vec::with_capacity(dconf.primer_len_max);
        kmer.extend_from_slice(&seq[lhs..index]);
        kmer = reverse_complement(&kmer); // Reverse is used as push is .O(1) and .insert(0) is O(n)

        // Remove gaps from the kmer
        kmer = kmer
            .into_iter()
            .filter(|b| *b != b'-' && *b != b' ')
            .collect();

        // Should not happen. As if set base is gap, it will be caught above
        if kmer.len() == 0 {
            let c = kmer_count
                .entry(Err(DigestError::NoValidPrimer))
                .or_insert(0.0);
            *c += 1.0;
            continue;
        }

        match check_kmer(&kmer) {
            KmerCheck::ATCGOnly => {
                // No ambiguous bases
                let results = walk_left(seq, lhs, index, tm::Oligo::new(kmer), dconf);
                let num_results = results.len() as f64;
                for r in results {
                    let count = kmer_count.entry(r).or_insert(0.0);
                    *count += 1.0 / num_results;
                }
            }
            KmerCheck::ContainsAmbiguities => {
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
                            let results = walk_left(seq, lhs, index, tm::Oligo::new(ek), dconf);
                            let num_results = results.len() as f64;
                            for r in results {
                                let count = kmer_count.entry(r).or_insert(0.0);
                                *count += 1.0 / (num_results * num_kmers as f64);
                            }
                        }
                    }
                }
            }
            KmerCheck::ContainsInvalidBases => {
                let c = kmer_count
                    .entry(Err(DigestError::InvalidBase))
                    .or_insert(0.0);
                *c += 1.0;
            }
            KmerCheck::ContainsNs => {
                let c = kmer_count.entry(Err(DigestError::ContainsN)).or_insert(0.0);
                *c += 1.0;
            }
        }
    }
    // Un-reverse the kmer
    let mut un_reversed: HashMap<Result<Vec<u8>, DigestError>, f64> = HashMap::new();
    for (k, v) in kmer_count.into_iter() {
        match k {
            Ok(mut kmer) => {
                kmer = reverse_complement(&kmer);
                un_reversed.insert(Ok(kmer), v);
            }
            Err(dr) => {
                un_reversed.insert(Err(dr), v);
            }
        }
    }

    // Finally return opto kmers
    opto_binding(un_reversed, dconf)
}

fn digest_f_at_index(
    seqs: &Vec<&[u8]>,
    index: usize,
    dconf: &DigestConfig,
) -> Result<FKmer, IndexResult> {
    let kmer_count: HashMap<Result<Vec<u8>, DigestError>, f64> =
        digest_f_to_count(seqs, index, dconf);

    // Process the results
    let dks: Vec<DigestionKmers> = process_seqs(kmer_count, dconf);

    // If any errors are found, return the error
    for dk in dks.iter() {
        match &dk.status {
            // Pass
            Some(IndexResult::ThermoResult(ThermoResult::Pass)) => {}
            Some(IndexResult::Pass()) => {}
            Some(IndexResult::DigestError(DigestError::EndOfSequence)) => {} // Ignore EOS errors
            // Ignore N
            Some(IndexResult::DigestError(DigestError::ContainsN)) => {
                if dconf.ignore_n {
                    continue;
                } else {
                    return Err(IndexResult::DigestError(DigestError::ContainsN));
                }
            }
            // Fail
            Some(indexresult) => {
                return Err(indexresult.clone());
            }
            None => {
                return Err(IndexResult::DigestError(DigestError::NoValidPrimer));
            }
        }
    }
    // Filter EOS
    let dks: Vec<DigestionKmers> = dks.into_iter().filter(|dk| dk.seq.is_some()).collect();

    let seqs = dks.into_iter().map(|dk| dk.seq.unwrap()).collect();
    if primaldimer::do_pool_interact_u8(&seqs, &seqs, dconf.dimerscore) {
        return Err(IndexResult::ThermoResult(ThermoResult::PrimerDimer));
    }
    if seqs.len() == 0 {
        return Err(IndexResult::DigestError(DigestError::NoValidPrimer));
    }
    // Create the FKmer
    Ok(FKmer::new(seqs, index, None))
}

pub fn digest_f_dk(seq_array: &Vec<&[u8]>, dconf: &DigestConfig) -> Vec<Vec<DigestionKmers>> {
    let indexes: Vec<usize> = (0..seq_array[0].len() + 1).collect();
    // Check that all sequences are the same length
    for seq in seq_array.iter() {
        if seq.len() != seq_array[0].len() {
            panic!("Sequences are not the same length");
        }
    }

    let progress_bar = ProgressBar::new(indexes.len() as u64);
    progress_bar.set_message("fprimer digestion");
    progress_bar.set_style(
        ProgressStyle::default_bar()
            .template("{msg} [{elapsed}] {wide_bar:.cyan/blue} {pos:>7}/{len:7} {eta}")
            .unwrap(),
    );

    let results: Vec<Vec<DigestionKmers>> = indexes
        .par_iter()
        .progress_with(progress_bar)
        .map(|i| {
            let fkc = digest_f_to_count(&seq_array, *i, dconf);
            process_seqs(fkc, dconf)
        })
        .collect();

    results
}

pub fn digest_f_primer(
    seq_array: &Vec<&[u8]>,
    dconf: &DigestConfig,
    findexes: Option<Vec<usize>>,
) -> Vec<Result<FKmer, IndexResult>> {
    // Allow custom indexes
    let indexes: Vec<usize> = match findexes {
        // Could add checks
        Some(i) => i,
        None => (0..seq_array[0].len() + 1).collect(),
    };
    // Check that all sequences are the same length
    for seq in seq_array.iter() {
        if seq.len() != seq_array[0].len() {
            panic!("Sequences are not the same length");
        }
    }

    let progress_bar = ProgressBar::new(indexes.len() as u64);
    progress_bar.set_message("fprimer digestion");
    progress_bar.set_style(
        ProgressStyle::default_bar()
            .template("{msg} [{elapsed}] {wide_bar:.cyan/blue} {pos:>7}/{len:7} {eta}")
            .unwrap(),
    );

    let results: Vec<Result<FKmer, IndexResult>> = indexes
        .par_iter()
        .progress_with(progress_bar)
        .map(|i| digest_f_at_index(&seq_array, *i, dconf))
        .collect();

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_digest_r_to_count() {
        let dconf = DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        );
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT".as_bytes()];

        let digested = digest_r_to_count(&seqs, 30, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 1);

        // Check sequence
        let exp_seq: Result<Vec<u8>, DigestError> =
            Ok("TCTACAAGAGATCGAAAGTTGGTTGGT".as_bytes().to_vec());

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
    fn test_digest_r_at_index_ps3() {
        let dconf = DigestConfig::new(
            None,
            None,
            None,
            Some(0.3),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        );
        let seqs =
            vec!["CCAATGGTGCAAAAGGTATAATCANTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT".as_bytes()];

        let digested = digest_r_at_index(&seqs, 25, &dconf);
        println!("{:?}", digested);
        // Check num of seqs
        assert_eq!(digested.is_ok(), true);

        // Check index
        let rkmer = digested.unwrap();
        assert_eq!(rkmer.start(), 25);

        // Match seq
        assert_eq!(
            rkmer.seqs(),
            vec![
                std::str::from_utf8("TGATTATACCTTTTGCACCATTGGACATTA".as_bytes())
                    .unwrap()
                    .to_string()
            ]
        );
    }

    #[test]
    fn test_digest_r_to_count_ambs() {
        let dconf = DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        );
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTRCGATCTCTTGTAGATCT".as_bytes()];

        let digested = digest_r_to_count(&seqs, 30, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 2);
        // Check count of ambiguous base
        for (_key, count) in digested.iter() {
            assert_eq!(count, &0.5, "Count: {}", count);
        }

        // Check sequence
        let seqs =
            vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTRCYGATCTCTTGTAGATCT".as_bytes()];
        let digested = digest_r_to_count(&seqs, 30, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 4);
        // Check count of ambiguous base
        for (_key, count) in digested.iter() {
            assert_eq!(count, &0.25, "Count: {}", count);
        }
    }

    #[test]
    fn test_digest_f_to_count() {
        let dconf = DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        );
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT".as_bytes()];

        let digested = digest_f_to_count(&seqs, 40, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 1);

        // Check sequence
        let exp_seq: Result<Vec<u8>, DigestError> =
            Ok("TTCCCAGGTAACAAACCAACCAAC".as_bytes().to_vec());

        assert_eq!(digested.contains_key(&exp_seq), true);
        // Check count
        assert_eq!(digested.get(&exp_seq).unwrap(), &1.0);
    }

    #[test]
    fn test_digest_f_to_count_ps3() {
        let dconf = DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        );
        let seqs =
            vec!["CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT".as_bytes()];

        let digested = digest_f_to_count(&seqs, 60, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 1);

        // Check sequence
        let exp_seq = Ok("GTCCAATGGTGCAAAAGGTATAATCATTAAT".as_bytes().to_vec());

        assert_eq!(digested.contains_key(&exp_seq), true);
        // Check count
        assert_eq!(digested.get(&exp_seq).unwrap(), &1.0);
    }

    #[test]
    fn test_digest_f_to_count_gap() {
        let dconf = DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        );
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAA-TTTCGATCTCTTGTAGATCT".as_bytes()];

        let digested = digest_f_to_count(&seqs, 40, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 1);
        // Check sequence
        assert_eq!(digested.contains_key(&Err(DigestError::GapOnSetBase)), true);
    }

    #[test]
    fn test_digest_f_to_count_ambs() {
        let dconf = DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        );
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAARTTTCGATCTCTTGTAGATCT".as_bytes()];

        let digested = digest_f_to_count(&seqs, 40, &dconf);

        // Check num of seqs
        assert_eq!(digested.len(), 2);
        // Check count of ambiguous base
        for (_key, count) in digested.iter() {
            assert_eq!(count, &0.5, "Count: {}", count);
        }
        // Check sequences
        let expected = vec![
            "CTTCCCAGGTAACAAACCAACCAAT".as_bytes().to_vec(),
            "TTCCCAGGTAACAAACCAACCAAC".as_bytes().to_vec(),
        ];
        for e in expected.into_iter() {
            assert!(digested.contains_key(&Ok(e)))
        }

        // Check sequence
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCRARTTTCGATCTCTTGTAGATCT".as_bytes()];
        let digested = digest_f_to_count(&seqs, 40, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 4);
        // Check count of ambiguous base
        for (_key, count) in digested.iter() {
            assert_eq!(count, &0.25, "Count: {}", count);
        }
    }

    #[test]
    fn test_digest_f_to_count_ambs_extend() {
        let dconf = DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        );
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCRARTTTCGATCTCTTGTAGATCT".as_bytes()];
        // Check that appending amb bases works as expected
        let digested = digest_f_to_count(&seqs, 60, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 4);

        // Check freqs
        for (_key, count) in digested.iter() {
            assert_eq!(*count, 1.0 / digested.len() as f64, "Count: {}", count);
        }
    }

    #[test]
    fn test_digest_f_to_count_wl() {
        let dconf = DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        );
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAA-TTTCGATCTCTTGTAGATCT".as_bytes()];

        let digested = digest_f_to_count(&seqs, 4, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 1);
        // Check sequence
        assert_eq!(
            digested.contains_key(&Err(DigestError::WalkedOutLeft)),
            true
        );
    }

    #[test]
    fn test_digest_f_to_count_contains_n() {
        let dconf = DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        );
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAANTTTCGATCTCTTGTAGATCT".as_bytes()];

        let digested: HashMap<Result<Vec<u8>, DigestError>, f64> =
            digest_f_to_count(&seqs, 40, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 1);
        // Check sequence
        assert_eq!(digested.contains_key(&Err(DigestError::ContainsN)), true);
    }
    #[test]
    fn test_digest_f_to_count_invalid_base() {
        let dconf = DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        );
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAA!TTTCGATCTCTTGTAGATCT".as_bytes()];

        let digested = digest_f_to_count(&seqs, 40, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 1);
        // Check sequence
        assert_eq!(digested.contains_key(&Err(DigestError::InvalidBase)), true);
    }

    #[test]
    fn test_opto_binding_tm_trunc() {
        let dconf = DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        );

        // See is to long. Therefore it should be truncated
        let seq1: Vec<u8> = b"CTTTTTAGTAGGTATAACCACAGCAGTTAAAAACTGATGCC".to_vec();
        let seq_count = 10.0;

        let mut count_map: HashMap<Result<Vec<u8>, DigestError>, f64> = HashMap::new();
        count_map.insert(Ok(seq1.to_vec()), seq_count);

        // Test TM
        let processed_map = opto_binding(count_map.clone(), &dconf);
        assert!(processed_map.len() == 1);
        // Check truncated from 5' end
        assert_eq!(
            *processed_map.keys().next().unwrap(),
            Ok(seq1[1..].to_vec())
        );
        assert_eq!(*processed_map.values().next().unwrap(), seq_count);
    }
    #[test]
    fn test_opto_binding_tm_no_trunc() {
        let dconf = DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        );

        // Seq is to short. Therefore it should not be truncated
        let seq1: Vec<u8> = b"CTTTTTAGTAGGTAT".to_vec();
        let seq_count = 10.0;

        let mut count_map: HashMap<Result<Vec<u8>, DigestError>, f64> = HashMap::new();
        count_map.insert(Ok(seq1.to_vec()), seq_count);

        // Test TM
        let processed_map = opto_binding(count_map.clone(), &dconf);
        assert!(processed_map.len() == 1);
        // Check not truncated
        assert_eq!(*processed_map.keys().next().unwrap(), Ok(seq1));
        assert_eq!(*processed_map.values().next().unwrap(), seq_count);
    }
    #[test]
    fn test_opto_binding_an_trunc() {
        let mut dconf = DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        );
        dconf.thermo_type = ThermoType::ANNEALING;
        dconf.primer_annealing_prop = Some(10.0);

        // See is to long. Therefore it should be truncated
        let seq1: Vec<u8> = b"CTTTTTAGTAGGTATAACCACAGCAGTTAAAAACTGATGCC".to_vec();
        let seq_count = 10.0;

        let mut count_map: HashMap<Result<Vec<u8>, DigestError>, f64> = HashMap::new();
        count_map.insert(Ok(seq1.to_vec()), seq_count);

        // Test annealing
        let processed_map = opto_binding(count_map.clone(), &dconf);
        assert!(processed_map.len() == 1);
        // Check truncated from 5' end
        assert_eq!(
            *processed_map.keys().next().unwrap(),
            Ok(seq1[1..].to_vec())
        );
        assert_eq!(*processed_map.values().next().unwrap(), seq_count);
    }
    #[test]
    fn test_opto_binding_an() {
        let mut dconf = DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        );
        dconf.thermo_type = ThermoType::ANNEALING;
        dconf.primer_annealing_prop = Some(10.0);

        // See is to long. Therefore it should be truncated
        let seq1: Vec<u8> = b"CTTTTTAGTAGGTAT".to_vec();
        let seq_count = 10.0;

        let mut count_map: HashMap<Result<Vec<u8>, DigestError>, f64> = HashMap::new();
        count_map.insert(Ok(seq1.to_vec()), seq_count);

        // Test annealing
        let processed_map = opto_binding(count_map.clone(), &dconf);
        assert!(processed_map.len() == 1);
        // Check truncated from 5' end
        assert_eq!(*processed_map.keys().next().unwrap(), Ok(seq1));
        assert_eq!(*processed_map.values().next().unwrap(), seq_count);
    }
}
