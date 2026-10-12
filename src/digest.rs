use crate::config::{DigestConfig, ThermoType};
use crate::kmer::{FKmer, KmerReport, RKmer, LOW_FREQ};
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

/// thal conditions for the hairpin check: the digest's salt and DNA
/// concentrations, matching primalscheme3's config (37 °C, max_loop 30).
fn thal_args() -> primer3_thal_rs::ThalArgs {
    primer3_thal_rs::ThalArgs {
        mv: K_MM,
        dv: DIVALENT_CONC,
        dntp: DNTP_CONC,
        dna_conc: DNA_NM,
        ..Default::default()
    }
}

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
    NotChecked(),
}

/// Log status of a kmer kept by the digest.
pub static PASSED: IndexResult = IndexResult::Pass();
/// Log status of a kmer returned by `return_all` that a normal digest rejects.
pub static RETURNED_FAILING: IndexResult = IndexResult::ThermoResult(ThermoResult::Fail);

pub struct DigestionKmers {
    pub seq: Option<Vec<u8>>,
    pub status: Option<IndexResult>,
    pub count: f64,
    /// `return_all` only: every check the sequence fails.
    pub failures: Vec<ThermoResult>,
    /// `return_all` only: below `min_freq` (a normal digest drops these).
    pub low_freq: bool,
}
impl DigestionKmers {
    pub fn new(seq: Option<Vec<u8>>, status: Option<IndexResult>, count: f64) -> DigestionKmers {
        DigestionKmers {
            seq: seq,
            status: status,
            count: count,
            failures: Vec::new(),
            low_freq: false,
        }
    }
    pub fn calc_freq(&self, total: f64) -> f64 {
        self.count / total
    }
    pub fn _thermo_check(&mut self, dconf: &DigestConfig) -> &Option<IndexResult> {
        // If not check skip all checks
        if !dconf.thermo_check {
            self.status = Some(IndexResult::NotChecked());
            return &self.status;
        }
        // Returns the None if valid primer
        let result = match &self.seq {
            // return_all: record every failure; the status is still the first.
            Some(seq) if dconf.return_all => {
                self.failures = thermo_failures(seq, dconf, true);
                match self.failures.first() {
                    None => IndexResult::Pass(),
                    Some(failed) => IndexResult::ThermoResult(failed.clone()),
                }
            }
            Some(seq) => match thermo_check(&seq, dconf) {
                ThermoResult::Pass => IndexResult::Pass(),
                failed => IndexResult::ThermoResult(failed),
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

/// The first check `kmer` fails, or `Pass`.
pub fn thermo_check(kmer: &[u8], dconf: &DigestConfig) -> ThermoResult {
    thermo_failures(kmer, dconf, false)
        .into_iter()
        .next()
        .unwrap_or(ThermoResult::Pass)
}

/// Every check `kmer` fails, in check order (GC, homopolymer, Tm/annealing,
/// hairpin). With `all` false, stops at the first failure.
pub fn thermo_failures(kmer: &[u8], dconf: &DigestConfig, all: bool) -> Vec<ThermoResult> {
    let mut failures = Vec::new();
    macro_rules! fail {
        ($result:expr) => {{
            failures.push($result);
            if !all {
                return failures;
            }
        }};
    }

    // Check GC content
    let gc = gc_content(kmer);
    if gc < dconf.primer_gc_min {
        fail!(ThermoResult::LowGC);
    }
    if gc > dconf.primer_gc_max {
        fail!(ThermoResult::HighGC);
    }

    // Check homopolymer
    if max_homopolymer(kmer) > dconf.max_homopolymers {
        fail!(ThermoResult::Homopolymer);
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
                fail!(ThermoResult::HighTm);
            }
            if tm < dconf.primer_tm_min {
                fail!(ThermoResult::LowTm);
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

            if prop > dconf.primer_annealing_prop.unwrap() + dconf.annealing_diff {
                fail!(ThermoResult::HighAnnealing);
            }
            if prop < dconf.primer_annealing_prop.unwrap() - dconf.annealing_diff {
                fail!(ThermoResult::LowAnnealing);
            }
        }
    }

    // Hairpin: fail if its Tm exceeds the limit (as primalscheme3's forms_hairpin)
    match primer3_thal_rs::calc_hairpin(kmer, &thal_args(), false) {
        Ok(hairpin) if hairpin.tm > dconf.primer_hairpin_th_max => fail!(ThermoResult::Hairpin),
        Ok(_) => {}
        // Only for > 60 nt, which primer_len_max prevents
        Err(_) => fail!(ThermoResult::Fail),
    }

    failures
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

    // Apply a freq filter (return_all marks rather than drops)
    let total: f64 = digested.iter().map(|d| d.count).sum();
    if dconf.return_all {
        for d in digested.iter_mut() {
            d.low_freq = d.calc_freq(total) < dconf.min_freq;
        }
    } else {
        digested.retain(|d| d.calc_freq(total) >= dconf.min_freq);
    }

    digested
}

/// Name of a status as reported to Python, e.g. `HighGC` or `ContainsN`.
fn status_name(status: &IndexResult) -> String {
    match status {
        IndexResult::ThermoResult(t) => format!("{t:?}"),
        IndexResult::DigestError(e) => format!("{e:?}"),
        IndexResult::Pass() => "Pass".to_string(),
        IndexResult::NotChecked() => "NotChecked".to_string(),
    }
}

/// Collect a `return_all` report for one index: every sequence with its
/// failures, the kmer's dimer status, and positions without a sequence.
/// Errors only if no sequence was produced at all.
fn report_all(dks: Vec<DigestionKmers>, dconf: &DigestConfig) -> Result<KmerReport, IndexResult> {
    let mut report = KmerReport {
        seqs: Vec::new(),
        counts: Vec::new(),
        statuses: Vec::new(),
        kmer_statuses: Vec::new(),
        errors: Vec::new(),
        dimers: Vec::new(),
    };
    let mut first_error = None;
    for dk in dks {
        let tag = |mut tags: Vec<String>| {
            if dk.low_freq {
                tags.push(LOW_FREQ.to_string());
            }
            tags
        };
        match (&dk.seq, &dk.status) {
            (Some(_), _) => {
                let failures = dk.failures.iter().map(|f| format!("{f:?}")).collect();
                report.statuses.push(tag(failures));
                report.counts.push(dk.count);
                report.seqs.push(dk.seq.unwrap());
            }
            // Ignored by a normal digest too
            (None, Some(IndexResult::DigestError(DigestError::EndOfSequence))) => {}
            (None, Some(IndexResult::DigestError(DigestError::ContainsN))) if dconf.ignore_n => {}
            (None, status) => {
                let status = status
                    .clone()
                    .unwrap_or(IndexResult::DigestError(DigestError::NoValidPrimer));
                report.errors.push((tag(vec![status_name(&status)]), dk.count));
                first_error.get_or_insert(status);
            }
        }
    }
    if report.seqs.is_empty() {
        return Err(first_error.unwrap_or(IndexResult::DigestError(DigestError::NoValidPrimer)));
    }
    // Dimer check over the sequences a normal digest would keep, recording
    // every interacting pair (as indices into report.seqs)
    if dconf.thermo_check {
        let kept: Vec<usize> = (0..report.seqs.len())
            .filter(|&i| !report.statuses[i].iter().any(|t| t == LOW_FREQ))
            .collect();
        let kept_seqs: Vec<Vec<u8>> = kept.iter().map(|&i| report.seqs[i].clone()).collect();
        report.dimers = primaldimer::pool_dimer_pairs(&kept_seqs, dconf.dimerscore)
            .into_iter()
            .map(|(a, b)| (kept[a], kept[b]))
            .collect();
        if !report.dimers.is_empty() {
            report.kmer_statuses.push(format!("{:?}", ThermoResult::PrimerDimer));
        }
    }
    Ok(report)
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

pub fn dedup_substr(
    mut seq_counts: HashMap<Result<Vec<u8>, DigestError>, f64>,
) -> HashMap<Result<Vec<u8>, DigestError>, f64> {
    // Check if the truncated version of the sequence is present.
    let keys: Vec<Vec<u8>> = seq_counts
        .keys()
        .filter_map(|k| k.as_ref().ok().cloned())
        .collect();

    for seq in keys {
        if seq.len() < 2 {
            continue;
        }
        // Check if the [1:] slice exists
        let sub = seq[1..].to_vec();
        let sub_key = Ok(sub);

        if seq_counts.contains_key(&sub_key) {
            // If it does, add the original count to sub sequence
            if let Some(val) = seq_counts.remove(&Ok(seq)) {
                if let Some(count) = seq_counts.get_mut(&sub_key) {
                    *count += val;
                }
            }
        }
    }
    return seq_counts;
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
    let opto_seq_count = opto_binding(rc_kmer_count, dconf);
    // Deduplicate
    dedup_substr(opto_seq_count)
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
    if dconf.return_all {
        return report_all(dks, dconf).map(|report| RKmer::from_report(report, index));
    }

    for dk in dks.iter() {
        match &dk.status {
            // Pass
            Some(IndexResult::ThermoResult(ThermoResult::Pass)) => {}
            Some(IndexResult::Pass()) => {}
            Some(IndexResult::DigestError(DigestError::EndOfSequence)) => {} // Ignore EOS errors
            // Maybe allow non-thermo checked
            Some(IndexResult::NotChecked()) => {
                // If thermo check is not enabled, we can ignore this
                if dconf.thermo_check {
                    return Err(IndexResult::DigestError(DigestError::NoValidPrimer));
                }
            }
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
    let counts = dks.iter().map(|dk| dk.count).collect::<Vec<f64>>();
    let seqs = dks.into_iter().map(|dk| dk.seq.unwrap()).collect();

    // Skip dimer checks is no-thermo
    if dconf.thermo_check && primaldimer::do_pool_interact_u8(&seqs, &seqs, dconf.dimerscore) {
        return Err(IndexResult::ThermoResult(ThermoResult::PrimerDimer));
    }

    if seqs.len() == 0 {
        return Err(IndexResult::DigestError(DigestError::NoValidPrimer));
    }

    Ok(RKmer::new(seqs, index, Some(counts)))
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

    let opto_seq_count = opto_binding(un_reversed, dconf);
    // Deduplicate
    dedup_substr(opto_seq_count)
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
    if dconf.return_all {
        return report_all(dks, dconf).map(|report| FKmer::from_report(report, index));
    }

    // If any errors are found, return the error
    for dk in dks.iter() {
        match &dk.status {
            // Pass
            Some(IndexResult::ThermoResult(ThermoResult::Pass)) => {}
            Some(IndexResult::Pass()) => {}
            Some(IndexResult::DigestError(DigestError::EndOfSequence)) => {} // Ignore EOS errors
            // Maybe allow non-thermo checked
            Some(IndexResult::NotChecked()) => {
                // If thermo check is not enabled, we can ignore this
                if dconf.thermo_check {
                    return Err(IndexResult::DigestError(DigestError::NoValidPrimer));
                }
            }
            // Maybe Ignore N
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

    let counts = dks.iter().map(|dk| dk.count).collect::<Vec<f64>>();
    let seqs = dks.into_iter().map(|dk| dk.seq.unwrap()).collect();
    // Skip dimer checks is no-thermo
    if dconf.thermo_check && primaldimer::do_pool_interact_u8(&seqs, &seqs, dconf.dimerscore) {
        return Err(IndexResult::ThermoResult(ThermoResult::PrimerDimer));
    }
    if seqs.len() == 0 {
        return Err(IndexResult::DigestError(DigestError::NoValidPrimer));
    }

    Ok(FKmer::new(seqs, index, Some(counts)))
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

    /// A config whose GC, Tm and homopolymer limits pass anything, so only
    /// the hairpin check can fail.
    fn hairpin_only_config(primer_hairpin_th_max: f64) -> DigestConfig {
        DigestConfig::new(
            None,
            None,
            Some(1.0),
            Some(0.0),
            Some(1000.0),
            Some(-1000.0),
            None,
            None,
            None,
            None,
            None,
            Some(100),
            None,
            None,
            None,
            None,
            Some(primer_hairpin_th_max),
            None,
        )
    }

    #[test]
    fn test_hairpin_tm_matches_primer3_py() {
        // primer3.calc_hairpin(seq, mv_conc=100, dv_conc=2, dntp_conc=0.8, dna_conc=15)
        let cases = [
            ("CCCCCATCCGATCAGGGGG", 67.96550347919458),
            ("ACGTACGTTTTTACGTACGT", 69.16336781546374),
            ("GTAAAACGACGGCCAGTGAAT", 0.0),
        ];
        for (seq, expected) in cases {
            let res = primer3_thal_rs::calc_hairpin(seq, &thal_args(), false).unwrap();
            assert!((res.tm - expected).abs() < 1e-6, "{seq}: {} vs {expected}", res.tm);
        }
    }

    #[test]
    fn test_thermo_check_hairpin() {
        let seq = b"CCCCCATCCGATCAGGGGG"; // hairpin Tm 67.97 °C
        assert_eq!(
            thermo_check(seq, &hairpin_only_config(55.0)),
            ThermoResult::Hairpin
        );
        assert_eq!(
            thermo_check(seq, &hairpin_only_config(100.0)),
            ThermoResult::Pass
        );
        // Boundaries: fails only when the Tm (67.9655 °C) is strictly above
        // the limit, as primalscheme3's `tm > primer_hairpin_th_max`.
        let tm = primer3_thal_rs::calc_hairpin(seq, &thal_args(), false).unwrap().tm;
        assert_eq!(thermo_check(seq, &hairpin_only_config(tm - 0.01)), ThermoResult::Hairpin);
        assert_eq!(thermo_check(seq, &hairpin_only_config(tm)), ThermoResult::Pass);
        assert_eq!(thermo_check(seq, &hairpin_only_config(tm + 0.01)), ThermoResult::Pass);
        // No hairpin at all (Tm 0) passes any positive limit.
        assert_eq!(
            thermo_check(b"GTAAAACGACGGCCAGTGAAT", &hairpin_only_config(1.0)),
            ThermoResult::Pass
        );
    }

    #[test]
    fn test_thermo_failures_reports_every_check() {
        // 74% GC and a 68 °C hairpin
        let seq = b"CCCCCATCCGATCAGGGGG";
        let mut dconf = hairpin_only_config(55.0);
        dconf.primer_gc_max = 0.5;
        assert_eq!(
            thermo_failures(seq, &dconf, true),
            vec![ThermoResult::HighGC, ThermoResult::Hairpin]
        );
        // The normal check still stops at the first failure
        assert_eq!(thermo_check(seq, &dconf), ThermoResult::HighGC);
        assert_eq!(thermo_failures(seq, &dconf, false), vec![ThermoResult::HighGC]);
        assert!(thermo_failures(b"GTAAAACGACGGCCAGTGAAT", &hairpin_only_config(55.0), true).is_empty());
    }

    /// Three aligned variants of a random genome with SNPs, an N run and an
    /// ambiguity code, so indexes hit thermo failures, digest errors and
    /// low-frequency sequences.
    fn variant_alignment() -> Vec<Vec<u8>> {
        let mut state: u64 = 42;
        let mut next = move || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as usize
        };
        let genome: Vec<u8> = (0..1500).map(|_| b"ACGT"[next() % 4]).collect();
        let mut seqs = vec![genome.clone(); 8];
        for seq in seqs.iter_mut().skip(1) {
            for _ in 0..12 {
                let i = next() % seq.len();
                seq[i] = b"ACGT"[next() % 4];
            }
        }
        seqs[3][700..706].copy_from_slice(b"NNNNNN");
        seqs[5][900] = b'R';
        seqs
    }

    #[test]
    fn test_return_all_matches_normal_digest() {
        // Default config, with some low-frequency sequences
        check_return_all_parity(|c| c.min_freq = 0.2);
    }

    #[test]
    fn test_return_all_parity_ignore_n() {
        check_return_all_parity(|c| c.ignore_n = true);
    }

    #[test]
    fn test_return_all_parity_annealing() {
        check_return_all_parity(|c| {
            c.thermo_type = ThermoType::ANNEALING;
            c.primer_annealing_prop = Some(30.0);
            c.annealing_diff = 20.0;
        });
    }

    #[test]
    fn test_return_all_parity_thermo_off() {
        check_return_all_parity(|c| c.thermo_check = false);
    }

    /// With `return_all`, every index agrees with a normal digest: `passed()`
    /// exactly when the normal digest keeps the kmer, with the same
    /// sequences and counts once low-frequency sequences are dropped.
    fn check_return_all_parity(tweak: fn(&mut DigestConfig)) {
        let seqs = variant_alignment();
        let seq_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let mut normal = DigestConfig::create_default();
        tweak(&mut normal);
        let mut all = DigestConfig::create_default();
        tweak(&mut all);
        all.return_all = true;

        let (mut kept, mut returned, mut flagged, mut low) = (0, 0, 0, 0);
        let mut dimer_kmers = 0;
        for index in 0..=seqs[0].len() {
            for forward in [true, false] {
                let (n, a) = if forward {
                    (
                        digest_f_at_index(&seq_refs, index, &normal).map(|k| (k.seqs(), k.counts().clone())),
                        digest_f_at_index(&seq_refs, index, &all).map(|k| {
                            assert_eq!(
                                !k.dimers().is_empty(),
                                k.kmer_statuses().contains(&"PrimerDimer".to_string()),
                                "index {index}"
                            );
                            dimer_kmers += !k.dimers().is_empty() as usize;
                            (k.seqs(), k.counts().clone(), k.statuses().unwrap(), k.passed())
                        }),
                    )
                } else {
                    (
                        digest_r_at_index(&seq_refs, index, &normal).map(|k| (k.seqs(), k.counts().clone())),
                        digest_r_at_index(&seq_refs, index, &all).map(|k| {
                            assert_eq!(
                                !k.dimers().is_empty(),
                                k.kmer_statuses().contains(&"PrimerDimer".to_string()),
                                "index {index}"
                            );
                            dimer_kmers += !k.dimers().is_empty() as usize;
                            (k.seqs(), k.counts().clone(), k.statuses().unwrap(), k.passed())
                        }),
                    )
                };
                if let Ok((_, _, statuses, _)) = &a {
                    returned += 1;
                    low += statuses.iter().filter(|t| t.iter().any(|x| x == LOW_FREQ)).count();
                    if statuses.iter().any(|t| t.len() > 1) {
                        flagged += 1;
                    }
                }
                match (n, a) {
                    (Ok((n_seqs, n_counts)), Ok((a_seqs, a_counts, statuses, passed))) => {
                        kept += 1;
                        assert!(passed, "index {index}: normal kept, return_all failed");
                        // Same kmer once low-frequency sequences are dropped
                        let keep: Vec<bool> = statuses.iter().map(|t| !t.iter().any(|x| x == LOW_FREQ)).collect();
                        let a_seqs: Vec<String> = a_seqs.into_iter().zip(&keep).filter(|(_, k)| **k).map(|(s, _)| s).collect();
                        let a_counts: Vec<f64> = a_counts.unwrap().into_iter().zip(&keep).filter(|(_, k)| **k).map(|(c, _)| c).collect();
                        assert_eq!(n_seqs, a_seqs, "index {index}");
                        assert_eq!(n_counts.unwrap(), a_counts, "index {index}");
                    }
                    (Ok(_), Err(e)) => panic!("index {index}: normal kept, return_all rejected {e:?}"),
                    (Err(_), Ok((_, _, _, passed))) => assert!(!passed, "index {index}: normal rejected, return_all passed"),
                    (Err(_), Err(_)) => {}
                }
            }
        }
        eprintln!("return_all parity: {kept} kept normally, {returned} returned, {flagged} with a multi-failure sequence, {low} low-frequency sequences, {dimer_kmers} with dimer pairs");
        assert!(kept > 0 && returned >= kept, "kept {kept}, returned {returned}");
        if all.thermo_check {
            assert!(returned > kept, "nothing was rejected");
            assert!(flagged > 0, "no sequence failed more than one check");
            assert!(dimer_kmers > 0, "no dimer kmers");
        }
        if all.min_freq > 0.0 {
            assert!(low > 0, "no low-frequency sequences");
        }
    }

    /// A config that passes any sequence: every check is wide open.
    fn permissive_config() -> DigestConfig {
        let mut dconf = hairpin_only_config(1000.0);
        dconf.max_homopolymers = 1000;
        dconf
    }

    #[test]
    fn test_thermo_failures_each_check_in_isolation() {
        let plain = b"ACGTACGATCGATGCATGCA"; // 50% GC, no long homopolymer
        let homopolymer = b"ACGTAAAAAAGATGCATGCA";
        let hairpin = b"CCCCCATCCGATCAGGGGG"; // hairpin Tm 68 °C
        type Tweak = fn(&mut DigestConfig);
        let cases: [(&str, &[u8], Tweak, ThermoResult); 9] = [
            ("low gc", plain, |c| c.primer_gc_min = 0.9, ThermoResult::LowGC),
            ("high gc", plain, |c| c.primer_gc_max = 0.1, ThermoResult::HighGC),
            ("homopolymer", homopolymer, |c| c.max_homopolymers = 3, ThermoResult::Homopolymer),
            ("high tm", plain, |c| c.primer_tm_max = 0.0, ThermoResult::HighTm),
            ("low tm", plain, |c| c.primer_tm_min = 1000.0, ThermoResult::LowTm),
            ("high annealing", plain, |c| {
                c.thermo_type = ThermoType::ANNEALING;
                c.primer_annealing_prop = Some(-1000.0);
                c.annealing_diff = 0.0;
            }, ThermoResult::HighAnnealing),
            ("low annealing", plain, |c| {
                c.thermo_type = ThermoType::ANNEALING;
                c.primer_annealing_prop = Some(1000.0);
                c.annealing_diff = 0.0;
            }, ThermoResult::LowAnnealing),
            ("hairpin", hairpin, |c| c.primer_hairpin_th_max = 55.0, ThermoResult::Hairpin),
            // Hairpin calculation refuses > 60 nt; primer_len_max normally prevents it
            ("too long", &[b'A'; 61], |_| {}, ThermoResult::Fail),
        ];
        for (name, seq, tweak, expected) in cases {
            let mut dconf = permissive_config();
            tweak(&mut dconf);
            assert_eq!(thermo_failures(seq, &dconf, true), vec![expected.clone()], "{name}");
            assert_eq!(thermo_failures(seq, &dconf, false), vec![expected.clone()], "{name}");
            assert_eq!(thermo_check(seq, &dconf), expected, "{name}");
        }
        // Nothing fails under the permissive config
        for seq in [&plain[..], &homopolymer[..], &hairpin[..]] {
            assert!(thermo_failures(seq, &permissive_config(), true).is_empty());
            assert_eq!(thermo_check(seq, &permissive_config()), ThermoResult::Pass);
        }
        // Annealing checks pass within the window
        let mut dconf = permissive_config();
        dconf.thermo_type = ThermoType::ANNEALING;
        dconf.primer_annealing_prop = Some(50.0);
        dconf.annealing_diff = 1000.0;
        assert!(thermo_failures(plain, &dconf, true).is_empty());
    }

    #[test]
    fn test_thermo_failures_all_in_check_order() {
        // Fails GC, homopolymer, Tm and hairpin at once
        let seq = b"CCCCCATCCGATCAGGGGG";
        let mut dconf = permissive_config();
        dconf.primer_gc_max = 0.1;
        dconf.max_homopolymers = 3;
        dconf.primer_tm_min = 1000.0;
        dconf.primer_hairpin_th_max = 55.0;
        let all = vec![
            ThermoResult::HighGC,
            ThermoResult::Homopolymer,
            ThermoResult::LowTm,
            ThermoResult::Hairpin,
        ];
        assert_eq!(thermo_failures(seq, &dconf, true), all);
        assert_eq!(thermo_failures(seq, &dconf, false), vec![ThermoResult::HighGC]);
        assert_eq!(thermo_check(seq, &dconf), ThermoResult::HighGC);
    }

    #[test]
    fn test_digestion_kmer_thermo_check_modes() {
        let seq = b"CCCCCATCCGATCAGGGGG".to_vec();
        let mut dconf = permissive_config();
        dconf.primer_gc_max = 0.1;
        dconf.primer_hairpin_th_max = 55.0;

        // Normal: first failure only, no failure list
        let mut dk = DigestionKmers::new(Some(seq.clone()), None, 1.0);
        dk._thermo_check(&dconf);
        assert_eq!(dk.status, Some(IndexResult::ThermoResult(ThermoResult::HighGC)));
        assert!(dk.failures.is_empty());

        // return_all: status is the first failure, failures has them all
        dconf.return_all = true;
        let mut dk = DigestionKmers::new(Some(seq.clone()), None, 1.0);
        dk._thermo_check(&dconf);
        assert_eq!(dk.status, Some(IndexResult::ThermoResult(ThermoResult::HighGC)));
        assert_eq!(dk.failures, vec![ThermoResult::HighGC, ThermoResult::Hairpin]);

        // A passing sequence
        let mut dk = DigestionKmers::new(Some(seq.clone()), None, 1.0);
        dk._thermo_check(&permissive_config());
        assert_eq!(dk.status, Some(IndexResult::Pass()));

        // Thermo checks off
        let mut off = permissive_config();
        off.thermo_check = false;
        off.return_all = true;
        let mut dk = DigestionKmers::new(Some(seq.clone()), None, 1.0);
        dk._thermo_check(&off);
        assert_eq!(dk.status, Some(IndexResult::NotChecked()));
        assert!(dk.failures.is_empty());

        // No sequence
        let mut dk = DigestionKmers::new(None, None, 1.0);
        dk._thermo_check(&dconf);
        assert_eq!(dk.status, Some(IndexResult::DigestError(DigestError::NoValidPrimer)));
    }

    #[test]
    fn test_process_seqs_min_freq() {
        let counts = || {
            HashMap::from([
                (Ok(b"ACGTACGATCGATGCATGCA".to_vec()), 9.0),
                (Ok(b"ACGTACGATCGATGCATGCT".to_vec()), 1.0),
                (Err(DigestError::ContainsN), 0.5),
            ])
        };
        let mut dconf = permissive_config();
        dconf.min_freq = 0.2;
        // Normal: below min_freq is dropped
        let kept = process_seqs(counts(), &dconf);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].count, 9.0);
        // return_all: kept and marked
        dconf.return_all = true;
        let all = process_seqs(counts(), &dconf);
        assert_eq!(all.len(), 3);
        for dk in &all {
            assert_eq!(dk.low_freq, dk.count < 2.1, "count {}", dk.count);
        }
    }

    fn dk(seq: Option<&[u8]>, status: Option<IndexResult>, count: f64, failures: Vec<ThermoResult>, low_freq: bool) -> DigestionKmers {
        let mut d = DigestionKmers::new(seq.map(|s| s.to_vec()), status, count);
        d.failures = failures;
        d.low_freq = low_freq;
        d
    }

    #[test]
    fn test_report_all_branches() {
        let mut dconf = permissive_config();
        dconf.return_all = true;
        let err = |e| Some(IndexResult::DigestError(e));
        let dks = vec![
            dk(Some(b"ACGTACGATCGATGCATGCA"), Some(IndexResult::Pass()), 5.0, vec![], false),
            dk(Some(b"CCCCCATCCGATCAGGGGG"), Some(IndexResult::ThermoResult(ThermoResult::HighGC)), 3.0,
               vec![ThermoResult::HighGC, ThermoResult::Hairpin], false),
            dk(Some(b"ACGTACGATCGATGCATGCT"), Some(IndexResult::Pass()), 1.0, vec![], true),
            dk(None, err(DigestError::EndOfSequence), 1.0, vec![], false), // skipped, as normal
            dk(None, err(DigestError::WalkedOutLeft), 2.0, vec![], false),
            dk(None, err(DigestError::ContainsN), 1.0, vec![], true),
            dk(None, None, 1.0, vec![], false), // no status: NoValidPrimer
        ];
        let report = report_all(dks, &dconf).unwrap();
        assert_eq!(report.counts, vec![5.0, 3.0, 1.0]);
        assert_eq!(
            report.statuses,
            vec![vec![], vec!["HighGC".to_string(), "Hairpin".to_string()], vec!["LowFreq".to_string()]]
        );
        assert_eq!(
            report.errors,
            vec![
                (vec!["WalkedOutLeft".to_string()], 2.0),
                (vec!["ContainsN".to_string(), "LowFreq".to_string()], 1.0),
                (vec!["NoValidPrimer".to_string()], 1.0),
            ]
        );

        // ContainsN is skipped when N is ignored
        dconf.ignore_n = true;
        let report = report_all(
            vec![
                dk(Some(b"ACGTACGATCGATGCATGCA"), Some(IndexResult::Pass()), 1.0, vec![], false),
                dk(None, err(DigestError::ContainsN), 1.0, vec![], false),
            ],
            &dconf,
        )
        .unwrap();
        assert!(report.errors.is_empty());

        // No sequence at all: the first error, as a normal digest reports
        let none = report_all(
            vec![
                dk(None, err(DigestError::EndOfSequence), 1.0, vec![], false),
                dk(None, err(DigestError::MaxWalk), 1.0, vec![], false),
                dk(None, err(DigestError::WalkedOutLeft), 1.0, vec![], false),
            ],
            &dconf,
        );
        assert_eq!(none.err(), Some(IndexResult::DigestError(DigestError::MaxWalk)));
        assert_eq!(
            report_all(vec![], &dconf).err(),
            Some(IndexResult::DigestError(DigestError::NoValidPrimer))
        );
    }

    #[test]
    fn test_report_all_dimer_check() {
        let selfcomp = b"AAAAAAAAAAGCGCGCGCGC"; // 3' self-complementary
        let plain = b"GTAAAACGACGGCCAGTGAAT"; // M13: dimerises with neither
        let mut dconf = permissive_config();
        dconf.return_all = true;
        let run = |dconf: &DigestConfig, low: bool| {
            report_all(
                vec![
                    dk(Some(plain), Some(IndexResult::Pass()), 1.0, vec![], false),
                    dk(Some(selfcomp), Some(IndexResult::Pass()), 1.0, vec![], low),
                ],
                dconf,
            )
            .unwrap()
        };
        let report = run(&dconf, false);
        assert_eq!(report.kmer_statuses, vec!["PrimerDimer".to_string()]);
        assert_eq!(report.dimers, vec![(1, 1)]);
        // A low-frequency sequence is left out of the dimer check
        let report = run(&dconf, true);
        assert!(report.kmer_statuses.is_empty() && report.dimers.is_empty());
        // No dimer check when thermo checks are off
        dconf.thermo_check = false;
        let report = run(&dconf, false);
        assert!(report.kmer_statuses.is_empty() && report.dimers.is_empty());
    }

    #[test]
    fn test_default_hairpin_limit() {
        assert_eq!(DigestConfig::create_default().primer_hairpin_th_max, 55.0);
    }

    #[test]
    fn test_digest_r_to_count() {
        let dconf = DigestConfig::create_default();
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
        let dconf = DigestConfig::create_default();
        let seqs =
            vec!["CCAATGGTGCAAAAGGTATAATCANTAATGTCCAAGTGGTGCAAGAAGGTATAATCATTAATGT".as_bytes()];

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
            vec![std::str::from_utf8("ACCTTCTTGCACCACTTGGACATTA".as_bytes())
                .unwrap()
                .to_string()]
        );
    }

    #[test]
    fn test_digest_r_to_count_ambs() {
        // Tests ambs are handled in the initial kmer slice
        let dconf = DigestConfig::create_default();
        let seqs = vec!["ATTABAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTACGATCTCTTGTAGTCT".as_bytes()];
        //TTATAAGGTTTATACCTTCCCAGGTAACAAAC
        //TTACAAGGTTTATACCTTCCCAGGTAACA
        //TTAGAAGGTTTATACCTTCCCAGGTAACA
        let digested = digest_r_to_count(&seqs, 1, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 3);
        // Check count of ambiguous base
        for (_key, count) in digested.iter() {
            assert_eq!(count, &(1.0 / 3.0), "Count: {}", count);
        }
        for seq in digested.keys().into_iter() {
            println!("{:?}", str::from_utf8(seq.as_ref().unwrap()))
        }
        // Check Sequence
        let expected = vec![
            "GTTTGTTACCTGGGAAGGTATAAACCTTATAA".as_bytes().to_vec(),
            "TGTTACCTGGGAAGGTATAAACCTTGTAA".as_bytes().to_vec(),
            "TGTTACCTGGGAAGGTATAAACCTTCTAA".as_bytes().to_vec(),
        ];
        for e in expected.into_iter() {
            assert!(digested.contains_key(&Ok(e)),)
        }
    }
    #[test]
    fn test_digest_r_to_count_ambs_ext() {
        let dconf = DigestConfig::create_default();
        let seqs =
            vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTRCYGATCTCTTGTAGATCT".as_bytes()];
        // rc"AGATCTACAAGAGATCRGYAAAGTTGGTTGGTTTGTTACCTGGGAAGGTATAAACCTTTAAT"
        let digested = digest_r_to_count(&seqs, 30, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 4);
        // Check count of ambiguous base
        for (_key, count) in digested.iter() {
            assert_eq!(count, &0.25, "Count: {}", count);
        }

        for seq in digested.keys().into_iter() {
            println!("{:?}", str::from_utf8(seq.as_ref().unwrap()))
        }
        // Check Sequence
        let expected = vec![
            "ACAAGAGATCGGTAAAGTTGGTTGGT".as_bytes().to_vec(),
            "CAAGAGATCAGCAAAGTTGGTTGGT".as_bytes().to_vec(),
            "AGAGATCGGCAAAGTTGGTTGGT".as_bytes().to_vec(),
            "TCTACAAGAGATCAGTAAAGTTGGTTGGT".as_bytes().to_vec(),
        ];
        for e in expected.into_iter() {
            assert!(digested.contains_key(&Ok(e)),)
        }
    }

    #[test]
    fn test_digest_f_to_count() {
        let dconf = DigestConfig::create_default();
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
        let dconf = DigestConfig::create_default();
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
        let dconf = DigestConfig::create_default();
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAA-TTTCGATCTCTTGTAGATCT".as_bytes()];

        let digested = digest_f_to_count(&seqs, 40, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 1);
        // Check sequence
        assert_eq!(digested.contains_key(&Err(DigestError::GapOnSetBase)), true);
    }

    #[test]
    fn test_digest_f_to_count_ambs() {
        let dconf = DigestConfig::create_default();
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAABTTTCGATCTCTTGTAGATCT".as_bytes()];

        let digested = digest_f_to_count(&seqs, 40, &dconf);

        // Check num of seqs
        assert_eq!(digested.len(), 3);
        // Check count of ambiguous base
        for (_key, count) in digested.iter() {
            assert_eq!(count, &(1.0 / 3.0), "Count: {}", count);
        }
        println!("{:?}", digested);
        // Check sequences
        let expected = vec![
            "CTTCCCAGGTAACAAACCAACCAAT".as_bytes().to_vec(),
            "TTCCCAGGTAACAAACCAACCAAC".as_bytes().to_vec(),
            "CTTCCCAGGTAACAAACCAACCAAG".as_bytes().to_vec(),
        ];
        for e in expected.into_iter() {
            assert!(digested.contains_key(&Ok(e)),)
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
        let dconf = DigestConfig::create_default();
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCBARTTTCGATCTCTTGTAGATCT".as_bytes()];
        // Check that appending amb bases works as expected
        let digested = digest_f_to_count(&seqs, 60, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 6);

        // Check freqs
        for (_key, count) in digested.iter() {
            assert_eq!(*count, 1.0 / digested.len() as f64, "Count: {}", count);
        }
        for seq in digested.keys().into_iter() {
            println!("{:?}", str::from_utf8(seq.as_ref().unwrap()))
        }
        // Check Sequence
        let expected = vec![
            "CCAACCTAGTTTCGATCTCTTGTAGATCT".as_bytes().to_vec(),
            "AACCGAGTTTCGATCTCTTGTAGATCT".as_bytes().to_vec(),
            "ACCAACCTAATTTCGATCTCTTGTAGATCT".as_bytes().to_vec(),
            "CAACCGAATTTCGATCTCTTGTAGATCT".as_bytes().to_vec(),
            "AACCCAGTTTCGATCTCTTGTAGATCT".as_bytes().to_vec(),
            "CAACCCAATTTCGATCTCTTGTAGATCT".as_bytes().to_vec(),
        ];
        for e in expected.into_iter() {
            assert!(digested.contains_key(&Ok(e)),)
        }
    }

    #[test]
    fn test_digest_f_to_count_wl() {
        let dconf = DigestConfig::create_default();
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
        let dconf = DigestConfig::create_default();
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
        let dconf = DigestConfig::create_default();
        let seqs = vec!["ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAA!TTTCGATCTCTTGTAGATCT".as_bytes()];

        let digested = digest_f_to_count(&seqs, 40, &dconf);
        // Check num of seqs
        assert_eq!(digested.len(), 1);
        // Check sequence
        assert_eq!(digested.contains_key(&Err(DigestError::InvalidBase)), true);
    }

    #[test]
    fn test_opto_binding_tm_trunc() {
        let dconf = DigestConfig::create_default();

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
        let dconf = DigestConfig::create_default();

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
        let mut dconf = DigestConfig::create_default();
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
        let mut dconf = DigestConfig::create_default();
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
    #[test]
    fn test_opto_binding_dedeup() {
        // Test to show that elements made to be duplicates by truncation are correctly combined.
        let mut dconf = DigestConfig::create_default();
        dconf.annealing_temp_c = 65.0;
        dconf.primer_annealing_prop = Some(0.0); // Very low annealing forces truncation.
        dconf.thermo_type = ThermoType::ANNEALING;

        let seqs = vec![
            "GGGTTGTGATGGTGGCAGTTTGTA".as_bytes().to_vec(),
            "TGGTTGTGATGGTGGCAGTTTGTA".as_bytes().to_vec(),
        ];

        let mut count_map: HashMap<Result<Vec<u8>, DigestError>, f64> = HashMap::new();
        for seq in seqs.into_iter() {
            count_map.insert(Ok(seq), 1.0);
        }
        // Opto the sequences
        let processed_map = opto_binding(count_map.clone(), &dconf);

        for seq in processed_map.keys().into_iter() {
            println!("{:?}", str::from_utf8(seq.as_ref().unwrap()))
        }

        // Check wanted seqs are present
        let wanted_seqs = vec!["GGTTGTGATGGTGGCAGTTTGTA".as_bytes().to_vec()];
        for wseq in wanted_seqs.into_iter() {
            assert_eq!(processed_map.contains_key(&Ok(wseq.clone())), true);
            assert_eq!(processed_map.get(&Ok(wseq)), Some(&2.0));
        }
    }
    #[test]
    fn test_opto_binding_dedeup2() {
        // This example ends in seq1 being trucated but seq2 not. Therefore, seq1 ends up as a subseq of seq2.
        let mut dconf = DigestConfig::create_default();
        dconf.annealing_temp_c = 65.0;
        dconf.primer_annealing_prop = Some(30.0);
        dconf.thermo_type = ThermoType::ANNEALING;

        let seqs = vec![
            "GGGTTGTGATGGTGGCAGTTTGTA".as_bytes().to_vec(),
            "TGGTTGTGATGGTGGCAGTTTGTA".as_bytes().to_vec(),
        ];

        let wanted_seqs = vec![
            "GGTTGTGATGGTGGCAGTTTGTA".as_bytes().to_vec(),
            "TGGTTGTGATGGTGGCAGTTTGTA".as_bytes().to_vec(),
        ];

        let mut count_map: HashMap<Result<Vec<u8>, DigestError>, f64> = HashMap::new();
        for seq in seqs.into_iter() {
            count_map.insert(Ok(seq), 1.0);
        }
        // Opto the sequences
        let processed_map = opto_binding(count_map.clone(), &dconf);

        for seq in processed_map.keys().into_iter() {
            println!("{:?}", str::from_utf8(seq.as_ref().unwrap()))
        }
        // Check wanted seqs are present
        for wseq in wanted_seqs.into_iter() {
            assert_eq!(processed_map.contains_key(&Ok(wseq.clone())), true);

            assert_eq!(processed_map.get(&Ok(wseq)), Some(&1.0));
        }
    }

    #[test]
    fn test_dedup_substr() {
        let seqs = vec![
            "GGTTGTGATGGTGGCAGTTTGTA".as_bytes().to_vec(),
            "TGGTTGTGATGGTGGCAGTTTGTA".as_bytes().to_vec(),
        ];

        let mut count_map: HashMap<Result<Vec<u8>, DigestError>, f64> = HashMap::new();
        for seq in seqs.into_iter() {
            count_map.insert(Ok(seq), 1.0);
        }

        // Run dedup_substr
        count_map = dedup_substr(count_map);

        // Check outputs
        let wanted_seqs = vec!["GGTTGTGATGGTGGCAGTTTGTA".as_bytes().to_vec()];
        assert_eq!(count_map.len(), 1);
        for wseq in wanted_seqs.into_iter() {
            assert_eq!(count_map.contains_key(&Ok(wseq.clone())), true);
            assert_eq!(count_map.get(&Ok(wseq)), Some(&2.0));
        }
    }
}
