use primer3_thal_rs::ThalArgs;

/// primalscheme3's default reaction conditions for the digest: 100 mM
/// monovalent, 2 mM divalent, 0.8 mM dNTP, 15 nM DNA and 0.8 M formamide
/// (as the digest has always used for Tm), at 37 °C with `max_loop` 30.
pub fn conditions() -> ThalArgs {
    ThalArgs {
        mv: 100.0,
        dv: 2.0,
        dntp: 0.8,
        dna_conc: 15.0,
        formamide_conc: 0.8,
        ..ThalArgs::default()
    }
}

pub enum ThermoType {
    TM,
    ANNEALING,
}

pub struct DigestConfig {
    pub primer_len_min: usize,
    pub primer_len_max: usize,
    pub primer_gc_max: f64,
    pub primer_gc_min: f64,
    // Tms
    pub primer_tm_max: f64,
    pub primer_tm_min: f64,
    // Annealing target
    pub primer_annealing_prop: Option<f64>,
    pub annealing_temp_c: f64,
    pub annealing_diff: f64,

    // Thermo mode
    pub thermo_type: ThermoType,
    pub max_homopolymers: usize,
    pub max_walk: usize,
    pub min_freq: f64,
    pub ignore_n: bool,
    pub dimerscore: f64,

    // Thermo check
    pub thermo_check: bool,
    // Max hairpin Tm (°C); kmers forming a more stable hairpin fail
    pub primer_hairpin_th_max: f64,
    // Return every kmer with per-sequence statuses instead of only passing kmers
    pub return_all: bool,
    /// Reaction conditions for primer Tm and annealing.
    pub thal: ThalArgs,
}

impl DigestConfig {
    pub fn new(
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
        annealing_diff: Option<f64>,
        // thermo_type
        thermo_type: Option<ThermoType>,

        max_walk: Option<usize>,
        max_homopolymers: Option<usize>,
        min_freq: Option<f64>,
        ignore_n: Option<bool>,
        dimerscore: Option<f64>,
        thermo_check: Option<bool>,
        primer_hairpin_th_max: Option<f64>,
        return_all: Option<bool>,
        thal: Option<ThalArgs>,
    ) -> DigestConfig {
        DigestConfig {
            primer_len_min: primer_len_min.unwrap_or(19),
            primer_len_max: primer_len_max.unwrap_or(34),
            primer_gc_max: primer_gc_max.unwrap_or(0.55),
            primer_gc_min: primer_gc_min.unwrap_or(0.35),
            // Tm
            primer_tm_max: primer_tm_max.unwrap_or(62.5),
            primer_tm_min: primer_tm_min.unwrap_or(59.5),
            // Annealing
            primer_annealing_prop: primer_annealing_prop,
            annealing_temp_c: annealing_temp_c.unwrap_or(65.0),
            annealing_diff: annealing_diff.unwrap_or(10.0),
            // Thermo
            thermo_type: thermo_type.unwrap_or(ThermoType::TM),

            max_homopolymers: max_homopolymers.unwrap_or(5),
            max_walk: max_walk.unwrap_or(80),
            min_freq: min_freq.unwrap_or(0.0),
            ignore_n: ignore_n.unwrap_or(false),
            dimerscore: dimerscore.unwrap_or(-26.0),
            thermo_check: thermo_check.unwrap_or(true),
            primer_hairpin_th_max: primer_hairpin_th_max.unwrap_or(55.0),
            return_all: return_all.unwrap_or(false),
            thal: thal.unwrap_or_else(conditions),
        }
    }
    /// Conditions for the hairpin check: [`DigestConfig::thal`] without
    /// DMSO or formamide, as primer3's `calc_hairpin` (and so
    /// primalscheme3's hairpin check), which has no such correction.
    pub fn hairpin_thal(&self) -> ThalArgs {
        ThalArgs {
            dmso_conc: 0.0,
            formamide_conc: 0.0,
            ..self.thal
        }
    }
    pub fn create_default() -> DigestConfig {
        DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None, None, None, None,
            None, None, None, None, None,
        )
    }
}
