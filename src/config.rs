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

    // Thermo mode
    pub thermo_type: ThermoType,

    pub max_homopolymers: usize,
    pub max_walk: usize,
    pub min_freq: f64,
    pub ignore_n: bool,
    pub dimerscore: f64,
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
        // thermo_type
        thermo_type: Option<ThermoType>,

        max_walk: Option<usize>,
        max_homopolymers: Option<usize>,
        min_freq: Option<f64>,
        ignore_n: Option<bool>,
        dimerscore: Option<f64>,
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
            // Thermo
            thermo_type: thermo_type.unwrap_or(ThermoType::TM),

            max_homopolymers: max_homopolymers.unwrap_or(5),
            max_walk: max_walk.unwrap_or(80),
            min_freq: min_freq.unwrap_or(0.0),
            ignore_n: ignore_n.unwrap_or(false),
            dimerscore: dimerscore.unwrap_or(-26.0),
        }
    }
}
