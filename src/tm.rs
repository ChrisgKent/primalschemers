static SANTA_LUCIA_2004_DH: [[i32; 5]; 5] = [
    [76, 84, 78, 72, 72],  // AA, AC, AG, AT, AN
    [85, 80, 106, 78, 78], // CA, CC, CG, CT, CN
    [82, 98, 80, 84, 80],  // GA, GC, GG, GT, GN
    [72, 82, 85, 76, 72],  // TA, TC, TG, TT, TN
    [72, 80, 78, 72, 72],  // NA, NC, NG, NT, NN
];
// dS *-0.1 cal/k*mol
static SANTA_LUCIA_2004_DS: [[i32; 5]; 5] = [
    [213, 224, 210, 204, 224], // AA, AC, AG, AT, AN
    [227, 199, 272, 210, 272], // CA, CC, CG, CT, CN
    [222, 244, 199, 224, 244], // GA, GC, GG, GT, GN
    [213, 222, 227, 213, 227], // TA, TC, TG, TT, TN
    [168, 210, 220, 215, 220], // NA, NC, NG, NT, NN
];

static SANTA_LUCIA_1998_DH: [[i32; 5]; 5] = [
    [79, 84, 78, 72, 72],  /* AA, AC, AG, AT, AN; */
    [85, 80, 106, 78, 78], /* CA, CC, CG, CT, CN; */
    [82, 98, 80, 84, 80],  /* GA, GC, GG, GT, GN; */
    [72, 82, 85, 79, 72],  /* TA, TC, TG, TT, TN; */
    [72, 80, 78, 72, 72],
]; /* NA, NC, NG, NT, NN; */

/* dS *-0.1 cal/k*mol */
static SANTA_LUCIA_1998_DS: [[i32; 5]; 5] = [
    [222, 224, 210, 204, 224], /* AA, AC, AG, AT, AN; */
    [227, 199, 272, 210, 272], /* CA, CC, CG, CT, CN; */
    [222, 244, 199, 224, 244], /* GA, GC, GG, GT, GN; */
    [213, 222, 227, 222, 227], /* TA, TC, TG, TT, TN; */
    [168, 210, 220, 215, 220],
]; /* NA, NC, NG, NT, NN; */

// dG *-0.001 cal/mol
// static SANTA_LUCIA_1998_DG: [[i32; 5]; 5] = [
//     [1000, 1440, 1280, 880, 880],   // AA, AC, AG, AT, AN
//     [1450, 1840, 2170, 1280, 1450], // CA, CC, CG, CT, CN
//     [1300, 2240, 1840, 1440, 1300], // GA, GC, GG, GT, GN
//     [580, 1300, 1450, 1000, 580],   // TA, TC, TG, TT, TN
//     [580, 1300, 1280, 880, 580],    // NA, NC, NG, NT, NN
// ];

#[allow(dead_code)]
fn santa_lucia_1998_ds(nn: &[u8; 2]) -> i32 {
    match nn[0] {
        b'A' => match nn[1] {
            b'A' => 222,
            b'C' => 224,
            b'G' => 210,
            b'T' => 204,
            _ => 224,
        },
        b'C' => match nn[1] {
            b'A' => 227,
            b'C' => 199,
            b'G' => 272,
            b'T' => 210,
            _ => 272,
        },
        b'G' => match nn[1] {
            b'A' => 222,
            b'C' => 244,
            b'G' => 199,
            b'T' => 224,
            _ => 244,
        },
        b'T' => match nn[1] {
            b'A' => 213,
            b'C' => 222,
            b'G' => 227,
            b'T' => 222,
            _ => 227,
        },
        _ => 0,
    }
}

#[allow(dead_code)]
fn santa_lucia_1998_dh(nn: &[u8; 2]) -> i32 {
    match nn[0] {
        b'A' => match nn[1] {
            b'A' => 79,
            b'C' => 84,
            b'G' => 78,
            b'T' => 72,
            _ => 0,
        },
        b'C' => match nn[1] {
            b'A' => 85,
            b'C' => 80,
            b'G' => 106,
            b'T' => 78,
            _ => 0,
        },
        b'G' => match nn[1] {
            b'A' => 82,
            b'C' => 98,
            b'G' => 80,
            b'T' => 84,
            _ => 0,
        },
        b'T' => match nn[1] {
            b'A' => 72,
            b'C' => 82,
            b'G' => 85,
            b'T' => 79,
            _ => 0,
        },
        _ => 0,
    }
}

fn santa_lucia_2004_ds(nn: &[u8; 2]) -> i32 {
    match nn[0] {
        b'A' => match nn[1] {
            b'A' => 213,
            b'C' => 224,
            b'G' => 210,
            b'T' => 204,
            _ => 0,
        },
        b'C' => match nn[1] {
            b'A' => 227,
            b'C' => 199,
            b'G' => 272,
            b'T' => 210,
            _ => 0,
        },
        b'G' => match nn[1] {
            b'A' => 222,
            b'C' => 244,
            b'G' => 199,
            b'T' => 224,
            _ => 0,
        },
        b'T' => match nn[1] {
            b'A' => 213,
            b'C' => 222,
            b'G' => 227,
            b'T' => 213,
            _ => 0,
        },
        _ => 0,
    }
}

fn santa_lucia_2004_dh(nn: &[u8; 2]) -> i32 {
    match nn[0] {
        b'A' => match nn[1] {
            b'A' => 76,
            b'C' => 84,
            b'G' => 78,
            b'T' => 72,
            _ => 0,
        },
        b'C' => match nn[1] {
            b'A' => 85,
            b'C' => 80,
            b'G' => 106,
            b'T' => 78,
            _ => 0,
        },
        b'G' => match nn[1] {
            b'A' => 82,
            b'C' => 98,
            b'G' => 80,
            b'T' => 84,
            _ => 0,
        },
        b'T' => match nn[1] {
            b'A' => 72,
            b'C' => 82,
            b'G' => 85,
            b'T' => 76,
            _ => 0,
        },
        _ => 0,
    }
}

static T_KELVIN: f64 = 273.15;

#[derive(Clone, Debug)]
pub enum TmMethod {
    SantaLucia1998,
    SantaLucia2004,
}

use std::fmt;

impl fmt::Display for TmMethod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TmMethod::SantaLucia2004 => write!(f, "SantaLucia2004"),
            TmMethod::SantaLucia1998 => write!(f, "SantaLucia1998"),
            // other variants
        }
    }
}

static MAX_OLIGO_LEN: usize = 50;

#[derive(Clone, Debug)]
pub struct Oligo {
    pub seq: Vec<u8>,
    dh_array: Vec<i32>,
    ds_array: Vec<i32>,
}

impl Oligo {
    pub fn new(seq: Vec<u8>) -> Self {
        let mut dh_array = Vec::with_capacity(MAX_OLIGO_LEN);
        let mut ds_array = Vec::with_capacity(MAX_OLIGO_LEN);

        // Populate dh and ds with nn values
        for i in 0..seq.len() - 1 {
            dh_array.push(santa_lucia_2004_dh(&[seq[i], seq[i + 1]]));
            ds_array.push(santa_lucia_2004_ds(&[seq[i], seq[i + 1]]));
        }

        Oligo {
            seq,
            dh_array,
            ds_array,
        }
    }

    pub fn add_base(&mut self, base: u8) {
        if self.seq.len() >= MAX_OLIGO_LEN {
            return;
        }
        self.seq.push(base);

        if self.seq.len() > 1 {
            self.dh_array
                .push(santa_lucia_2004_dh(&[self.seq[self.seq.len() - 2], base]));
            self.ds_array
                .push(santa_lucia_2004_ds(&[self.seq[self.seq.len() - 2], base]));
        }
    }

    pub fn calc_tm(
        &self,
        dna_nm: f64,
        k_mm: f64,
        divalent_conc: f64,
        dntp_conc: f64,
        dmso_conc: f64,
        dmso_fact: f64,
        formamide_conc: f64,
    ) -> f64 {
        //
        let di_to_mo = divalent_to_monovalent(divalent_conc, dntp_conc).unwrap();
        let sym = symmetry_utf8(&self.seq);
        let gc_count = self.seq.iter().filter(|&&x| x == b'C' || x == b'G').count();

        // SantaLucia2004
        let mut ds = 57;
        let mut dh = -2;

        if sym {
            ds += 14;
        }
        // Terminal penalty
        match self.seq[0] {
            b'A' | b'T' => {
                ds += -69;
                dh += -22;
            }
            _ => {}
        }

        // Sum the arrays
        ds += self.ds_array.iter().sum::<i32>();
        dh += self.dh_array.iter().sum::<i32>();

        // End Terminal Pen
        match self.seq[self.seq.len() - 1] {
            b'A' | b'T' => {
                ds += -69;
                dh += -22;
            }
            _ => {}
        }

        // Convert
        let delta_h = dh as f64 * -100.0;
        let mut delta_s = ds as f64 * -0.1;

        // Salt correction using santalucia
        let adj_k_mm = k_mm + di_to_mo;
        delta_s += 0.368 * ((self.seq.len() - 1) as f64) * (adj_k_mm / 1000.0).ln();

        let dna_sym_adj = match sym {
            true => 1000000000.0,
            false => 4000000000.0,
        };

        let mut tm = delta_h / (delta_s + 1.987 * (dna_nm / dna_sym_adj).ln()) - T_KELVIN;
        tm -= dmso_conc * dmso_fact;
        tm += (0.453 * ((gc_count / self.seq.len()) as f64) - 2.88) * formamide_conc;

        // Calc bound, and return if asked
        // match annealing_temp {
        //     Some(temp) => {
        //         let ddg = delta_h - (temp + T_KELVIN) * delta_s;
        //         let ka = ((-ddg) / (1.987 * (temp + T_KELVIN))).exp();
        //         let bound = (1.0 / (1.0 + (1.0 / ((dna_nm / dna_sym_adj) * ka)).sqrt())) * 100.0;
        //         return bound;
        //     }
        //     None => {}
        // }

        tm
    }

    pub fn calc_annealing(
        &self,
        annealing_temp: f64,
        dna_nm: f64,
        k_mm: f64,
        divalent_conc: f64,
        dntp_conc: f64,
    ) -> f64 {
        //
        let di_to_mo = divalent_to_monovalent(divalent_conc, dntp_conc).unwrap();
        let sym = symmetry_utf8(&self.seq);

        // SantaLucia2004
        let mut ds = 57;
        let mut dh = -2;

        if sym {
            ds += 14;
        }
        // Terminal penalty
        match self.seq[0] {
            b'A' | b'T' => {
                ds += -69;
                dh += -22;
            }
            _ => {}
        }

        // Sum the arrays
        ds += self.ds_array.iter().sum::<i32>();
        dh += self.dh_array.iter().sum::<i32>();

        // End Terminal Pen
        match self.seq[self.seq.len() - 1] {
            b'A' | b'T' => {
                ds += -69;
                dh += -22;
            }
            _ => {}
        }

        // Convert
        let delta_h = dh as f64 * -100.0;
        let mut delta_s = ds as f64 * -0.1;

        // Salt correction using santalucia
        let adj_k_mm = k_mm + di_to_mo;
        delta_s += 0.368 * ((self.seq.len() - 1) as f64) * (adj_k_mm / 1000.0).ln();

        let dna_sym_adj = match sym {
            true => 1000000000.0,
            false => 4000000000.0,
        };

        let ddg = delta_h - (annealing_temp + T_KELVIN) * delta_s;
        let ka = ((-ddg) / (1.987 * (annealing_temp + T_KELVIN))).exp();
        let bound = (1.0 / (1.0 + (1.0 / ((dna_nm / dna_sym_adj) * ka)).sqrt())) * 100.0;

        bound
    }
}

fn divalent_to_monovalent(divalent: f64, dntp: f64) -> Result<f64, &'static str> {
    if divalent == 0.0 {
        return Ok(0.0);
    }
    if divalent < 0.0 || dntp < 0.0 {
        return Err("OLIGOTM_ERROR");
    }
    let adjusted_divalent = if divalent < dntp { dntp } else { divalent };
    Ok(120.0 * (adjusted_divalent - dntp).sqrt())
}

fn symmetry_utf8(seq: &[u8]) -> bool {
    let seq_len = seq.len();
    if seq_len % 2 == 1 {
        return false;
    }
    let mp = seq_len / 2;
    for i in 0..mp {
        let s = seq[i];
        let e = seq[seq_len - 1 - i];
        if (s == b'A' && e != b'T')
            || (s == b'T' && e != b'A')
            || (e == b'A' && s != b'T')
            || (e == b'T' && s != b'A')
        {
            return false;
        }
        if (s == b'C' && e != b'G')
            || (s == b'G' && e != b'C')
            || (e == b'C' && s != b'G')
            || (e == b'G' && s != b'C')
        {
            return false;
        }
    }
    true
}

fn symmetry(seq: &Vec<usize>) -> bool {
    let seq_len = seq.len();
    if seq_len % 2 == 1 {
        return false;
    }
    let mp = seq_len / 2;
    for i in 0..mp {
        let s = seq[i];
        let e = seq[seq_len - 1 - i];
        if (s == 0 && e != 3) || (s == 3 && e != 0) || (e == 0 && s != 3) || (e == 3 && s != 0) {
            return false;
        }
        if (s == 1 && e != 2) || (s == 2 && e != 1) || (e == 1 && s != 2) || (e == 2 && s != 1) {
            return false;
        }
    }
    true
}

pub fn oligo_thermo(
    seq_array: &[u8],
    dna_nm: f64,
    k_mm: f64,
    divalent_conc: f64,
    dntp_conc: f64,
    dmso_conc: f64,
    dmso_fact: f64,
    formamide_conc: f64,
    annealing_temp: f64,
    tm_method: TmMethod,
) -> (f64, f64) {
    // Todo Salt correction

    let di_to_mo = divalent_to_monovalent(divalent_conc, dntp_conc).unwrap();

    let mut gc_count = 0;
    // Convert to intarray:
    let oligo_int: Vec<usize> = seq_array
        .into_iter()
        .map(|b| match b {
            b'A' => 0,
            b'C' => {
                gc_count += 1;
                1
            }
            b'G' => {
                gc_count += 1;
                2
            }
            b'T' => 3,
            _ => 10,
        })
        .collect();

    let sym = symmetry(&oligo_int);

    let mut ds = 0;
    let mut dh = 0;

    match tm_method {
        TmMethod::SantaLucia1998 => {
            if sym {
                ds += 14
            }
            // Terminal penalty
            match oligo_int[0] {
                0 | 3 => {
                    ds += -41;
                    dh += -23;
                }
                1 | 2 => {
                    ds += 28;
                    dh += -1;
                }
                _ => {}
            }
            // Sum the pairs up
            for i in 0..oligo_int.len() - 1 {
                ds += SANTA_LUCIA_1998_DS[oligo_int[i]][oligo_int[i + 1]];
                dh += SANTA_LUCIA_1998_DH[oligo_int[i]][oligo_int[i + 1]];
            }

            // End Terminal Pen
            match oligo_int[oligo_int.len() - 1] {
                0 | 3 => {
                    ds += -41;
                    dh += -23;
                }
                1 | 2 => {
                    ds += 28;
                    dh += -1;
                }
                _ => {}
            }
        }
        TmMethod::SantaLucia2004 => {
            ds += 57;
            dh += -2;

            if sym {
                ds += 14;
            }
            match oligo_int[0] {
                0 | 3 => {
                    ds += -69;
                    dh += -22;
                }
                _ => {}
            }
            // Sum the pairs up
            for i in 0..oligo_int.len() - 1 {
                ds += SANTA_LUCIA_2004_DS[oligo_int[i]][oligo_int[i + 1]];
                dh += SANTA_LUCIA_2004_DH[oligo_int[i]][oligo_int[i + 1]];
            }
            match oligo_int[oligo_int.len() - 1] {
                0 | 3 => {
                    ds += -69;
                    dh += -22;
                }
                _ => {}
            }
        }
    }

    // Convert to
    let delta_h = dh as f64 * -100.0;
    let mut delta_s = ds as f64 * -0.1;

    // Return values
    let mut bound = 0.0;

    // Salt correction using santalucia
    let adj_k_mm = k_mm + di_to_mo;
    delta_s += 0.368 * ((oligo_int.len() - 1) as f64) * (adj_k_mm / 1000.0).ln();

    let dna_sym_adj = match sym {
        true => 1000000000.0,
        false => 4000000000.0,
    };

    let mut tm = delta_h / (delta_s + 1.987 * (dna_nm / dna_sym_adj).ln()) - T_KELVIN;
    tm -= dmso_conc * dmso_fact;
    tm += (0.453 * ((gc_count / oligo_int.len()) as f64) - 2.88) * formamide_conc;

    // Calc bound
    if annealing_temp > 0.0 {
        let ddg = delta_h - (annealing_temp + T_KELVIN) * delta_s;
        let ka = ((-ddg) / (1.987 * (annealing_temp + T_KELVIN))).exp();
        bound = (1.0 / (1.0 + (1.0 / ((dna_nm / dna_sym_adj) * ka)).sqrt())) * 100.0;
    }

    //

    (tm, bound)
}

fn _encode_base(b: u8) -> usize {
    // Converts ascii base to index
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

pub fn oligo_tm_utf8(
    sequence: &[u8],
    dna_nm: f64,
    k_mm: f64,
    divalent_conc: f64,
    dntp_conc: f64,
    dmso_conc: f64,
    dmso_fact: f64,
    formamide_conc: f64,
    tm_method: TmMethod,
) -> f64 {
    oligo_thermo(
        sequence,
        dna_nm,
        k_mm,
        divalent_conc,
        dntp_conc,
        dmso_conc,
        dmso_fact,
        formamide_conc,
        -10.0,
        tm_method,
    )
    .0
}

pub fn oligo_annealing_utf8(
    sequence: &[u8],
    dna_nm: f64,
    k_mm: f64,
    divalent_conc: f64,
    dntp_conc: f64,
    dmso_conc: f64,
    dmso_fact: f64,
    formamide_conc: f64,
    annealing_temp: f64,
    tm_method: TmMethod,
) -> f64 {
    oligo_thermo(
        sequence,
        dna_nm,
        k_mm,
        divalent_conc,
        dntp_conc,
        dmso_conc,
        dmso_fact,
        formamide_conc,
        annealing_temp,
        tm_method,
    )
    .1
}

#[cfg(test)]
mod tests {
    use crate::seqfuncs::reverse_complement;

    use super::*;

    #[test]
    fn test_sl_1998_ds() {
        for n1 in [b'A', b'C', b'G', b'T'].iter() {
            for n2 in [b'A', b'C', b'G', b'T'].iter() {
                let nn = [*n1 as u8, *n2 as u8];
                let r = santa_lucia_1998_ds(&nn);
                let q = SANTA_LUCIA_1998_DS[_encode_base(nn[0])][_encode_base(nn[1])];
                assert_eq!(r, q);
            }
        }

        let nn = [b'A', b'A'];
        let r = santa_lucia_1998_ds(&nn);
        let q = SANTA_LUCIA_1998_DS[0][0];
        assert_eq!(r, q);
    }
    #[test]
    fn test_sl_1998_dh() {
        for n1 in [b'A', b'C', b'G', b'T'].iter() {
            for n2 in [b'A', b'C', b'G', b'T'].iter() {
                let nn = [*n1 as u8, *n2 as u8];
                let r = santa_lucia_1998_dh(&nn);
                let q = SANTA_LUCIA_1998_DH[_encode_base(nn[0])][_encode_base(nn[1])];
                assert_eq!(r, q);
            }
        }

        let nn = [b'A', b'A'];
        let r = santa_lucia_1998_dh(&nn);
        let q = SANTA_LUCIA_1998_DH[0][0];
        assert_eq!(r, q);
    }

    #[test]
    fn test_sl_2004_ds() {
        for n1 in [b'A', b'C', b'G', b'T'].iter() {
            for n2 in [b'A', b'C', b'G', b'T'].iter() {
                let nn = [*n1 as u8, *n2 as u8];
                let r = santa_lucia_2004_ds(&nn);
                let q = SANTA_LUCIA_2004_DS[_encode_base(nn[0])][_encode_base(nn[1])];
                assert_eq!(r, q);
            }
        }

        let nn = [b'A', b'A'];
        let r = santa_lucia_2004_ds(&nn);
        let q = SANTA_LUCIA_2004_DS[0][0];
        assert_eq!(r, q);
    }
    #[test]
    fn test_sl_2004_dh() {
        for n1 in [b'A', b'C', b'G', b'T'].iter() {
            for n2 in [b'A', b'C', b'G', b'T'].iter() {
                let nn = [*n1 as u8, *n2 as u8];
                let r = santa_lucia_2004_dh(&nn);
                let q = SANTA_LUCIA_2004_DH[_encode_base(nn[0])][_encode_base(nn[1])];
                assert_eq!(r, q);
            }
        }

        let nn = [b'A', b'A'];
        let r = santa_lucia_2004_dh(&nn);
        let q = SANTA_LUCIA_2004_DH[0][0];
        assert_eq!(r, q);
    }
    #[test]
    fn test_oligo_vs_oligotm() {
        let seq = "TCATTGTATCCTCACATAACTCTCCCAAA".as_bytes();

        let oligo = Oligo::new(seq.to_vec());
        let tm = oligo_tm_utf8(
            &seq,
            15.0,
            100.0,
            2.0,
            0.8,
            0.0,
            0.0,
            0.8,
            TmMethod::SantaLucia2004,
        );

        let tm2 = oligo.calc_tm(15.0, 100.0, 2.0, 0.8, 0.0, 0.0, 0.8);
        assert_eq!(tm, tm2);
    }
    #[test]
    fn test_oligo_tm_rev() {
        let kmer = "TCATTGTATCCTCACATAACTCTCCCAAA".as_bytes().to_vec();
        let tm = oligo_tm_utf8(
            &kmer,
            15.0,
            100.0,
            2.0,
            0.8,
            0.0,
            0.0,
            0.8,
            TmMethod::SantaLucia2004,
        );

        let rc_kmer = reverse_complement(&kmer);
        let tm_rev = oligo_tm_utf8(
            &rc_kmer,
            15.0,
            100.0,
            2.0,
            0.8,
            0.0,
            0.0,
            0.8,
            TmMethod::SantaLucia2004,
        );

        assert_eq!(tm, tm_rev);
    }
}
