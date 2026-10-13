//! Primer Tm and annealing on thal's nearest-neighbour model
//! ([`primer3_thal_rs::nn`]), so they agree with the hairpin and dimer
//! checks: an oligo's Tm equals thal's for the oligo paired with its
//! reverse complement.

use primer3_thal_rs::nn::{self, Duplex};
use primer3_thal_rs::{ThalArgs, ThalError};

/// An oligo grown one base at a time during the digest walk, keeping its
/// running nearest-neighbour sums so each extension is O(1).
#[derive(Clone, Debug)]
pub struct Oligo {
    pub seq: Vec<u8>,
    duplex: Duplex,
}

impl Oligo {
    /// Errors on any base other than A/C/G/T.
    pub fn new(seq: Vec<u8>) -> Result<Self, ThalError> {
        let duplex = Duplex::from_seq(&seq)?;
        Ok(Oligo { seq, duplex })
    }

    /// Append a base at the 3' end. Errors, leaving the oligo unchanged, on
    /// any base other than A/C/G/T.
    pub fn add_base(&mut self, base: u8) -> Result<(), ThalError> {
        self.duplex.push(base)?;
        self.seq.push(base);
        Ok(())
    }

    /// Tm (°C) at `args`, as [`oligo_tm`].
    pub fn tm(&self, args: &ThalArgs) -> f64 {
        self.duplex.tm(args, nn::is_self_complementary(&self.seq))
    }

    /// Percentage bound at `temp_c` °C, as [`oligo_annealing`].
    pub fn annealing(&self, args: &ThalArgs, temp_c: f64) -> f64 {
        self.duplex
            .annealing(args, nn::is_self_complementary(&self.seq), temp_c)
    }
}

/// Tm (°C) of `seq` with its complement at `args`, including any DMSO or
/// formamide correction. Errors on an empty sequence or a base other than
/// A/C/G/T.
pub fn oligo_tm(seq: &[u8], args: &ThalArgs) -> Result<f64, ThalError> {
    nn::duplex_tm(seq, args)
}

/// Percentage of `seq` bound to its complement at `temp_c` °C (50 at the Tm
/// without DMSO/formamide). Errors as [`oligo_tm`].
pub fn oligo_annealing(seq: &[u8], args: &ThalArgs, temp_c: f64) -> Result<f64, ThalError> {
    nn::duplex_annealing(seq, args, temp_c)
}

/// [`oligo_annealing`] at each of `temps` (°C), in order. Errors as
/// [`oligo_tm`].
pub fn oligo_annealing_profile(
    seq: &[u8],
    args: &ThalArgs,
    temps: &[f64],
) -> Result<Vec<f64>, ThalError> {
    let oligo = Oligo::new(seq.to_vec())?;
    if oligo.seq.is_empty() {
        return Err(ThalError::EmptySequence);
    }
    Ok(temps.iter().map(|&t| oligo.annealing(args, t)).collect())
}

/// Percentage of `seq1` bound to `seq2` (both 5'->3', mismatches allowed) at
/// `temp_c` °C: thal's heterodimer ΔH/ΔS through the same binding equation
/// as [`oligo_annealing`], which it equals when `seq2` is `seq1`'s reverse
/// complement. 0 if no structure forms. Errors as thal's heterodimer.
pub fn hetero_annealing(
    seq1: &[u8],
    seq2: &[u8],
    args: &ThalArgs,
    temp_c: f64,
) -> Result<f64, ThalError> {
    let res = primer3_thal_rs::calc_heterodimer(seq1, seq2, args, false)?;
    if !res.structure_found {
        return Ok(0.0);
    }
    let symmetric = nn::is_self_complementary(seq1) && nn::is_self_complementary(seq2);
    Ok(nn::bound_percent(res.dh, res.ds, args, symmetric, temp_c))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digest::conditions;
    use crate::seqfuncs::reverse_complement;

    #[test]
    fn test_oligo_tm_rev() {
        let kmer = b"TCATTGTATCCTCACATAACTCTCCCAAA".to_vec();
        let tm = oligo_tm(&kmer, &conditions()).unwrap();
        let tm_rev = oligo_tm(&reverse_complement(&kmer), &conditions()).unwrap();
        assert!((tm - tm_rev).abs() < 1e-9);
    }
}

/// Every Tm/annealing here must equal thal's for the same duplex.
#[cfg(test)]
mod consistency_tests {
    use super::*;
    use crate::digest::conditions;
    use crate::seqfuncs::reverse_complement;
    use primer3_thal_rs::{calc_heterodimer, ThalArgs};

    fn random_seqs(n: usize, mut state: u64) -> Vec<Vec<u8>> {
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as usize
        };
        (0..n)
            .map(|_| (0..15 + next() % 20).map(|_| b"ACGT"[next() % 4]).collect())
            .collect()
    }

    fn all_conditions() -> Vec<ThalArgs> {
        let d = ThalArgs::default();
        vec![
            conditions(),
            ThalArgs {
                formamide_conc: 0.0,
                ..conditions()
            },
            d,
            ThalArgs {
                dmso_conc: 5.0,
                formamide_conc: 1.5,
                ..d
            },
        ]
    }

    fn seqs() -> Vec<Vec<u8>> {
        let mut seqs = random_seqs(300, 17);
        for p in [
            "ACGTACGTACGTACGTACGT",
            "GAATTCGAATTCGAATTC",
            "AAAAAAAAAAAAAAAAAAAA",
        ] {
            seqs.push(p.as_bytes().to_vec());
        }
        seqs
    }

    #[test]
    fn test_oligo_tm_equals_thal_duplex_tm() {
        for args in all_conditions() {
            for seq in seqs() {
                let thal = calc_heterodimer(&seq, reverse_complement(&seq), &args, false)
                    .unwrap()
                    .tm;
                let tm = oligo_tm(&seq, &args).unwrap();
                let walk = Oligo::new(seq.clone()).unwrap().tm(&args);
                let label = String::from_utf8_lossy(&seq);
                assert!((tm - thal).abs() < 1e-9, "{label}: {tm} vs thal {thal}");
                assert_eq!(walk, tm, "{label}: walk Oligo differs");
            }
        }
    }

    #[test]
    fn test_oligo_add_base_matches_whole_sequence() {
        let args = conditions();
        for seq in seqs() {
            let mut oligo = Oligo::new(seq[..2].to_vec()).unwrap();
            for &b in &seq[2..] {
                oligo.add_base(b).unwrap();
            }
            assert_eq!(oligo.seq, seq);
            assert_eq!(oligo.tm(&args), oligo_tm(&seq, &args).unwrap());
            for t in [55.0, 65.0, 70.0] {
                assert_eq!(
                    oligo.annealing(&args, t),
                    oligo_annealing(&seq, &args, t).unwrap()
                );
            }
        }
    }

    #[test]
    fn test_annealing_is_half_bound_at_tm() {
        // Formamide shifts Tm, not binding, so compare without it
        let args = ThalArgs {
            formamide_conc: 0.0,
            ..conditions()
        };
        for seq in seqs() {
            let tm = oligo_tm(&seq, &args).unwrap();
            assert!((oligo_annealing(&seq, &args, tm).unwrap() - 50.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_long_oligos_are_not_truncated() {
        // The old Oligo silently stopped growing at 50 nt
        let seq = random_seqs(1, 9)[0].repeat(4);
        assert!(seq.len() > 50);
        let mut oligo = Oligo::new(seq[..1].to_vec()).unwrap();
        for &b in &seq[1..] {
            oligo.add_base(b).unwrap();
        }
        assert_eq!(oligo.seq.len(), seq.len());
        assert_eq!(
            oligo.tm(&conditions()),
            oligo_tm(&seq, &conditions()).unwrap()
        );
    }

    #[test]
    fn test_invalid_bases_are_errors() {
        assert!(oligo_tm(b"", &conditions()).is_err());
        assert!(oligo_tm(b"ACGTNACGT", &conditions()).is_err());
        assert!(Oligo::new(b"ACGR".to_vec()).is_err());
        let mut oligo = Oligo::new(b"ACGT".to_vec()).unwrap();
        assert!(oligo.add_base(b'-').is_err());
        assert_eq!(oligo.seq, b"ACGT");
    }

    #[test]
    fn test_annealing_profile_equals_annealing_at_each_temperature() {
        let temps: Vec<f64> = (40..=80).map(f64::from).collect();
        for args in all_conditions() {
            for seq in seqs() {
                let profile = oligo_annealing_profile(&seq, &args, &temps).unwrap();
                assert_eq!(profile.len(), temps.len());
                for (&t, &a) in temps.iter().zip(&profile) {
                    assert_eq!(a, oligo_annealing(&seq, &args, t).unwrap());
                }
            }
        }
        assert!(oligo_annealing_profile(b"", &conditions(), &temps).is_err());
        assert!(oligo_annealing_profile(b"ACGN", &conditions(), &temps).is_err());
        assert!(oligo_annealing_profile(b"ACGT", &conditions(), &[])
            .unwrap()
            .is_empty());
    }

    #[test]
    fn test_hetero_annealing_with_complement_equals_annealing() {
        for args in all_conditions() {
            for seq in seqs() {
                let rc = reverse_complement(&seq);
                for t in [50.0, 60.0, 70.0] {
                    let hetero = hetero_annealing(&seq, &rc, &args, t).unwrap();
                    let duplex = oligo_annealing(&seq, &args, t).unwrap();
                    assert!(
                        (hetero - duplex).abs() < 1e-9,
                        "{}: {hetero} vs {duplex}",
                        String::from_utf8_lossy(&seq)
                    );
                }
            }
        }
    }

    #[test]
    fn test_hetero_annealing_mismatches_and_errors() {
        let args = conditions();
        let seq = b"ACCAACGATGGTGTGTCCAT";
        let mut site = reverse_complement(seq);
        let perfect = hetero_annealing(seq, &site, &args, 55.0).unwrap();
        site[10] = if site[10] == b'A' { b'C' } else { b'A' }; // one mismatch
        let mismatched = hetero_annealing(seq, &site, &args, 55.0).unwrap();
        assert!(mismatched < perfect && mismatched > 0.0);
        // No structure binds nothing
        assert_eq!(
            hetero_annealing(b"AAAAAAAAAA", b"AAAAAAAAAA", &args, 55.0),
            Ok(0.0)
        );
        // thal's input errors pass through
        assert!(hetero_annealing(b"", b"ACGT", &args, 55.0).is_err());
    }
}
