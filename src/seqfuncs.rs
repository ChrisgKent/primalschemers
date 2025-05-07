pub fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        _ => base,
    }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.into_iter().rev().map(|&b| complement_base(b)).collect()
}

pub fn atcg_only(kmer: &[u8]) -> bool {
    // Will return true if the kmer contains only A, T, C, or G
    for base in kmer.iter() {
        match base {
            b'A' | b'C' | b'G' | b'T' => {}
            _ => return false,
        }
    }
    true
}

#[derive(Debug, PartialEq)]
pub enum KmerCheck {
    ATCGOnly,
    ContainsAmbiguities,
    ContainsInvalidBases,
    ContainsNs,
}

pub fn check_kmer(kmer: &[u8]) -> KmerCheck {
    // Any invalid bases will return ContainsInvalidBases
    // Any N will return ContainsNs
    // Any ambiguity base will return ContainsAmbiguities
    // If the kmer is ATCG only, return ATCGOnly

    let mut contains_ambs = false;
    let mut contains_ns = false;

    for base in kmer.iter() {
        // Check for non iupac and non gap bases
        if !is_iupac_base(base) {
            // If the base is a gap or space, continue
            match base {
                b'-' | b' ' => {
                    continue;
                }
                _ => return KmerCheck::ContainsInvalidBases,
            }
        }
        // Check for N
        if base == &b'N' {
            contains_ns = true;
            continue;
        }
        // Check for amb
        if is_amb(base) {
            contains_ambs = true;
            continue;
        }
    }
    // Return the appropriate KmerCheck
    match (contains_ambs, contains_ns) {
        (_, true) => KmerCheck::ContainsNs,
        (false, false) => KmerCheck::ATCGOnly,
        (true, false) => KmerCheck::ContainsAmbiguities,
    }
}

pub fn contains_ambs(kmer: &[u8]) -> bool {
    kmer.iter().any(|&base| is_amb(&base))
}
pub fn is_amb(base: &u8) -> bool {
    matches!(
        base,
        b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'D' | b'H' | b'V' | b'B'
    )
}

pub fn max_homopolymer(kmer: &[u8]) -> usize {
    let mut max = 0;
    let mut current = 0;
    let mut last_base = b'.';

    for base in kmer.iter() {
        if *base == last_base {
            current += 1;
        } else {
            current = 1;
            last_base = *base;
        }
        if current > max {
            max = current;
        }
    }
    max
}

pub fn is_iupac_base(base: &u8) -> bool {
    match base {
        b'A' | b'C' | b'G' | b'T' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'D' | b'H'
        | b'V' | b'B' | b'N' => true,
        _ => false,
    }
}

pub fn is_allowed_base(base: &u8) -> bool {
    // iupac, ' ', or '-'
    is_iupac_base(base) || *base == b'-' || *base == b' '
}

pub fn contains_allowed_bases(kmer: &[u8]) -> bool {
    kmer.iter().all(|&base| is_allowed_base(&base))
}

pub fn expand_amb_base(amb: u8) -> Option<Vec<u8>> {
    match amb {
        b'R' => Some(vec![b'A', b'G']),
        b'Y' => Some(vec![b'C', b'T']),
        b'S' => Some(vec![b'G', b'C']),
        b'W' => Some(vec![b'A', b'T']),
        b'K' => Some(vec![b'G', b'T']),
        b'M' => Some(vec![b'A', b'C']),
        b'D' => Some(vec![b'A', b'G', b'T']),
        b'H' => Some(vec![b'A', b'C', b'T']),
        b'V' => Some(vec![b'A', b'C', b'G']),
        b'B' => Some(vec![b'C', b'G', b'T']),
        // b'N' => Some(vec![b'A', b'C', b'G', b'T']),
        b'A' => Some(vec![b'A']),
        b'C' => Some(vec![b'C']),
        b'G' => Some(vec![b'G']),
        b'T' => Some(vec![b'T']),

        _ => None,
    }
}

pub fn expand_amb_sequence(seq: &[u8]) -> Option<Vec<Vec<u8>>> {
    let mut expanded: Vec<Vec<u8>> = Vec::new();
    for base in seq.iter() {
        let new_bases = match expand_amb_base(*base) {
            Some(bases) => bases,
            None => return None,
        };
        if expanded.len() == 0 {
            for b in new_bases.iter() {
                expanded.push(vec![*b]);
            }
        } else {
            let mut new_expanded: Vec<Vec<u8>> = Vec::new();
            for b in new_bases.iter() {
                for e in expanded.iter() {
                    let mut new_e = e.clone();
                    new_e.push(*b);
                    new_expanded.push(new_e);
                }
            }
            expanded = new_expanded;
        }
    }
    Some(expanded)
}

pub fn gc_content(kmer: &[u8]) -> f64 {
    // handle empty string
    if kmer.len() == 0 {
        return 0.0;
    }
    kmer.iter().filter(|&b| b == &b'G' || b == &b'C').count() as f64 / kmer.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expand_ambs_base() {
        assert_eq!(expand_amb_base(b'R'), Some(vec![b'A', b'G']));
        assert_eq!(expand_amb_base(b'Y'), Some(vec![b'C', b'T']));
        assert_eq!(expand_amb_base(b'S'), Some(vec![b'G', b'C']));
        assert_eq!(expand_amb_base(b'W'), Some(vec![b'A', b'T']));
        assert_eq!(expand_amb_base(b'K'), Some(vec![b'G', b'T']));
        assert_eq!(expand_amb_base(b'M'), Some(vec![b'A', b'C']));
        assert_eq!(expand_amb_base(b'D'), Some(vec![b'A', b'G', b'T']));
        assert_eq!(expand_amb_base(b'H'), Some(vec![b'A', b'C', b'T']));
        assert_eq!(expand_amb_base(b'V'), Some(vec![b'A', b'C', b'G']));
        assert_eq!(expand_amb_base(b'B'), Some(vec![b'C', b'G', b'T']));

        assert_eq!(expand_amb_base(b'A'), Some(vec![b'A']));
        assert_eq!(expand_amb_base(b'C'), Some(vec![b'C']));
        assert_eq!(expand_amb_base(b'G'), Some(vec![b'G']));
        assert_eq!(expand_amb_base(b'T'), Some(vec![b'T']));
        assert_eq!(expand_amb_base(b'X'), None);
    }

    #[test]
    fn test_expand_amb_sequence() {
        let seq = b"ATCG";
        let expanded = expand_amb_sequence(seq);
        assert_eq!(expanded.unwrap().len(), 1);

        let seq = b"ATCGR";
        let expanded = expand_amb_sequence(seq);
        assert_eq!(expanded.unwrap().len(), 2);

        let seq = b"RATCGY";
        let expanded = expand_amb_sequence(seq);
        assert_eq!(expanded.unwrap().len(), 4);

        let seq = b"ATCGN";
        let expanded = expand_amb_sequence(seq);
        assert_eq!(expanded, None);
    }

    #[test]
    fn test_gc_content() {
        let kmer = b"ATCG";
        assert_eq!(gc_content(kmer), 0.5);

        let kmer = b"ATCGGG";
        assert_eq!(gc_content(kmer), 0.6666666666666666);

        let kmer = b"GGGGGG";
        assert_eq!(gc_content(kmer), 1.0);

        let kmer = b"CCCC";
        assert_eq!(gc_content(kmer), 1.0);

        let kmer = b"AAAA";
        assert_eq!(gc_content(kmer), 0.0);

        let kmer = b"CGAGAACATTACCCATATGATAAGAGATTGT";
        assert_eq!(gc_content(kmer), 0.3548387096774194);

        let kmer = b"";
        assert_eq!(gc_content(kmer), 0.0);
    }

    #[test]
    fn test_max_homopolymer() {
        let kmer = b"ATCG";
        assert_eq!(max_homopolymer(kmer), 1);

        let kmer = b"";
        assert_eq!(max_homopolymer(kmer), 0);

        let kmer = b"ATCGGG";
        assert_eq!(max_homopolymer(kmer), 3);

        let kmer = b"GGGGGG";
        assert_eq!(max_homopolymer(kmer), 6);

        let kmer = b"CCCC";
        assert_eq!(max_homopolymer(kmer), 4);

        let kmer = b"AAAA";
        assert_eq!(max_homopolymer(kmer), 4);

        let kmer = b"CGAGAACATTACCCATATGATAAGAGATTGT";
        assert_eq!(max_homopolymer(kmer), 3);
    }

    #[test]
    fn test_contains_allowed_bases() {
        assert_eq!(contains_allowed_bases(b"ATCG"), true);

        assert_eq!(contains_allowed_bases(b"ATCGN"), true);

        assert_eq!(contains_allowed_bases(b"ATCGX"), false);

        assert_eq!(contains_allowed_bases(b"ATCG-"), true);

        assert_eq!(contains_allowed_bases(b"ATCGH"), true);

        assert_eq!(contains_allowed_bases(b"ATCG."), false);
    }

    #[test]
    fn test_check_kmer() {
        assert_eq!(check_kmer(b"ATCG"), KmerCheck::ATCGOnly);
        assert_eq!(check_kmer(b"ATCG "), KmerCheck::ATCGOnly);
        assert_eq!(check_kmer(b"ATCG-"), KmerCheck::ATCGOnly);
        assert_eq!(check_kmer(b"ATCG-N"), KmerCheck::ContainsNs);
        assert_eq!(check_kmer(b"ATCGX-"), KmerCheck::ContainsInvalidBases);
        assert_eq!(check_kmer(b"ATCGH"), KmerCheck::ContainsAmbiguities);
        assert_eq!(check_kmer(b"ATCG."), KmerCheck::ContainsInvalidBases);
        assert_eq!(check_kmer(b"ATCGB"), KmerCheck::ContainsAmbiguities);

        // Check Priority
        assert_eq!(check_kmer(b"ATCGNRY "), KmerCheck::ContainsNs);
        assert_eq!(check_kmer(b"-ATCGNR."), KmerCheck::ContainsInvalidBases);
    }
}
