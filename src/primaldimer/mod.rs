mod scores;
use itertools::Itertools;
use scores::nn_dg_scores;

static BONUS_ARRAY: [f64; 10] = [
    1.11217618,
    0.55187469,
    1.01582516,
    1.03180592,
    -2.76687727,
    -0.81903133,
    0.93596145,
    2.32758405,
    3.24507248,
    0.80416919,
];
// 0 = PENALTY_DOUBLE_MISMATCH
// 1 = PENALTY_LEFT_OVERHANG_MISMATCH
// 2 = PENALTY_RIGHT_OVERHANG_MISMATCH
// 3 = BONUS_ALL_MATCH
// 4 = BONUS_3P_MATCH_GC
// 5 = BONUS_3P_MATCH_AT
// 6 = SCORE_3P_MISMATCH
// 7 = LONGEST_MATCH_COEF
// 8 = MATCH_PROP_COEF
// 9 = BUBBLE_COEF

// 5' ATAATCAATCTAGACTATCGTATTTGCCTCC
// 5' CGTGATATTTTCTATCAATGGGGAAATTATTACG

// 5' ATAATCAATCTAGACTATCGTATTTGCCTCC
//                             ^i
//                          3' GCATTATTAAAGGGGTAACTATCTTTTATAGTGC
//                                  ^j

//base_to_u8 = {"A": 65, "T": 84, "C": 67, "G": 71}

pub fn calc_core_region(seq1: &[u8], seq2: &[u8], i: usize, j: usize) -> f64 {
    let mut dg_score = 0.;

    // Check for overhang on the right side
    if i < seq1.len() - 1 {
        match nn_dg_scores(&seq1[i..i + 2], &seq2[j..j + 2]) {
            Some(score) => dg_score += score,
            None => dg_score += BONUS_ARRAY[6],
        }
    } else {
        dg_score += BONUS_ARRAY[2];
    }

    // Check for overhang on the left side
    if j > 0 {
        match nn_dg_scores(&seq1[i..i + 2], &seq2[j..j + 2]) {
            Some(score) => dg_score += score,
            None => dg_score += BONUS_ARRAY[1],
        }
    } else {
        dg_score += BONUS_ARRAY[2];
    }

    return dg_score;
}

// base_to_encode = {"A": 0, "T": 3, "C": 1, "G": 2}
pub fn encode_base(sequence: &str) -> Vec<usize> {
    let encoded_base = sequence
        .as_bytes()
        .iter()
        .map(|base| match *base {
            65 => 0,
            84 => 3,
            67 => 1,
            71 => 2,
            _ => panic!("NON STANDRD BASE found in {}", sequence),
        })
        .collect();
    return encoded_base;
}

pub fn decode_base(encoded_base: &[usize]) -> String {
    let decoded_base = encoded_base
        .iter()
        .map(|base| match *base {
            0 => "A",
            3 => "T",
            1 => "C",
            2 => "G",
            _ => panic!("NON STANDRD BASE found in {:?}", encoded_base),
        })
        .collect::<Vec<&str>>()
        .join("");
    return decoded_base;
}

fn calc_dangling_ends_stability(seq1: &[u8], seq2: &[u8], mapping: &Vec<(usize, usize)>) -> f64 {
    let mut dg_score = 0.;

    // Look for overhang on the right side
    let (seq2_i, seq1_i) = mapping[mapping.len() - 1];

    match scores::seq2_overhang_dg(&seq1[seq1_i], &seq2[seq2_i], &seq2[seq2_i + 1]) {
        Some(score) => dg_score += score,
        None => dg_score += BONUS_ARRAY[2],
    }

    // Look for overhang on the leftside
    let (seq2_i, seq1_i) = mapping[0];

    if seq1_i > 0 {
        match scores::seq1_overhang_dg(&seq1[seq1_i], &seq2[seq2_i], &seq1[seq1_i - 1]) {
            Some(score) => dg_score += score,
            None => dg_score += BONUS_ARRAY[1],
        }
    } else if seq2_i > 0 {
        match scores::seq2_overhang_dg(&seq1[seq1_i], &seq2[seq2_i], &seq2[seq2_i - 1]) {
            Some(score) => dg_score += score,
            None => dg_score += BONUS_ARRAY[1],
        }
    }

    return dg_score;
}

fn calc_nn_thermo(seq1: &[u8], seq2: &[u8], mapping: &Vec<(usize, usize)>) -> f64 {
    let mut dg_score: f64 = 0.;
    for (seq2_i, seq1_i) in mapping.iter() {
        match nn_dg_scores(&seq1[*seq1_i..*seq1_i + 2], &seq2[*seq2_i..*seq2_i + 2]) {
            Some(score) => dg_score += score,
            None => dg_score += BONUS_ARRAY[6],
        }
    }
    return dg_score;
}

fn calc_extension(seq1: &[u8], match_bool: &Vec<bool>) -> Option<f64> {
    // Guard for no matches in final two 3' bases
    if !match_bool[match_bool.len() - 2..].iter().any(|f| *f) {
        return None;
    }

    let mut score: f64 = 0.;

    let kmer_3p_bool: Vec<(usize, &bool)> = match_bool
        .iter()
        .rev()
        .enumerate()
        .filter(|(index, _bool)| *index < 4)
        .collect();

    // Look at the last 4 bases in the match

    for (index, match_bool) in kmer_3p_bool.iter() {
        let seq1_index = seq1.len() - 1 - index;
        // Only count matches
        if **match_bool {
            // Add match score
            match seq1[seq1_index] {
                b'G' | b'C' => score += 3. * (1. / (index + 1) as f64), // CG match
                b'A' | b'T' => score += 2. * (1. / (index + 1) as f64), // AT match
                _ => continue,
            }
        }
    }

    if kmer_3p_bool.iter().all(|(_index, bool)| **bool) {
        score += 2.;
    }

    return Some(-score);
}

fn apply_bonus(match_bool: &Vec<bool>) -> f64 {
    // Find the longest continous match
    let mut current_match = 0;
    let mut longest_match = 0;

    for m in match_bool.iter() {
        match m {
            true => current_match += 1,
            false => current_match = 0,
        }
        if current_match > longest_match {
            longest_match = current_match
        }
    }
    let mut score = 0.;

    // Find proportion of matches
    score += -((0.8
        - (match_bool.iter().filter(|b| **b).count() as f64 / match_bool.len() as f64))
        * BONUS_ARRAY[8]);

    // Group the match bool
    let grouped_match_bool: Vec<(bool, usize)> = match_bool
        .iter()
        .chunk_by(|bool| **bool)
        .into_iter()
        .map(|(bool, iter)| (bool, iter.count()))
        .collect();

    // Work out the longest match
    let longest_match = grouped_match_bool
        .iter()
        .filter(|(bool, _count)| *bool)
        .map(|(_bool, count)| count)
        .max();

    match longest_match {
        Some(max) => score += -(*max as f64 * BONUS_ARRAY[7]),
        None => (),
    }

    // Resolve bubbles
    for (match_bool, count) in grouped_match_bool.iter() {
        if !*match_bool && count > &2 {
            score += -((*count as f64 - 2.) * BONUS_ARRAY[0]) * BONUS_ARRAY[9]
        }
    }

    return score;
}

pub fn calc_at_offset(seq1: &[u8], seq2: &[u8], offset: i32) -> Option<f64> {
    // Create the mapping
    let mut mapping: Vec<(usize, usize)> = Vec::with_capacity(seq1.len().max(seq2.len()));

    for x in 0..seq1.len() {
        let seq2_index = x as i32 + offset;
        if seq2_index >= 0 && (seq2_index as usize) < seq2.len() - 1 {
            mapping.push((seq2_index as usize, x))
        }
    }

    // Create the match_bool
    let match_bool = mapping
        .iter()
        .map(|(seq2i, seq1i)| scores::match_array(seq1[*seq1i], seq2[*seq2i]))
        .collect();

    let mut dg_score = calc_dangling_ends_stability(&seq1, &seq2, &mapping);

    match calc_extension(seq1, &match_bool) {
        Some(score) => dg_score += score,
        None => return None,
    };

    // Apply longest match, and match proportion
    dg_score += apply_bonus(&match_bool);

    // Remove the end element of mapping before giving to NN
    mapping.pop();
    dg_score += calc_nn_thermo(seq1, seq2, &mapping);

    return Some(dg_score);
}

pub fn does_seq1_extend(seq1: &[u8], seq2: &[u8], t: f64) -> bool {
    let mut seq2_rev = seq2.to_owned();
    seq2_rev.reverse();

    for offset in -(seq1.len() as i32 - 2)..(seq2.len() as i32) - (seq1.len() as i32) {
        match calc_at_offset(&seq1, &seq2_rev, offset) {
            Some(score) => {
                //println!("{}", score);
                if score <= t {
                    return true;
                }
            }
            None => (),
        }
    }
    return false;
}

pub fn does_seq1_extend_no_alloc(seq1: &[u8], seq2_rev: &[u8], t: f64) -> bool {
    for offset in -(seq1.len() as i32 - 2)..(seq2_rev.len() as i32) - (seq1.len() as i32) {
        match calc_at_offset(&seq1, &seq2_rev, offset) {
            Some(score) => {
                //println!("{}", score);
                if score <= t {
                    return true;
                }
            }
            None => (),
        }
    }
    return false;
}

pub fn do_seqs_interact(seq1: &str, seq2: &str, t: f64) -> bool {
    let s1 = seq1.as_bytes();
    let s2 = seq2.as_bytes();

    return does_seq1_extend(&s1, &s2, t) | does_seq1_extend(&s2, &s1, t);
}

pub fn do_pools_interact(pool1: Vec<&str>, pool2: Vec<&str>, t: f64) -> bool {
    // Encode the pools
    let pool1_encoded: Vec<Vec<u8>> = pool1.iter().map(|s| s.as_bytes().to_vec()).collect();
    let pool2_encoded: Vec<Vec<u8>> = pool2.iter().map(|s| s.as_bytes().to_vec()).collect();

    // Will look for interactions between every seq in pool1 and pool2
    for (s1, s2) in pool1_encoded.iter().cartesian_product(pool2_encoded.iter()) {
        if does_seq1_extend(&s1, &s2, t) | does_seq1_extend(&s2, &s1, t) {
            return true;
        }
    }
    return false;
}

pub fn do_pool_interact_u8(seqs1: &Vec<Vec<u8>>, seqs2: &Vec<Vec<u8>>, t: f64) -> bool {
    // Create the reverse complement of the sequences
    let mut seqs1_rev: Vec<Vec<u8>> = seqs1.iter().map(|s| s.to_vec()).collect();
    for s in seqs1_rev.iter_mut() {
        s.reverse();
    }
    let mut seqs2_rev: Vec<Vec<u8>> = seqs2.iter().map(|s| s.to_vec()).collect();
    for s in seqs2_rev.iter_mut() {
        s.reverse();
    }

    // Check for interactions
    for seq1i in 0..seqs1.len() {
        for seq2i in 0..seqs2.len() {
            if does_seq1_extend_no_alloc(&seqs1[seq1i], &seqs2_rev[seq2i], t)
                || does_seq1_extend_no_alloc(&seqs2[seq2i], &seqs1_rev[seq1i], t)
            {
                return true;
            }
        }
    }
    false
}

pub fn do_pool_interact_u8_slice(seqs1: &Vec<&[u8]>, seqs2: &Vec<&[u8]>, t: f64) -> bool {
    // Create the reverse complement of the sequences
    let mut seqs1_rev: Vec<Vec<u8>> = seqs1.iter().map(|s| s.to_vec()).collect();
    for s in seqs1_rev.iter_mut() {
        s.reverse();
    }
    let mut seqs2_rev: Vec<Vec<u8>> = seqs2.iter().map(|s| s.to_vec()).collect();
    for s in seqs2_rev.iter_mut() {
        s.reverse();
    }

    // Check for interactions
    for seq1i in 0..seqs1.len() {
        for seq2i in 0..seqs2.len() {
            if does_seq1_extend_no_alloc(&seqs1[seq1i], &seqs2_rev[seq2i], t)
                || does_seq1_extend_no_alloc(&seqs2[seq2i], &seqs1_rev[seq1i], t)
            {
                return true;
            }
        }
    }
    false
}
