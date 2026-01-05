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

fn calc_dangling_ends_stability(
    seq1: &[u8],
    seq2: &[u8],
    start_x: usize,
    end_x: usize,
    offset: i32,
) -> f64 {
    let mut dg_score = 0.;

    // Look for overhang on the right side
    let last_x = end_x - 1;
    let last_s2_idx = (last_x as i32 + offset) as usize;

    match scores::seq2_overhang_dg(&seq1[last_x], &seq2[last_s2_idx], &seq2[last_s2_idx + 1]) {
        Some(score) => dg_score += score,
        None => dg_score += BONUS_ARRAY[2],
    }

    // Look for overhang on the leftside
    let first_x = start_x;
    let first_s2_idx = (first_x as i32 + offset) as usize;

    if first_x > 0 {
        match scores::seq1_overhang_dg(&seq1[first_x], &seq2[first_s2_idx], &seq1[first_x - 1]) {
            Some(score) => dg_score += score,
            None => dg_score += BONUS_ARRAY[1],
        }
    } else if first_s2_idx > 0 {
        match scores::seq2_overhang_dg(&seq1[first_x], &seq2[first_s2_idx], &seq2[first_s2_idx - 1])
        {
            Some(score) => dg_score += score,
            None => dg_score += BONUS_ARRAY[1],
        }
    }

    return dg_score;
}

fn calc_nn_thermo(seq1: &[u8], seq2: &[u8], start_x: usize, end_x: usize, offset: i32) -> f64 {
    let mut dg_score: f64 = 0.;
    for x in start_x..end_x {
        let s2_idx = (x as i32 + offset) as usize;
        match nn_dg_scores(&seq1[x..x + 2], &seq2[s2_idx..s2_idx + 2]) {
            Some(score) => dg_score += score,
            None => dg_score += BONUS_ARRAY[6],
        }
    }
    return dg_score;
}

fn calc_extension(
    seq1: &[u8],
    seq2: &[u8],
    start_x: usize,
    end_x: usize,
    offset: i32,
) -> Option<f64> {
    // Note: Early exit check is done in calc_at_offset.
    // We proceed to calculate score.

    let mut score: f64 = 0.;
    let len = end_x - start_x;
    let check_len = len.min(4);
    let mut all_match = true;

    for i in 0..check_len {
        let x = end_x - 1 - i;
        let s2_idx = (x as i32 + offset) as usize;
        let is_match = scores::match_array(seq1[x], seq2[s2_idx]);

        if !is_match {
            all_match = false;
        }

        if is_match {
            // Add match score
            match seq1[x] {
                b'G' | b'C' => score += 3. * (1. / (i + 1) as f64), // CG match
                b'A' | b'T' => score += 2. * (1. / (i + 1) as f64), // AT match
                _ => continue,
            }
        }
    }

    if all_match && (check_len == 4 || len < 4) {
        score += 2.;
    }

    return Some(-score);
}

fn apply_bonus(seq1: &[u8], seq2: &[u8], start_x: usize, end_x: usize, offset: i32) -> f64 {
    let mut current_match = 0;
    let mut longest_match = 0;
    let mut match_count = 0;
    let total_len = end_x - start_x;

    // For bubbles
    let mut current_run_val = false;
    let mut current_run_count = 0;
    let mut bubble_score = 0.;

    for x in start_x..end_x {
        let s2_idx = (x as i32 + offset) as usize;
        let is_match = scores::match_array(seq1[x], seq2[s2_idx]);

        // Longest match logic
        if is_match {
            current_match += 1;
            match_count += 1;
        } else {
            current_match = 0;
        }
        if current_match > longest_match {
            longest_match = current_match;
        }

        // Bubble logic
        if x == start_x {
            current_run_val = is_match;
            current_run_count = 1;
        } else {
            if is_match == current_run_val {
                current_run_count += 1;
            } else {
                // Process previous run
                if !current_run_val && current_run_count > 2 {
                    bubble_score +=
                        -((current_run_count as f64 - 2.) * BONUS_ARRAY[0]) * BONUS_ARRAY[9];
                }
                // Reset
                current_run_val = is_match;
                current_run_count = 1;
            }
        }
    }
    // Process last run for bubbles
    if total_len > 0 {
        if !current_run_val && current_run_count > 2 {
            bubble_score += -((current_run_count as f64 - 2.) * BONUS_ARRAY[0]) * BONUS_ARRAY[9];
        }
    }

    let mut score = 0.;

    // Find proportion of matches
    if total_len > 0 {
        score += -((0.8 - (match_count as f64 / total_len as f64)) * BONUS_ARRAY[8]);
    }

    // Longest match
    score += -(longest_match as f64 * BONUS_ARRAY[7]);

    // Resolve bubbles
    score += bubble_score;

    return score;
}

pub fn calc_at_offset(seq1: &[u8], seq2: &[u8], offset: i32) -> Option<f64> {
    // Calculate bounds
    let start_x_signed = if offset < 0 { -offset } else { 0 };
    let start_x = start_x_signed as usize;

    // seq2_index < seq2.len() - 1
    // x + offset < seq2.len() - 1
    // x < seq2.len() - 1 - offset
    let max_seq2_idx = seq2.len().saturating_sub(1);
    let end_x_limit = max_seq2_idx as i32 - offset;

    if end_x_limit <= 0 {
        return None;
    }

    let end_x = seq1.len().min(end_x_limit as usize);

    if start_x >= end_x {
        return None;
    }

    // Early exit check (3' end)
    // Last element is at end_x - 1
    let last_x = end_x - 1;
    let last_s2_idx = (last_x as i32 + offset) as usize;
    let last_match = scores::match_array(seq1[last_x], seq2[last_s2_idx]);

    let second_last_match = if last_x > start_x {
        let prev_x = last_x - 1;
        let prev_s2_idx = (prev_x as i32 + offset) as usize;
        scores::match_array(seq1[prev_x], seq2[prev_s2_idx])
    } else {
        false
    };

    if !last_match && !second_last_match {
        return None;
    }

    // Now calculate scores without allocations

    // 1. Dangling ends
    let mut dg_score = calc_dangling_ends_stability(seq1, seq2, start_x, end_x, offset);

    // 2. Extension
    match calc_extension(seq1, seq2, start_x, end_x, offset) {
        Some(score) => dg_score += score,
        None => return None,
    };

    // 3. Bonus
    dg_score += apply_bonus(seq1, seq2, start_x, end_x, offset);

    // 4. NN Thermo
    // Note: original code did mapping.pop() before NN thermo.
    // So we pass end_x - 1 to NN thermo.
    if end_x > start_x {
        dg_score += calc_nn_thermo(seq1, seq2, start_x, end_x - 1, offset);
    }

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
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_valid_encode_base() {
        let seq = "ATCG";

        // base_to_encode = {"A": 0, "T": 3, "C": 1, "G": 2}
        assert_eq!(encode_base(seq), vec![0, 3, 1, 2])
    }
    #[test]
    #[should_panic]
    fn test_invalid_encode_base() {
        encode_base("z");
    }
    #[test]
    fn test_all_match() {
        // Set up values
        let seq1 = "ACGAT";
        let seq2 = "TGCTA";
        let offset = 0;

        let mut pred_score: f64 = 0.;

        // AC / TG match
        pred_score += scores::nn_dg_scores(b"AC", b"TG").unwrap_or(0.);
        // CG / GC match
        pred_score += scores::nn_dg_scores(b"CG", b"GC").unwrap_or(0.);
        // GA / CT match
        pred_score += scores::nn_dg_scores(b"GA", b"CT").unwrap_or(0.);
        // AT / TA match
        pred_score += scores::nn_dg_scores(b"AT", b"TA").unwrap_or(0.);

        assert_eq!(
            calc_nn_thermo(seq1.as_bytes(), seq2.as_bytes(), 0, seq1.len() - 1, offset),
            pred_score
        )
    }
    #[test]
    fn test_mismatch_with_offset() {
        // Set up values
        //   ACCTC
        //   |||.|
        // ACTGGTGCTAC
        let seq1 = "ACCTC";
        let seq2 = "ACTGGTGCTAC";
        let offset = 2;

        let mut pred_score = 0.;

        // AC / TG match
        pred_score += scores::nn_dg_scores(b"AC", b"TG").unwrap_or(0.);
        // CC / GG  match
        pred_score += scores::nn_dg_scores(b"CC", b"GG").unwrap_or(0.);
        // CT / GT mismatch
        pred_score += scores::nn_dg_scores(b"CT", b"GT").unwrap_or(0.);
        // TC / TG mismatch
        pred_score += scores::nn_dg_scores(b"TC", b"TG").unwrap_or(0.);

        assert_eq!(
            calc_nn_thermo(seq1.as_bytes(), seq2.as_bytes(), 0, seq1.len() - 1, offset),
            pred_score
        )
    }
    #[test]
    fn test_match_array() {
        // MATCHES
        // A / T
        assert!(scores::match_array(b'A', b'T'));
        // T / A
        assert!(scores::match_array(b'T', b'A'));
        // C / G
        assert!(scores::match_array(b'C', b'G'));
        // G / C
        assert!(scores::match_array(b'G', b'C'));
        // MISMATCHES
        // A / A
        assert_eq!(scores::match_array(b'A', b'A'), false);
        // A / C
        assert_eq!(scores::match_array(b'A', b'C'), false);
        // A / G
        assert_eq!(scores::match_array(b'A', b'G'), false);

        // T / T
        assert_eq!(scores::match_array(b'T', b'T'), false);
        // T / C
        assert_eq!(scores::match_array(b'T', b'C'), false);
        // T / G
        assert_eq!(scores::match_array(b'T', b'G'), false);

        // C / C
        assert_eq!(scores::match_array(b'C', b'C'), false);
        // C / A
        assert_eq!(scores::match_array(b'C', b'A'), false);
        // C / T
        assert_eq!(scores::match_array(b'C', b'T'), false);

        // G / G
        assert_eq!(scores::match_array(b'G', b'G'), false);
        // G / A
        assert_eq!(scores::match_array(b'G', b'A'), false);
        // G / T
        assert_eq!(scores::match_array(b'G', b'T'), false);
    }
    #[test]
    fn test_ensure_consistant_result() {
        // nCoV-2019_76_RIGHT_0 nCoV-2019_18_LEFT_0
        // score: -40.74 (-40.736826004)
        // 5'-ACACCTGTGCCTGTTAAACCAT-3' >
        //                ||||||||||
        //             3'-CAATTTGGTAATTGAACACCCATAAAGGT-5'

        let s1 = "ACACCTGTGCCTGTTAAACCAT"; //5'-3'
        let s2 = "CAATTTGGTAATTGAACACCCATAAAGGT"; //3'-5'
        let offset = -12;

        assert_eq!(
            super::calc_at_offset(s1.as_bytes(), s2.as_bytes(), offset),
            Some(-40.736826004)
        );
    }
    #[test]
    fn test_ensure_detection() {
        // nCoV-2019_76_RIGHT_0 nCoV-2019_18_LEFT_0
        // score: -40.74 (-40.736826004)
        // 5'-ACACCTGTGCCTGTTAAACCAT-3' >
        //                ||||||||||
        //             3'-CAATTTGGTAATTGAACACCCATAAAGGT-5'

        let s1 = "ACACCTGTGCCTGTTAAACCAT"; //5'-3'
        let s2 = "TGGAAATACCCACAAGTTAATGGTTTAAC"; //5'-3'
        let threshold = -27.0;

        assert!(super::does_seq1_extend(
            s1.as_bytes(),
            s2.as_bytes(),
            threshold,
        ));
    }
    #[test]
    fn test_encode_decode() {
        // Test round trip encoding and decoding
        let seq = "CTCTTGTAGATCTGTTCTCTAAACGAACTTT";
        assert_eq!(decode_base(&encode_base(seq)), seq);
    }
}
