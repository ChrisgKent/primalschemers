pub fn create_mapping_array(seq: &[u8], remap: bool) -> Vec<Option<usize>> {
    let mut mapping_array: Vec<Option<usize>> = vec![None; seq.len()];

    // Mapping array of just the index at each index
    if !remap {
        for i in 0..seq.len() {
            mapping_array[i] = Some(i)
        }
        return mapping_array;
    }
    let mut i = 0;
    for (index, base) in seq.iter().enumerate() {
        match base {
            b' ' | b'-' => (),
            _ => {
                mapping_array[index] = Some(i);
                i += 1;
            }
        }
    }
    mapping_array
}

pub fn create_ref_to_msa(mapping_array: &Vec<Option<usize>>) -> Vec<usize> {
    let mut ref_to_msa: Vec<usize> = vec![0; mapping_array.len()];
    let mut final_pos = 0;

    for (index, ref_index) in mapping_array.iter().enumerate() {
        match ref_index {
            Some(ref_i) => {
                ref_to_msa[*ref_i] = index;
                final_pos = *ref_i;
            }
            None => {}
        }
    }
    // As array is created with size of mapping array, we need to resize it to the final position
    ref_to_msa.resize(final_pos + 1, final_pos);

    ref_to_msa
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_mapping_array() {
        let seq = b"ACGT";
        let mapping_array = create_mapping_array(seq, true);
        assert_eq!(mapping_array, vec![Some(0), Some(1), Some(2), Some(3)]);

        let seq = b"AC-GT";
        let mapping_array = create_mapping_array(seq, true);
        assert_eq!(
            mapping_array,
            vec![Some(0), Some(1), None, Some(2), Some(3)]
        );

        let seq = b"-AC-GT-";
        let mapping_array = create_mapping_array(seq, true);
        assert_eq!(
            mapping_array,
            vec![None, Some(0), Some(1), None, Some(2), Some(3), None]
        );

        let seq = b" -AC-GT---";
        let mapping_array = create_mapping_array(seq, true);
        assert_eq!(
            mapping_array,
            vec![
                None,
                None,
                Some(0),
                Some(1),
                None,
                Some(2),
                Some(3),
                None,
                None,
                None
            ]
        );
    }

    #[test]
    fn test_create_mapping_array_no_remap() {
        let seq = b"ACGT";
        let mapping_array = create_mapping_array(seq, false);
        assert_eq!(mapping_array, vec![Some(0), Some(1), Some(2), Some(3)]);

        let seq = b"AC-GT";
        let mapping_array = create_mapping_array(seq, false);
        assert_eq!(
            mapping_array,
            vec![Some(0), Some(1), Some(2), Some(3), Some(4)]
        );

        let seq = b"-AC-GT-";
        let mapping_array = create_mapping_array(seq, false);
        assert_eq!(
            mapping_array,
            vec![
                Some(0),
                Some(1),
                Some(2),
                Some(3),
                Some(4),
                Some(5),
                Some(6)
            ]
        );

        let seq = b" -AC-GT---";
        let mapping_array = create_mapping_array(seq, false);
        assert_eq!(
            mapping_array,
            vec![
                Some(0),
                Some(1),
                Some(2),
                Some(3),
                Some(4),
                Some(5),
                Some(6),
                Some(7),
                Some(8),
                Some(9)
            ]
        );
    }

    #[test]
    fn test_create_ref_to_msa() {
        let seq = b"ACGT";
        let mapping_array = create_mapping_array(seq, true);
        let ref_to_msa = create_ref_to_msa(&mapping_array);
        assert_eq!(ref_to_msa, vec![0, 1, 2, 3]);

        let seq = b"AC-GT";
        let mapping_array = create_mapping_array(seq, true);
        let ref_to_msa = create_ref_to_msa(&mapping_array);
        assert_eq!(ref_to_msa, vec![0, 1, 3, 4]);

        let seq = b"-AC-GT-A";
        let mapping_array = create_mapping_array(seq, true);
        let ref_to_msa = create_ref_to_msa(&mapping_array);
        assert_eq!(ref_to_msa, vec![1, 2, 4, 5, 7,]);
    }

    #[test]
    fn test_round_trip_mapping() {
        let gap_seq = b"  ATGC-TGTGCGAT-CGTAGCTA---GCTAGC--TGTAGC---TAGCTGATCA--";
        let seq = gap_seq
            .iter()
            .filter(|b| **b != b' ' && **b != b'-')
            .copied()
            .collect::<Vec<u8>>();

        let mapping_array = create_mapping_array(gap_seq, true);
        let ref_to_msa = create_ref_to_msa(&mapping_array);

        for (ref_i, &msa_i) in ref_to_msa.iter().enumerate() {
            assert_eq!(seq[ref_i], gap_seq[msa_i]);
        }
    }
}
