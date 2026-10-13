use rayon::iter::{IntoParallelIterator, ParallelIterator};

pub fn fasta_reader(file: &str) -> (Vec<String>, Vec<String>) {
    try_fasta_reader(file).unwrap_or_else(|e| panic!("{e}"))
}

/// Headers and uppercased sequences of a FASTA file, or why it can't be
/// read: a missing file, or sequence before the first header.
pub fn try_fasta_reader(file: &str) -> Result<(Vec<String>, Vec<String>), String> {
    let mut seqs: Vec<Vec<String>> = Vec::new();
    let mut headers: Vec<String> = Vec::new();

    let contents =
        std::fs::read_to_string(file).map_err(|e| format!("Failed to read {file}: {e}"))?;
    for line in contents.lines() {
        if line.starts_with('>') {
            headers.push(line.to_string());
            seqs.push(Vec::new());
        } else if let Some(seq) = seqs.last_mut() {
            seq.push(line.to_string());
        } else if !line.trim().is_empty() {
            return Err(format!("{file}: sequence before the first FASTA header"));
        }
    }

    let seqs_final = seqs
        .iter_mut()
        .map(|seq_list| seq_list.concat().to_uppercase())
        .collect::<Vec<String>>();

    Ok((headers, seqs_final))
}

pub fn remove_end_insertions(mut seq_array: Vec<Vec<u8>>) -> Vec<Vec<u8>> {
    for seq in seq_array.iter_mut() {
        // Remove the right ends
        for base in seq.iter_mut() {
            match base {
                b'-' => {
                    *base = b' ';
                }
                _ => break,
            }
        }
        // Remove the left ends
        for base in seq.iter_mut().rev() {
            match base {
                b'-' => {
                    *base = b' ';
                }
                _ => break,
            }
        }
    }
    seq_array
}
pub fn calc_most_common_base(seq_array: &[Vec<u8>]) -> Vec<u8> {
    // If equal values. Will last in order of ACGT
    let indexes = 0..seq_array[0].len();
    let most_common_base = indexes
        .into_par_iter()
        .map(|basei| {
            let mut base_count: Vec<usize> = vec![0; 4];
            for seq in seq_array {
                match seq[basei] {
                    b'A' => base_count[0] += 1,
                    b'C' => base_count[1] += 1,
                    b'G' => base_count[2] += 1,
                    b'T' => base_count[3] += 1,
                    _ => (),
                }
            }
            let max_base = match base_count.iter().enumerate().max_by_key(|x| x.1).unwrap().0 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => panic!("Invalid base"),
            };
            max_base
        })
        .collect();

    most_common_base
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calc_most_common_base() {
        let seqs = vec![
            "AACATCGATCGATCGATCGATCGATCGATCGATCGATCGYTCGATCGATCGAT"
                .as_bytes()
                .to_vec(),
            "AAGCTCGATCGATCGATCGATCGATCGATCGATCGATCGYTCGATCGATCGAT"
                .as_bytes()
                .to_vec(),
            "  GGTCGATCGATCGATCGATCGATCGATCGATCGATCGYTCGATCGATCG  "
                .as_bytes()
                .to_vec(),
        ];

        let most_common_base = calc_most_common_base(&seqs);
        assert_eq!(most_common_base[..2], [b'A', b'A']);
        assert_eq!(most_common_base[2], b'G');
        assert_eq!(most_common_base[3], b'G');
    }
    #[test]
    fn test_remove_end_insertions() {
        let seqs = vec![
            "--CGATCGATCGATCGATCGATCGATCGATCGATCGATCGYTCGATCGATCGAT-GATCGATCGATCGATCGATCGATCGATCGATCGATCGAT--".as_bytes().to_vec(),
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATAGATCGAT-GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".as_bytes().to_vec(),
        ];

        let seqs_new = remove_end_insertions(seqs);
        assert_eq!(
            seqs_new,
            vec![
                "  CGATCGATCGATCGATCGATCGATCGATCGATCGATCGYTCGATCGATCGAT-GATCGATCGATCGATCGATCGATCGATCGATCGATCGAT  ".as_bytes().to_vec(),
                "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATAGATCGAT-GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".as_bytes().to_vec(),
            ]
        );
    }
}
