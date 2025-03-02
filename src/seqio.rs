use rayon::iter::{IntoParallelIterator, ParallelIterator};

pub fn fasta_reader(file: &str) -> (Vec<String>, Vec<String>) {
    let mut seqs: Vec<Vec<String>> = Vec::new();
    let mut headers: Vec<String> = Vec::new();

    let file = std::fs::read_to_string(file).expect("Failed to read file");
    for line in file.lines() {
        if line.starts_with('>') {
            headers.push(line.to_string());
            seqs.push(Vec::new());
        } else {
            seqs.last_mut().unwrap().push(line.to_string());
        }
    }

    let seqs_final = seqs
        .iter_mut()
        .map(|seq_list| seq_list.concat().to_uppercase())
        .collect::<Vec<String>>();

    (headers, seqs_final)
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
pub fn calc_most_common_base(seq_array: &Vec<Vec<u8>>) -> Vec<u8> {
    // If equal values. Will last in order of ACGT
    let indexes = 0..seq_array[0].len();
    let most_common_base = indexes
        .into_par_iter()
        .map(|basei| {
            let mut base_count: Vec<usize> = vec![0; 4];
            for seqi in 0..seq_array.len() {
                match seq_array[seqi][basei] {
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
