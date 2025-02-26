use digest::IndexResult;
use kmer::FKmer;

pub mod digest;
pub mod kmer;
pub mod seqfuncs;
pub mod tm;

fn fasta_reader(file: &str) -> (Vec<String>, Vec<String>) {
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
fn remove_end_insertions(mut seq_array: Vec<Vec<u8>>) -> Vec<Vec<u8>> {
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

fn calc_most_common_base(seq_array: &Vec<Vec<u8>>) -> Vec<u8> {
    let mut most_common_base: Vec<u8> = Vec::with_capacity(seq_array.len());
    for seqi in 0..seq_array.len() {
        let mut base_count: Vec<usize> = vec![0; 4];

        for basei in 0..seq_array[0].len() {
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

        most_common_base.push(max_base as u8);
    }
    most_common_base
}

fn main() {
    let (id, seqs) =
        fasta_reader("/Users/kentcg/schemes/artic-dengue/denv_all.filt.fasta.trimmed.together.aln");

    // let seqs = vec![
    //     "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGYTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string(),
    //     "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATAGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string(),
    // ];

    let mut seq_array: Vec<Vec<u8>> = seqs.iter().map(|seq| seq.as_bytes().to_vec()).collect();

    // Calculate the most common base
    let most_common_base = calc_most_common_base(&seq_array);

    seq_array = remove_end_insertions(seq_array);

    // Create thread pool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .unwrap();

    let seq_array_refs: Vec<&[u8]> = seq_array.iter().map(|seq| seq.as_slice()).collect();

    pool.install(|| {
        let digested_f = digest::digest_f_primer(&seq_array_refs);
        let digested_r = digest::digest_r_primer(&seq_array_refs);
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calc_most_common_base() {
        let seqs = vec![
            "AAGATCGATCGATCGATCGATCGATCGATCGATCGATCGYTCGATCGATCGAT"
                .as_bytes()
                .to_vec(),
            "AAGATCGATCGATCGATCGATCGATCGATCGATCGATCGYTCGATCGATCGAT"
                .as_bytes()
                .to_vec(),
            "  GATCGATCGATCGATCGATCGATCGATCGATCGATCGYTCGATCGATCG  "
                .as_bytes()
                .to_vec(),
        ];

        let most_common_base = calc_most_common_base(&seqs);
        assert_eq!(most_common_base[..2], [b'A', b'A']);
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
