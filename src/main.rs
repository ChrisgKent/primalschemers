use rayon::iter::{IntoParallelIterator, ParallelIterator};

pub mod config;
pub mod digest;
pub mod kmer;
pub mod mapping;
pub mod primaldimer;
pub mod seqfuncs;
pub mod seqio;
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

fn main() {
    let (_id, seqs) = fasta_reader("/Users/kentcg/primerschemes/primerschemes/artic-measles/400/v1.0.0/work/all_genomes.align.ds.align.repaired.fasta");

    // let seqs = vec![
    //     "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGYTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string(),
    //     "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATAGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string(),
    // ];

    let mut seq_array: Vec<Vec<u8>> = seqs.iter().map(|seq| seq.as_bytes().to_vec()).collect();

    // Calculate the most common base

    seq_array = remove_end_insertions(seq_array);

    // Create thread pool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .unwrap();

    let seq_array_refs: Vec<&[u8]> = seq_array.iter().map(|seq| seq.as_slice()).collect();

    let dconf = config::DigestConfig::new(
        None, None, None, None, None, None, None, None, None, None, None, None, None, None,
    );

    pool.install(|| {
        let digested_f: Vec<Result<kmer::FKmer, digest::IndexResult>> =
            digest::digest_f_primer(&seq_array_refs, &dconf, None);
        let digested_r: Vec<Result<kmer::RKmer, digest::IndexResult>> =
            digest::digest_r_primer(&seq_array_refs, &dconf, None);

        let fkmers: Vec<kmer::FKmer> = digested_f.into_par_iter().filter_map(Result::ok).collect();
        let rkmers: Vec<kmer::RKmer> = digested_r.into_par_iter().filter_map(Result::ok).collect();

        // Generate all primerpairs
        let fkmers_refs: Vec<&kmer::FKmer> = fkmers.iter().collect();
        let rkmers_refs: Vec<&kmer::RKmer> = rkmers.iter().collect();

        println!(
            "Generated {} forward and {} reverse kmers",
            fkmers_refs.len(),
            rkmers_refs.len()
        );
    });
}
