use pyo3::prelude::*;
use std::collections::HashMap;

use crate::{kmer, mapping, seqio};

#[pyclass]
struct MSA {
    name: String,
    path: String,
    msa_index: usize,

    // Calc on init
    _thread_pool: rayon::ThreadPool,
    _seq_array: Vec<Vec<u8>>,
    _uuid: String,
    _seq_id_to_index: HashMap<String, usize>,

    // Calc on eval
    _chrom_name: Option<String>,
    _mapping_array: Option<Vec<Option<usize>>>,
    _ref_to_msa_array: Option<Vec<usize>>,
    fkmers: Option<Vec<kmer::FKmer>>,
    rkmers: Option<Vec<kmer::RKmer>>,
}

impl MSA {
    fn new(name: &str, path: String, msa_index: usize, ncores: usize, uuid: String) -> Self {
        // Create the thread_pool
        let _thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(ncores)
            .build()
            .unwrap();

        // Read in the MSA
        let (_headers, seqs) = seqio::fasta_reader(&path);

        // Create seq_dict
        let _seq_id_to_index: HashMap<String, usize> = _headers
            .into_iter()
            .enumerate()
            .map(|(i, h)| (h, i))
            .collect();

        // parse the fasta to a seq array
        let mut _seq_array: Vec<Vec<u8>> = seqs
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect::<Vec<Vec<u8>>>();
        _seq_array = seqio::remove_end_insertions(_seq_array);

        MSA {
            name: name.to_string(),
            path: path,
            msa_index,
            _thread_pool,
            _seq_array,
            _uuid: uuid,
            _seq_id_to_index,
            _chrom_name: None,
            _mapping_array: None,
            _ref_to_msa_array: None,
            fkmers: None,
            rkmers: None,
        }
    }
}
