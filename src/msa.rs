use crate::config;
use crate::digest::IndexResult;
use crate::kmer::{FKmer, RKmer};
use crate::seqio;
use crate::{digest, mapping};

use pyo3::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::Path;

#[pyclass]
pub struct MSA {
    // Provided
    pub name: String,
    pub path: String,
    pub index: usize,
    _rust_log_file: Option<fs::File>,
    _log_verbosity: usize,

    // Calculated on init
    pub headers: Vec<String>,
    pub seqs: Vec<Vec<u8>>,
    pub _chromnames: String,
    pub _mapping_array: Vec<Option<usize>>,
    pub _ref_to_msa: Vec<usize>,
    dconf: config::DigestConfig,

    // Calculated on demand
    pub fkmers: Option<Vec<FKmer>>,
    pub rkmers: Option<Vec<RKmer>>,
}

#[pymethods]
impl MSA {
    #[new]
    #[pyo3(signature = (name, path, index, _log_verbosity, _rust_log_file=None))]
    pub fn new(
        name: String,
        path: String,
        index: usize,
        _log_verbosity: usize,
        _rust_log_file: Option<String>,
    ) -> MSA {
        let (headers, seqs) = seqio::fasta_reader(&path);
        let seq_array = seqs
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect::<Vec<Vec<u8>>>();

        let dconf = config::DigestConfig::new(
            None, None, None, None, None, None, None, None, None, None, None,
        );

        // Remove end insertions
        let seq_array = seqio::remove_end_insertions(seq_array);

        let mapping_array = mapping::create_mapping_array(&seqs[0].as_bytes());
        let ref_to_msa = mapping::create_ref_to_msa(&mapping_array);

        // Parse the headers
        let headers: Vec<String> = headers
            .into_iter()
            .map(|s| s.replace("/", "-").replace(">", ""))
            .collect();

        let _rust_log_file = match _rust_log_file {
            Some(log_file) => {
                let file: fs::File = fs::OpenOptions::new()
                    .write(true)
                    .create(true)
                    .truncate(true)
                    .open(&log_file)
                    .unwrap();
                Some(file)
            }
            None => None,
        };

        MSA {
            name,
            path,
            index,
            seqs: seq_array,
            _chromnames: headers[0].clone(),
            headers,
            _mapping_array: mapping_array,
            _ref_to_msa: ref_to_msa,
            fkmers: None,
            rkmers: None,
            _rust_log_file,
            _log_verbosity: _log_verbosity,
            dconf,
        }
    }

    pub fn digest_f(&mut self) -> () {
        let seq_slices = self
            .seqs
            .iter()
            .map(|s| s.as_slice())
            .collect::<Vec<&[u8]>>();

        let digested_f = digest::digest_f_primer(&seq_slices, &self.dconf, None);

        match self._log_verbosity {
            1 => {
                // Log sum stats
                let mut counter: HashMap<IndexResult, usize> = HashMap::new();
                for res in digested_f.iter() {
                    match res {
                        Ok(_) => {
                            let count = counter.entry(IndexResult::Pass()).or_insert(0);
                            *count += 1;
                        }
                        Err(e) => {
                            let count = counter.entry(e.clone()).or_insert(0);
                            *count += 1;
                        }
                    }
                }
                // Sort values and write to file
                let mut values = counter.into_iter().collect::<Vec<(IndexResult, usize)>>();
                values.sort_by(|a, b| b.1.cmp(&a.1));

                self._rust_log_file
                    .as_ref()
                    .unwrap()
                    .write_all(format!("fp:{:?}\n", values).as_bytes())
                    .unwrap();
            }
            _ => {} // Do nothing
        }

        // Assign the digested fkmers
        // self.fkmers = Some(digested_f.into_iter().filter_map(Result::ok).collect());
    }

    pub fn digest_r(&mut self) -> () {
        let seq_slices = self
            .seqs
            .iter()
            .map(|s| s.as_slice())
            .collect::<Vec<&[u8]>>();

        let digested_r = digest::digest_r_primer(&seq_slices, &self.dconf, None);

        match self._log_verbosity {
            1 => {
                // Log sum stats
                let mut counter: HashMap<IndexResult, usize> = HashMap::new();
                for (_i, res) in digested_r.iter().enumerate() {
                    match res {
                        Ok(_) => {
                            let count = counter.entry(IndexResult::Pass()).or_insert(0);
                            *count += 1;
                        }
                        Err(e) => {
                            let count = counter.entry(e.clone()).or_insert(0);
                            *count += 1;
                        }
                    }
                }
                // Sort values and write to file
                let mut values = counter.into_iter().collect::<Vec<(IndexResult, usize)>>();
                values.sort_by(|a, b| b.1.cmp(&a.1));

                self._rust_log_file
                    .as_ref()
                    .unwrap()
                    .write_all(format!("rp:{:?}\n", values).as_bytes())
                    .unwrap();
            }
            _ => {} // Do nothing
        }

        // Assign the digested rkmers
        // self.rkmers = Some(digested_r.into_iter().filter_map(Result::ok).collect());
    }

    pub fn digest(&mut self) -> () {
        // Digest the forward and reverse primers
        self.digest_f();
        self.digest_r();
    }

    pub fn write_msa_to_file(&self, path: &str) -> Result<(), std::io::Error> {
        let path = Path::new(path);

        let mut file = fs::OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(path)?;

        for (header, seq) in self.headers.iter().zip(&self.seqs) {
            writeln!(file, ">{}", header)?;

            for line in seq.chunks(60).collect::<Vec<&[u8]>>() {
                writeln!(
                    file,
                    "{}",
                    String::from_utf8(line.to_vec()).unwrap().replace(" ", "-")
                )?;
            }
        }
        Ok(())
    }
}
