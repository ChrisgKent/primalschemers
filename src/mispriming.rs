//! Off-target (mispriming) screening of primers' 3' terminal kmers.
//!
//! The primers are indexed, not the references: each reference is streamed
//! once past a hash of the query kmers (and their single-mismatch variants),
//! so a host-sized genome needs no more memory than its largest contig.
use crate::mapping::create_mapping_array;
use crate::seqfuncs::{expand_amb_base, expand_amb_sequence, is_amb};
use crate::seqio::try_fasta_reader;
use flate2::read::MultiGzDecoder;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};

/// Bases per parallel work unit when scanning a background contig
const CHUNK: usize = 1 << 20;
/// Read buffer for background FASTAs
const READ_BUFFER: usize = 1 << 20;
/// Most sequences an ambiguous window or primer is expanded into
const MAX_EXPANSION: usize = 16;

/// A primer binding site: `pos` is the reference position of the primer's
/// 3' end base. `plus` primers extend rightwards, minus primers leftwards.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Hit {
    pub ref_id: u32,
    pub contig_id: u32,
    pub pos: u32,
    pub plus: bool,
}

/// A hit as Python sees it: (ref_id, contig_id, pos, "+" | "-")
type PyHit = (u32, u32, u32, char);

impl Hit {
    fn to_py(self) -> PyHit {
        let strand = if self.plus { '+' } else { '-' };
        (self.ref_id, self.contig_id, self.pos, strand)
    }

    fn from_py((ref_id, contig_id, pos, strand): PyHit) -> PyResult<Hit> {
        let plus = match strand {
            '+' => true,
            '-' => false,
            _ => {
                return Err(PyValueError::new_err(format!(
                    "Invalid strand '{strand}', expected '+' or '-'"
                )))
            }
        };
        Ok(Hit {
            ref_id,
            contig_id,
            pos,
            plus,
        })
    }

    /// Would `self` (plus) and `other` (minus) amplify a product shorter
    /// than `max_size`?
    fn primes_product_with(&self, other: &Hit, max_size: u32) -> bool {
        self.plus
            && !other.plus
            && self.ref_id == other.ref_id
            && self.contig_id == other.contig_id
            && other.pos > self.pos
            && other.pos - self.pos < max_size
    }
}

fn from_py_hits(hits: Vec<PyHit>) -> PyResult<Vec<Hit>> {
    hits.into_iter().map(Hit::from_py).collect()
}

fn base_code(base: u8) -> Option<u64> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// 2-bit encoding of an A/C/G/T kmer of at most 32 bases
pub fn encode_kmer(seq: &[u8]) -> Option<u64> {
    seq.iter()
        .try_fold(0u64, |acc, &b| Some((acc << 2) | base_code(b)?))
}

/// The encoded `kmer` and every kmer one substitution away from it
pub fn single_mismatches(kmer: u64, k: usize) -> Vec<u64> {
    let mut variants = Vec::with_capacity(3 * k + 1);
    variants.push(kmer);
    for shift in (0..k).map(|i| 2 * i) {
        let current = (kmer >> shift) & 3;
        for alt in (0..4u64).filter(|&alt| alt != current) {
            variants.push(kmer ^ ((current ^ alt) << shift));
        }
    }
    variants
}

/// Reverse complement of an encoded kmer. A base's complement is 3 - its
/// code (A=0, C=1, G=2, T=3).
pub fn reverse_complement_code(code: u64, k: usize) -> u64 {
    (0..k).fold(0u64, |rc, i| (rc << 2) | (3 - ((code >> (2 * i)) & 3)))
}

/// The smaller of an encoded kmer and its reverse complement, and whether the
/// kmer itself is that one
fn canonical(code: u64, rc: u64) -> (u64, bool) {
    if code <= rc {
        (code, true)
    } else {
        (rc, false)
    }
}

fn mix(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

#[derive(Clone, Copy, Default)]
#[repr(align(64))]
struct Block([u64; 8]);

/// A blocked Bloom filter of kmer codes: each key sets 4 bits in one
/// cache-line block, so a lookup reads one line. At 10 bits a key it is a
/// few MB, which stays in cache where the query map does not.
struct KmerFilter {
    blocks: Vec<Block>,
}

impl KmerFilter {
    const BITS_PER_KEY: usize = 10;

    fn new<'a>(keys: impl ExactSizeIterator<Item = &'a u64>) -> KmerFilter {
        let n_blocks = (keys.len() * Self::BITS_PER_KEY).div_ceil(512).max(1);
        let mut filter = KmerFilter {
            blocks: vec![Block::default(); n_blocks],
        };
        for &key in keys {
            let (block, bits) = filter.locate(key);
            for bit in bits {
                filter.blocks[block].0[bit >> 6] |= 1 << (bit & 63);
            }
        }
        filter
    }

    fn locate(&self, key: u64) -> (usize, [usize; 4]) {
        let h = mix(key);
        let block = ((h as u128 * self.blocks.len() as u128) >> 64) as usize;
        let bit = |i: u32| ((h >> (9 * i)) & 511) as usize;
        (block, [bit(0), bit(1), bit(2), bit(3)])
    }

    fn may_contain(&self, key: u64) -> bool {
        let (block, bits) = self.locate(key);
        let words = &self.blocks[block].0;
        bits.iter()
            .all(|&bit| words[bit >> 6] & (1 << (bit & 63)) != 0)
    }
}

/// Number of sequences `seq` expands into, or None if a base can't expand
fn expansion_size(seq: &[u8]) -> Option<usize> {
    seq.iter().try_fold(1usize, |acc, &b| {
        Some(acc.saturating_mul(expand_amb_base(b)?.len()))
    })
}

/// The A/C/G/T sequences `seq` stands for, if there are at most
/// [`MAX_EXPANSION`] of them
fn bounded_expansion(seq: &[u8]) -> Option<Vec<Vec<u8>>> {
    match expansion_size(seq)? {
        n if n <= MAX_EXPANSION => expand_amb_sequence(seq),
        _ => None,
    }
}

/// Query kmers (and their variants) by canonical encoding, each with whether
/// the variant is in its canonical orientation. A reference window in the
/// same orientation puts the primer on the plus strand, the opposite on the
/// minus strand.
pub struct QueryIndex {
    k: usize,
    filter: KmerFilter,
    map: FxHashMap<u64, Vec<(u32, bool)>>,
}

impl QueryIndex {
    /// `queries` must all be `k` A/C/G/T bases
    pub fn new(queries: &[Vec<u8>], k: usize, fuzzy: bool) -> QueryIndex {
        let mut map: FxHashMap<u64, Vec<(u32, bool)>> = FxHashMap::default();
        for (query_id, query) in queries.iter().enumerate() {
            let Some(code) = encode_kmer(query) else {
                continue;
            };
            let variants = if fuzzy {
                single_mismatches(code, k)
            } else {
                vec![code]
            };
            for variant in variants {
                let rc = reverse_complement_code(variant, k);
                let (key, forward) = canonical(variant, rc);
                let entries = map.entry(key).or_default();
                entries.push((query_id as u32, forward));
                // A palindrome binds on both strands
                if variant == rc {
                    entries.push((query_id as u32, !forward));
                }
            }
        }
        QueryIndex {
            k,
            filter: KmerFilter::new(map.keys()),
            map,
        }
    }

    fn lookup(&self, code: u64, rc: u64, start: usize, emit: &mut impl FnMut(u32, usize, bool)) {
        let (key, forward) = canonical(code, rc);
        if !self.filter.may_contain(key) {
            return;
        }
        if let Some(queries) = self.map.get(&key) {
            for &(query_id, query_forward) in queries {
                emit(query_id, start, forward == query_forward);
            }
        }
    }

    /// Call `emit(query, window_start, plus)` for every window of `seq`
    /// matching a query. Windows with ambiguity codes are expanded; windows
    /// with an N or other unexpandable base are skipped.
    pub fn scan(&self, seq: &[u8], emit: &mut impl FnMut(u32, usize, bool)) {
        let k = self.k;
        let mask = if k == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };
        let rc_shift = 2 * (k - 1);
        // The window and its reverse complement, updated a base at a time
        let mut code = 0u64;
        let mut rc = 0u64;
        // Consecutive A/C/G/T bases ending at the current one
        let mut run = 0usize;
        let mut last_unexpandable: Option<usize> = None;

        for (i, &base) in seq.iter().enumerate() {
            match base_code(base) {
                Some(c) => {
                    code = ((code << 2) | c) & mask;
                    rc = (rc >> 2) | ((3 - c) << rc_shift);
                    run += 1;
                }
                None => {
                    run = 0;
                    if !is_amb(&base) {
                        last_unexpandable = Some(i);
                    }
                }
            }
            if i + 1 < k {
                continue;
            }
            let start = i + 1 - k;
            if run >= k {
                self.lookup(code, rc, start, emit);
            } else if last_unexpandable.is_none_or(|j| j < start) {
                for expanded in bounded_expansion(&seq[start..=i]).unwrap_or_default() {
                    if let Some(c) = encode_kmer(&expanded) {
                        self.lookup(c, reverse_complement_code(c, k), start, emit);
                    }
                }
            }
        }
    }
}

/// Matches as (query, hit), deduplicated
type QueryHits = FxHashSet<(u32, Hit)>;

/// An MSA's distinct rows and, per column, the position hits there report:
/// the reference coordinate when remapped (the nearest preceding one for a
/// column the reference has a gap in), else the column.
struct MsaRef {
    rows: Vec<Vec<u8>>,
    col_pos: Vec<u32>,
}

impl MsaRef {
    fn load(path: &str, remap: bool) -> Result<MsaRef, String> {
        let (_headers, seqs) = try_fasta_reader(path)?;
        let Some(first) = seqs.first() else {
            return Err(format!("{path}: no sequences in the MSA"));
        };
        let mut last = 0u32;
        let col_pos = create_mapping_array(first.as_bytes(), remap)
            .into_iter()
            .map(|mapped| {
                if let Some(pos) = mapped {
                    last = pos as u32;
                }
                last
            })
            .collect();
        let mut seen = FxHashSet::default();
        let rows = seqs
            .into_iter()
            .map(String::into_bytes)
            .filter(|row| seen.insert(row.clone()))
            .collect();
        Ok(MsaRef { rows, col_pos })
    }

    fn scan(&self, ref_id: u32, index: &QueryIndex) -> QueryHits {
        let k = index.k;
        self.rows
            .par_iter()
            .fold(QueryHits::default, |mut found, row| {
                let (bases, cols): (Vec<u8>, Vec<usize>) = row
                    .iter()
                    .enumerate()
                    .filter(|(_, &b)| b != b'-' && b != b' ')
                    .map(|(col, &b)| (b, col))
                    .unzip();
                index.scan(&bases, &mut |query_id, start, plus| {
                    let three_prime = if plus { start + k - 1 } else { start };
                    found.insert((
                        query_id,
                        Hit {
                            ref_id,
                            contig_id: 0,
                            pos: self.col_pos[cols[three_prime]],
                            plus,
                        },
                    ));
                });
                found
            })
            .reduce(QueryHits::default, merge)
    }
}

fn merge(mut a: QueryHits, b: QueryHits) -> QueryHits {
    if a.len() < b.len() {
        return merge(b, a);
    }
    a.extend(b);
    a
}

fn open_fasta(path: &str) -> std::io::Result<Box<dyn BufRead>> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::with_capacity(
            READ_BUFFER,
            MultiGzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(READ_BUFFER, file)))
    }
}

/// Call `f` with each record's name (the header's first word) and
/// uppercased sequence, reading one record at a time. Plain or gzipped.
/// Gaps (`-`, `.`) and whitespace are dropped, so an aligned FASTA reads as
/// its ungapped genomes.
pub fn for_each_record(path: &str, mut f: impl FnMut(String, Vec<u8>)) -> Result<(), String> {
    let read_err = |e: std::io::Error| format!("Failed to read {path}: {e}");
    let mut reader = open_fasta(path).map_err(read_err)?;
    let mut name: Option<String> = None;
    let mut seq = Vec::new();
    let mut buf = Vec::new();
    loop {
        buf.clear();
        if reader.read_until(b'\n', &mut buf).map_err(read_err)? == 0 {
            break;
        }
        let line = buf.trim_ascii_end();
        if let Some(header) = line.strip_prefix(b">") {
            if let Some(prev) = name.take() {
                f(prev, std::mem::take(&mut seq));
            }
            let header = String::from_utf8_lossy(header);
            name = Some(header.split_whitespace().next().unwrap_or("").to_string());
        } else if name.is_some() {
            push_bases(&mut seq, line);
        } else if !line.trim_ascii().is_empty() {
            return Err(format!("{path}: sequence before the first FASTA header"));
        }
    }
    if let Some(prev) = name {
        f(prev, seq);
    }
    Ok(())
}

fn is_gap_or_space(base: u8) -> bool {
    matches!(base, b'-' | b'.') || base.is_ascii_whitespace()
}

/// Append a sequence line's bases, uppercased, without gaps or whitespace
fn push_bases(seq: &mut Vec<u8>, line: &[u8]) {
    // A non-short-circuiting fold, so the common clean line is vectorised
    if line.iter().fold(false, |any, &b| any | is_gap_or_space(b)) {
        seq.extend(
            line.iter()
                .filter(|&&b| !is_gap_or_space(b))
                .map(u8::to_ascii_uppercase),
        );
    } else {
        let start = seq.len();
        seq.extend_from_slice(line);
        seq[start..].make_ascii_uppercase();
    }
}

/// Every match in a contig, scanned in parallel chunks. Chunks overlap by
/// k - 1 bases so each window falls in exactly one.
fn scan_contig(seq: &[u8], ref_id: u32, contig_id: u32, index: &QueryIndex) -> QueryHits {
    let k = index.k;
    (0..seq.len())
        .step_by(CHUNK)
        .collect::<Vec<_>>()
        .into_par_iter()
        .fold(QueryHits::default, |mut found, chunk_start| {
            let chunk_end = (chunk_start + CHUNK + k - 1).min(seq.len());
            index.scan(
                &seq[chunk_start..chunk_end],
                &mut |query_id, start, plus| {
                    let three_prime = chunk_start + start + if plus { k - 1 } else { 0 };
                    found.insert((
                        query_id,
                        Hit {
                            ref_id,
                            contig_id,
                            pos: three_prime as u32,
                            plus,
                        },
                    ));
                },
            );
            found
        })
        .reduce(QueryHits::default, merge)
}

/// Every match in a background FASTA, and its contig names in file order.
/// A reader thread reads the next contig while this one is scanned.
fn scan_background(
    path: &str,
    ref_id: u32,
    index: &QueryIndex,
) -> Result<(QueryHits, Vec<String>), String> {
    let (send, contigs) = std::sync::mpsc::sync_channel::<(String, Vec<u8>)>(1);
    std::thread::scope(|scope| {
        let reader = scope.spawn(move || {
            for_each_record(path, |name, seq| {
                // Only fails once the scan below has stopped receiving
                let _ = send.send((name, seq));
            })
        });
        let mut found = QueryHits::default();
        let mut names = Vec::new();
        for (name, seq) in contigs {
            found = merge(found, scan_contig(&seq, ref_id, names.len() as u32, index));
            names.push(name);
        }
        reader.join().expect("FASTA reader panicked")?;
        Ok((found, names))
    })
}

/// The references primers are screened against: MSAs first (ref ids
/// 0..n_msas, positions in the digest's coordinates, one contig each), then
/// background genomes (read from disk on each screen, one contig per record).
/// A primer's hits are cached by its 3' terminal kmer.
#[pyclass]
pub struct MisprimeDB {
    k: usize,
    fuzzy: bool,
    msas: Vec<MsaRef>,
    backgrounds: Vec<String>,
    ref_names: Vec<String>,
    contig_names: Vec<Vec<String>>,
    cache: FxHashMap<Vec<u8>, Vec<Hit>>,
}

impl MisprimeDB {
    pub fn new(
        msa_paths: Vec<String>,
        background_paths: Vec<String>,
        kmer_size: usize,
        fuzzy: bool,
        remap: bool,
    ) -> Result<MisprimeDB, String> {
        if !(1..=32).contains(&kmer_size) {
            return Err(format!("kmer_size must be 1-32, not {kmer_size}"));
        }
        let msas = msa_paths
            .iter()
            .map(|p| MsaRef::load(p, remap))
            .collect::<Result<Vec<_>, _>>()?;
        for path in &background_paths {
            File::open(path).map_err(|e| format!("Failed to read {path}: {e}"))?;
        }
        let ref_names = msa_paths.iter().chain(&background_paths).cloned().collect();
        let contig_names = msa_paths
            .iter()
            .map(|_| vec![String::new()])
            .chain(background_paths.iter().map(|_| Vec::new()))
            .collect();
        Ok(MisprimeDB {
            k: kmer_size,
            fuzzy,
            msas,
            backgrounds: background_paths,
            ref_names,
            contig_names,
            cache: FxHashMap::default(),
        })
    }

    /// A primer's uppercased 3' terminal kmer; None if it is shorter than
    /// the kmer, as it can't be screened
    fn terminal(&self, seq: &[u8]) -> Option<Vec<u8>> {
        let start = seq.len().checked_sub(self.k)?;
        Some(seq[start..].to_ascii_uppercase())
    }

    /// Screen every primer whose terminal kmer isn't cached yet, in one pass
    /// over the references. Returns how many terminal kmers were screened.
    pub fn screen(&mut self, seqs: &[Vec<u8>]) -> Result<usize, String> {
        let mut terminals: Vec<Vec<u8>> = Vec::new();
        let mut seen = FxHashSet::default();
        for terminal in seqs.iter().filter_map(|seq| self.terminal(seq)) {
            if !self.cache.contains_key(&terminal) && seen.insert(terminal.clone()) {
                terminals.push(terminal);
            }
        }
        if terminals.is_empty() {
            return Ok(0);
        }
        let mut queries = Vec::new();
        let mut owners = Vec::new();
        for (owner, terminal) in terminals.iter().enumerate() {
            for query in bounded_expansion(terminal).unwrap_or_default() {
                queries.push(query);
                owners.push(owner);
            }
        }
        let index = QueryIndex::new(&queries, self.k, self.fuzzy);

        let mut found = QueryHits::default();
        for (ref_id, msa) in self.msas.iter().enumerate() {
            found = merge(found, msa.scan(ref_id as u32, &index));
        }
        for (i, path) in self.backgrounds.iter().enumerate() {
            let ref_id = self.msas.len() + i;
            let (hits, names) = scan_background(path, ref_id as u32, &index)?;
            found = merge(found, hits);
            self.contig_names[ref_id] = names;
        }

        let mut per_terminal: Vec<Vec<Hit>> = vec![Vec::new(); terminals.len()];
        for (query_id, hit) in found {
            per_terminal[owners[query_id as usize]].push(hit);
        }
        let screened = terminals.len();
        for (terminal, mut hits) in terminals.into_iter().zip(per_terminal) {
            hits.sort_unstable();
            hits.dedup();
            self.cache.insert(terminal, hits);
        }
        Ok(screened)
    }

    /// Where a primer's 3' terminal kmer binds, screening it if needed.
    /// Empty for a primer shorter than the kmer.
    pub fn hits(&mut self, seq: &[u8]) -> Result<&[Hit], String> {
        self.screen(&[seq.to_vec()])?;
        Ok(self
            .terminal(seq)
            .map_or(&[], |terminal| self.cache[&terminal].as_slice()))
    }
}

#[pymethods]
impl MisprimeDB {
    /// Load the MSAs (with the digest's `remap` setting, so positions match
    /// the kmers') and check the background FASTAs (plain or gzipped) can be
    /// opened. ValueError if a file can't be read or kmer_size isn't 1-32.
    #[new]
    #[pyo3(signature = (msa_paths, background_paths, kmer_size, fuzzy, remap))]
    fn py_new(
        msa_paths: Vec<String>,
        background_paths: Vec<String>,
        kmer_size: usize,
        fuzzy: bool,
        remap: bool,
    ) -> PyResult<MisprimeDB> {
        MisprimeDB::new(msa_paths, background_paths, kmer_size, fuzzy, remap)
            .map_err(PyValueError::new_err)
    }

    /// Screen these primers' 3' terminal kmers in one pass over the
    /// references. Returns how many weren't already cached.
    #[pyo3(name = "screen")]
    fn py_screen(&mut self, py: Python<'_>, seqs: Vec<Vec<u8>>) -> PyResult<usize> {
        py.detach(|| self.screen(&seqs))
            .map_err(PyValueError::new_err)
    }

    /// Hits of each primer's 3' terminal kmer, as (ref_id, contig_id, pos,
    /// strand) with pos the primer's 3' end base. Unscreened primers are
    /// screened first, in one pass. Primers shorter than the kmer have none.
    #[pyo3(name = "hits")]
    fn py_hits(&mut self, py: Python<'_>, seqs: Vec<Vec<u8>>) -> PyResult<Vec<Vec<PyHit>>> {
        py.detach(|| self.screen(&seqs))
            .map_err(PyValueError::new_err)?;
        seqs.iter()
            .map(|seq| {
                let hits = self.hits(seq).map_err(PyValueError::new_err)?;
                Ok(hits.iter().map(|h| h.to_py()).collect())
            })
            .collect()
    }

    #[getter]
    fn kmer_size(&self) -> usize {
        self.k
    }

    #[getter]
    fn n_msas(&self) -> usize {
        self.msas.len()
    }

    /// Each reference's path, MSAs first
    fn ref_names(&self) -> Vec<String> {
        self.ref_names.clone()
    }

    /// A background's contig names, known once it has been screened
    fn contig_names(&self, ref_id: usize) -> Vec<String> {
        self.contig_names.get(ref_id).cloned().unwrap_or_default()
    }

    /// Distinct binding sites per reference over every cached primer
    fn site_counts(&self) -> Vec<usize> {
        let sites: FxHashSet<&Hit> = self.cache.values().flatten().collect();
        let mut counts = vec![0; self.ref_names.len()];
        for hit in sites {
            counts[hit.ref_id as usize] += 1;
        }
        counts
    }

    fn __len__(&self) -> usize {
        self.cache.len()
    }
}

#[derive(Default)]
struct Sites {
    plus: Vec<u32>,
    minus: Vec<u32>,
}

/// The binding sites of a pool's primers, as sorted positions per reference
/// contig and strand. A multiset: sites shared by primers are kept until
/// every primer holding them is removed.
#[pyclass]
#[derive(Default)]
pub struct PoolHits {
    sites: FxHashMap<(u32, u32), Sites>,
    len: usize,
}

impl PoolHits {
    pub fn add(&mut self, hits: &[Hit]) {
        for hit in hits {
            let sites = self.sites.entry((hit.ref_id, hit.contig_id)).or_default();
            let positions = if hit.plus {
                &mut sites.plus
            } else {
                &mut sites.minus
            };
            let i = positions.partition_point(|&p| p < hit.pos);
            positions.insert(i, hit.pos);
            self.len += 1;
        }
    }

    pub fn remove(&mut self, hits: &[Hit]) {
        for hit in hits {
            let Some(sites) = self.sites.get_mut(&(hit.ref_id, hit.contig_id)) else {
                continue;
            };
            let positions = if hit.plus {
                &mut sites.plus
            } else {
                &mut sites.minus
            };
            if let Ok(i) = positions.binary_search(&hit.pos) {
                positions.remove(i);
                self.len -= 1;
            }
        }
    }

    /// Does `hit` amplify a product shorter than `max_size` with a site
    /// already in the pool?
    fn primes_product(&self, hit: &Hit, max_size: u32) -> bool {
        let Some(sites) = self.sites.get(&(hit.ref_id, hit.contig_id)) else {
            return false;
        };
        if hit.plus {
            let i = sites.minus.partition_point(|&p| p <= hit.pos);
            sites.minus.get(i).is_some_and(|&p| p - hit.pos < max_size)
        } else {
            let i = sites.plus.partition_point(|&p| p < hit.pos);
            i > 0 && hit.pos - sites.plus[i - 1] < max_size
        }
    }

    /// The first product shorter than `max_size` adding `hits` would give:
    /// its reference, and whether the hits form it among themselves (`alone`)
    /// rather than with the pool. Two `expected` hits (the primers' intended
    /// sites) are allowed to pair with each other.
    pub fn find_product(
        &self,
        hits: &[Hit],
        max_size: u32,
        expected: &[Hit],
    ) -> Option<(u32, bool)> {
        if let Some(hit) = hits.iter().find(|h| self.primes_product(h, max_size)) {
            return Some((hit.ref_id, false));
        }
        hits.iter()
            .find(|f| {
                hits.iter().any(|r| {
                    f.primes_product_with(r, max_size)
                        && !(expected.contains(f) && expected.contains(r))
                })
            })
            .map(|f| (f.ref_id, true))
    }

    pub fn forms_product(&self, hits: &[Hit], max_size: u32, expected: &[Hit]) -> bool {
        self.find_product(hits, max_size, expected).is_some()
    }
}

#[pymethods]
impl PoolHits {
    #[new]
    fn py_new() -> PoolHits {
        PoolHits::default()
    }

    #[pyo3(name = "add")]
    fn py_add(&mut self, hits: Vec<PyHit>) -> PyResult<()> {
        self.add(&from_py_hits(hits)?);
        Ok(())
    }

    /// Remove one copy of each hit; absent hits are ignored
    #[pyo3(name = "remove")]
    fn py_remove(&mut self, hits: Vec<PyHit>) -> PyResult<()> {
        self.remove(&from_py_hits(hits)?);
        Ok(())
    }

    #[pyo3(name = "forms_product", signature = (hits, max_size, expected=Vec::new()))]
    fn py_forms_product(
        &self,
        hits: Vec<PyHit>,
        max_size: u32,
        expected: Vec<PyHit>,
    ) -> PyResult<bool> {
        Ok(self.forms_product(&from_py_hits(hits)?, max_size, &from_py_hits(expected)?))
    }

    /// The first product the hits would give, as (ref_id, alone), or None
    #[pyo3(name = "find_product", signature = (hits, max_size, expected=Vec::new()))]
    fn py_find_product(
        &self,
        hits: Vec<PyHit>,
        max_size: u32,
        expected: Vec<PyHit>,
    ) -> PyResult<Option<(u32, bool)>> {
        Ok(self.find_product(&from_py_hits(hits)?, max_size, &from_py_hits(expected)?))
    }

    fn __len__(&self) -> usize {
        self.len
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seqfuncs::reverse_complement;
    use std::io::Write;

    fn hit(ref_id: u32, pos: u32, plus: bool) -> Hit {
        Hit {
            ref_id,
            contig_id: 0,
            pos,
            plus,
        }
    }

    fn write_tmp(name: &str, contents: &str) -> String {
        let path = std::env::temp_dir().join(format!(
            "primalschemers-mispriming-{}-{name}",
            std::process::id()
        ));
        File::create(&path)
            .unwrap()
            .write_all(contents.as_bytes())
            .unwrap();
        path.to_string_lossy().into_owned()
    }

    fn scan_all(seq: &[u8], queries: &[&[u8]], fuzzy: bool) -> Vec<(u32, usize, bool)> {
        let queries: Vec<Vec<u8>> = queries.iter().map(|q| q.to_vec()).collect();
        let index = QueryIndex::new(&queries, queries[0].len(), fuzzy);
        let mut found = Vec::new();
        index.scan(seq, &mut |q, start, plus| found.push((q, start, plus)));
        found.sort();
        found.dedup();
        found
    }

    #[test]
    fn test_encode_kmer() {
        assert_eq!(encode_kmer(b"A"), Some(0));
        assert_eq!(encode_kmer(b"ACGT"), Some(0b00_01_10_11));
        assert_eq!(encode_kmer(b"ACNT"), None);
        assert_eq!(encode_kmer(&[b'T'; 32]), Some(u64::MAX));
    }

    #[test]
    fn test_reverse_complement_code() {
        for seq in [&b"ACGTTGCA"[..], b"AAAAC", b"GATTACA", &[b'G'; 32]] {
            assert_eq!(
                reverse_complement_code(encode_kmer(seq).unwrap(), seq.len()),
                encode_kmer(&reverse_complement(seq)).unwrap()
            );
        }
    }

    #[test]
    fn test_kmer_filter_has_no_false_negatives() {
        let keys: Vec<u64> = (0..100_000u64).map(|i| mix(i) >> 26).collect();
        let filter = KmerFilter::new(keys.iter());
        assert!(keys.iter().all(|&key| filter.may_contain(key)));
        let false_positives = (0..100_000u64)
            .map(|i| mix(i + (1 << 40)) >> 26)
            .filter(|key| filter.may_contain(*key))
            .count();
        assert!(false_positives < 2_000, "{false_positives}");
    }

    #[test]
    fn test_scan_palindrome_both_strands() {
        // ACGT is its own reverse complement
        assert_eq!(
            scan_all(b"TTACGTTT", &[b"ACGT"], false),
            vec![(0, 2, false), (0, 2, true)]
        );
    }

    #[test]
    fn test_single_mismatches() {
        let k = 4;
        let variants = single_mismatches(encode_kmer(b"ACGT").unwrap(), k);
        assert_eq!(variants.len(), 3 * k + 1);
        let unique: FxHashSet<u64> = variants.iter().copied().collect();
        assert_eq!(unique.len(), variants.len());
        assert!(unique.contains(&encode_kmer(b"TCGT").unwrap()));
        assert!(unique.contains(&encode_kmer(b"ACGA").unwrap()));
        assert!(!unique.contains(&encode_kmer(b"TCGA").unwrap()));
    }

    #[test]
    fn test_scan_strands() {
        //         0123456789
        let seq = b"TTACGGATTT";
        // The kmer itself: a plus primer
        assert_eq!(scan_all(seq, &[b"ACGGA"], false), vec![(0, 2, true)]);
        // Its reverse complement: a minus primer
        assert_eq!(scan_all(seq, &[b"TCCGT"], false), vec![(0, 2, false)]);
    }

    #[test]
    fn test_scan_fuzzy() {
        let seq = b"TTACGGATTT";
        assert!(scan_all(seq, &[b"ACCGA"], false).is_empty());
        assert_eq!(scan_all(seq, &[b"ACCGA"], true), vec![(0, 2, true)]);
        assert!(scan_all(seq, &[b"GGGGG"], true).is_empty());
    }

    #[test]
    fn test_scan_ambiguity_and_n() {
        // R expands to A/G, so the window matches
        assert_eq!(
            scan_all(b"TTRCGGATTT", &[b"ACGGA"], false),
            vec![(0, 2, true)]
        );
        // N never matches, and windows past it still do
        assert_eq!(
            scan_all(b"TTNCGGACGGA", &[b"ACGGA"], false),
            vec![(0, 6, true)]
        );
        // The rolling code recovers after an ambiguity
        assert_eq!(
            scan_all(b"RACGGATT", &[b"ACGGA"], false),
            vec![(0, 1, true)]
        );
    }

    #[test]
    fn test_msa_positions_across_gaps() {
        // The reference (row 0) has a gap at column 4 and row 1 doesn't match
        let path = write_tmp("msa.fasta", ">ref\nAACG-GATTT\n>alt\nAAC--GATTT\n");
        let index = QueryIndex::new(&[b"ACGGA".to_vec()], 5, false);
        let msa = MsaRef::load(&path, true).unwrap();
        let hits: Vec<Hit> = msa.scan(7, &index).into_iter().map(|(_, h)| h).collect();
        // ACGGA ends on column 6, reference position 5
        assert_eq!(hits, vec![hit(7, 5, true)]);

        let msa = MsaRef::load(&path, false).unwrap();
        let hits: Vec<Hit> = msa.scan(7, &index).into_iter().map(|(_, h)| h).collect();
        assert_eq!(hits, vec![hit(7, 6, true)]);
    }

    #[test]
    fn test_background_contigs_and_chunks() {
        let mut big = "T".repeat(CHUNK - 2);
        big.push_str("ACGGA"); // straddles the chunk boundary
        big.push_str(&"T".repeat(10));
        let path = write_tmp(
            "bg.fasta",
            &format!(">one desc\nTTACG\nGATT\n>two\n{big}\n"),
        );
        let index = QueryIndex::new(&[b"ACGGA".to_vec()], 5, false);
        let (found, names) = scan_background(&path, 3, &index).unwrap();
        let mut hits: Vec<Hit> = found.into_iter().map(|(_, h)| h).collect();
        hits.sort();
        assert_eq!(names, vec!["one", "two"]);
        assert_eq!(
            hits,
            vec![
                Hit {
                    ref_id: 3,
                    contig_id: 0,
                    pos: 6,
                    plus: true
                },
                Hit {
                    ref_id: 3,
                    contig_id: 1,
                    pos: (CHUNK + 2) as u32,
                    plus: true
                },
            ]
        );
    }

    #[test]
    fn test_gapped_background() {
        let path = write_tmp("gapped.fasta", ">aligned\nTT-AC.G GA\n--TTT\n");
        let index = QueryIndex::new(&[b"ACGGA".to_vec()], 5, false);
        let (found, _) = scan_background(&path, 0, &index).unwrap();
        let hits: Vec<Hit> = found.into_iter().map(|(_, h)| h).collect();
        // Found across the gaps, at its ungapped position
        assert_eq!(hits, vec![hit(0, 6, true)]);
    }

    #[test]
    fn test_gzipped_background() {
        use flate2::write::GzEncoder;
        let path = std::env::temp_dir().join(format!(
            "primalschemers-mispriming-{}-bg.fasta.gz",
            std::process::id()
        ));
        let mut gz = GzEncoder::new(File::create(&path).unwrap(), flate2::Compression::fast());
        gz.write_all(b">one\nttacggattt\n").unwrap();
        gz.finish().unwrap();
        let mut records = Vec::new();
        for_each_record(path.to_str().unwrap(), |name, seq| {
            records.push((name, seq))
        })
        .unwrap();
        assert_eq!(records, vec![("one".to_string(), b"TTACGGATTT".to_vec())]);
    }

    #[test]
    fn test_db_screens_terminal_kmer_and_caches() {
        let msa = write_tmp("db-msa.fasta", ">ref\nTTACGGATTTTTTTGGGCCCTT\n");
        let bg = write_tmp("db-bg.fasta", ">chr\nCCCCCACGGACCCCC\n");
        let mut db = MisprimeDB::new(vec![msa], vec![bg], 5, false, true).unwrap();
        // Only the last 5 bases count
        let primer = b"GGGGGACGGA".to_vec();
        assert_eq!(
            db.hits(&primer).unwrap(),
            &[
                hit(0, 6, true),
                Hit {
                    ref_id: 1,
                    contig_id: 0,
                    pos: 9,
                    plus: true
                },
            ]
        );
        assert_eq!(db.screen(&[primer]).unwrap(), 0);
        // Too short to screen
        assert!(db.hits(b"ACGG").unwrap().is_empty());
    }

    #[test]
    fn test_pool_hits_products() {
        let mut pool = PoolHits::default();
        pool.add(&[hit(0, 100, true), hit(0, 1000, false)]);
        // A minus site 100 after the pool's plus
        assert!(pool.forms_product(&[hit(0, 200, false)], 500, &[]));
        // Too far
        assert!(!pool.forms_product(&[hit(0, 700, false)], 500, &[]));
        // Upstream of the plus site: no product
        assert!(!pool.forms_product(&[hit(0, 99, false)], 500, &[]));
        // A plus site 100 before the pool's minus
        assert!(pool.forms_product(&[hit(0, 900, true)], 500, &[]));
        // Same position on another reference or contig
        assert!(!pool.forms_product(&[hit(1, 200, false)], 500, &[]));
        let other_contig = Hit {
            contig_id: 1,
            ..hit(0, 200, false)
        };
        assert!(!pool.forms_product(&[other_contig], 500, &[]));
    }

    #[test]
    fn test_pool_hits_new_pair_with_itself() {
        let pool = PoolHits::default();
        let (f, r) = (hit(0, 100, true), hit(0, 300, false));
        assert!(pool.forms_product(&[f, r], 500, &[]));
        // The pair's own intended sites are its amplicon
        assert!(!pool.forms_product(&[f, r], 500, &[f, r]));
        // An off-target site pairs with an intended one
        let off = hit(0, 250, true);
        assert!(pool.forms_product(&[f, r, off], 500, &[f, r]));
        assert_eq!(
            pool.find_product(&[f, r, off], 500, &[f, r]),
            Some((0, true))
        );
    }

    #[test]
    fn test_find_product_reports_reference() {
        let mut pool = PoolHits::default();
        pool.add(&[hit(2, 100, true)]);
        assert_eq!(
            pool.find_product(&[hit(2, 200, false)], 500, &[]),
            Some((2, false))
        );
        assert_eq!(pool.find_product(&[hit(1, 200, false)], 500, &[]), None);
    }

    #[test]
    fn test_pool_hits_multiset_remove() {
        let mut pool = PoolHits::default();
        let shared = hit(0, 100, true);
        pool.add(&[shared]);
        pool.add(&[shared]);
        pool.remove(&[shared]);
        assert_eq!(pool.len, 1);
        assert!(pool.forms_product(&[hit(0, 200, false)], 500, &[]));
        pool.remove(&[shared, shared]);
        assert_eq!(pool.len, 0);
        assert!(!pool.forms_product(&[hit(0, 200, false)], 500, &[]));
    }
}
