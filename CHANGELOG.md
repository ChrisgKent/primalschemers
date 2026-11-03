# Changelog

## 0.3.2

- New `MisprimeDB` and `PoolHits` for off-target (mispriming) screening of primers' 3' terminal kmers, replacing primalscheme3's Python `MatchDB`. The primers are indexed (with their single-mismatch variants when `fuzzy`) and each reference is streamed past them once, so background genomes (plain or gzipped FASTA, any size; gaps are dropped, so an alignment reads as its genomes) can be screened as well as the MSAs. `PoolHits` holds a pool's binding sites as a multiset and finds products by binary search, including between a new pair's own hits. `PoolHits.find_product` also says which reference a product is on and whether the pair forms it alone.

- Tm mode: the digest walk now extends to the midpoint of the Tm window (`(primer_tm_min + primer_tm_max) / 2`) instead of stopping at `primer_tm_min`, and keeps the length closest to it, as annealing mode does with `primer_annealing_prop`.
- A walk that ends before its target (`primer_len_max`, `max_walk`, the alignment edge or an N) keeps the sequence if it is still in the window (Tm ≥ `primer_tm_min`, or annealing ≥ `primer_annealing_prop - primer_annealing_diff`), so no site that passed before is lost.

## 0.3.1

- Fixed: the digest walk could return a primer one base longer than `primer_len_max`. Sites that passed only with that primer are now rejected with `ToLong`.
- `return_all`: a walk that passes `primer_len_max` keeps extending until it reaches the Tm/annealing target. That sequence is returned in `seqs()`, tagged `ToLong` in `statuses()`, instead of as a bare `ToLong` error. A normal digest is unchanged.
