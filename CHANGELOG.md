# Changelog

## 0.3.1

- Fixed: the digest walk could return a primer one base longer than `primer_len_max`. Sites that passed only with that primer are now rejected with `ToLong`.
- `return_all`: a walk that passes `primer_len_max` keeps extending until it reaches the Tm/annealing target. That sequence is returned in `seqs()`, tagged `ToLong` in `statuses()`, instead of as a bare `ToLong` error. A normal digest is unchanged.
