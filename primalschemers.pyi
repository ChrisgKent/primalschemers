from typing import Any, Dict, List, Optional, Tuple, Union

class FKmer:
    def __init__(
        self, seqs: List[bytes], end: int, counts: Optional[List[float]] = None
    ) -> None:
        """ValueError if counts and seqs differ in length."""
    def starts(self) -> List[int]: ...
    @property
    def end(self) -> int: ...
    def counts(self) -> Optional[List[float]]: ...
    def len(self) -> List[int]: ...
    def num_seqs(self) -> int: ...
    def seqs(self) -> List[str]: ...
    def seqs_bytes(self) -> List[bytes]: ...
    def region(self) -> Tuple[int, int]: ...
    def to_bed(self, chrom: str, amplicon_prefix: str, pool: int) -> str: ...
    def remap(self, end: int) -> None: ...
    def gc(self) -> List[float]: ...
    def counts_sum(self) -> Optional[float]: ...
    def statuses(self) -> Optional[List[List[str]]]: ...
    def kmer_statuses(self) -> List[str]: ...
    def errors(self) -> List[Tuple[List[str], float]]: ...
    def dimers(self) -> List[Tuple[int, int]]: ...
    def passed(self) -> bool: ...

class RKmer:
    def __init__(
        self, seqs: List[bytes], start: int, counts: Optional[List[float]] = None
    ) -> None:
        """ValueError if counts and seqs differ in length."""
    def seqs(self) -> List[str]: ...
    @property
    def start(self) -> int: ...
    def counts(self) -> Optional[List[float]]: ...
    def ends(self) -> List[int]: ...
    def lens(self) -> List[int]: ...
    def num_seqs(self) -> int: ...
    def seqs_bytes(self) -> List[bytes]: ...
    def region(self) -> Tuple[int, int]: ...
    def to_bed(self, chrom: str, amplicon_prefix: str, pool: int) -> str: ...
    def remap(self, start: int) -> None: ...
    def gc(self) -> List[float]: ...
    def counts_sum(self) -> Optional[float]: ...
    def statuses(self) -> Optional[List[List[str]]]: ...
    def kmer_statuses(self) -> List[str]: ...
    def errors(self) -> List[Tuple[List[str], float]]: ...
    def dimers(self) -> List[Tuple[int, int]]: ...
    def passed(self) -> bool: ...

class Digester:
    def __init__(self, msa_path: str, ncores: int, remap: bool) -> None:
        """ValueError for an unreadable or empty MSA, or sequences of
        different lengths."""
    def create_seq_slice(self) -> List[bytes]: ...
    def digest(
        self,
        findexes: Optional[List[int]] = None,
        rindexes: Optional[List[int]] = None,
        primer_len_min: Optional[int] = None,
        primer_len_max: Optional[int] = None,
        primer_gc_max: Optional[float] = None,
        primer_gc_min: Optional[float] = None,
        primer_tm_max: Optional[float] = None,
        primer_tm_min: Optional[float] = None,
        primer_annealing_prop: Optional[float] = None,
        primer_annealing_diff: Optional[float] = None,
        annealing_temp_c: Optional[float] = None,
        max_walk: Optional[int] = None,
        max_homopolymers: Optional[int] = None,
        min_freq: Optional[float] = None,
        ignore_n: Optional[bool] = None,
        dimerscore: Optional[float] = None,
        thermo_check: bool = False,
        primer_hairpin_th_max: Optional[float] = None,
        return_all: bool = False,
        mv_conc: Optional[float] = None,
        dv_conc: Optional[float] = None,
        dntp_conc: Optional[float] = None,
        dna_conc: Optional[float] = None,
        formamide_conc: Optional[float] = None,
) -> Tuple[List[FKmer], List[RKmer], List[str]]:
        """Digest the MSA (releases the GIL). Concentrations (mM; DNA nM) and
        formamide (mol/l) set the Tm/annealing conditions, defaulting to
        primalscheme3's (100, 2, 0.8, 15, formamide 0.8). The hairpin check
        uses the same salts without formamide, as primer3's calc_hairpin.
        Annealing mode is used whenever primer_annealing_prop is set."""
    @property
    def _seq_array(self) -> List[bytes]: ...
    @property
    def _mapping_array(self) -> List[Optional[int]]: ...
    @property
    def _ref_to_msa_array(self) -> List[int]: ...

def digest_seq(
    msa_path: str,
    ncores: int,
    remap: bool,
    findexes: Optional[List[int]] = None,
    rindexes: Optional[List[int]] = None,
    primer_len_min: Optional[int] = None,
    primer_len_max: Optional[int] = None,
    primer_gc_max: Optional[float] = None,
    primer_gc_min: Optional[float] = None,
    primer_tm_max: Optional[float] = None,
    primer_tm_min: Optional[float] = None,
    primer_annealing_prop: Optional[float] = None,
    primer_annealing_diff: Optional[float] = None,
    annealing_temp_c: Optional[float] = None,
    max_walk: Optional[int] = None,
    max_homopolymers: Optional[int] = None,
    min_freq: Optional[float] = None,
    ignore_n: Optional[bool] = None,
    dimerscore: Optional[float] = None,
    thermo_check: bool = True,
    primer_hairpin_th_max: Optional[float] = None,
    return_all: bool = False,
    mv_conc: Optional[float] = None,
    dv_conc: Optional[float] = None,
    dntp_conc: Optional[float] = None,
    dna_conc: Optional[float] = None,
    formamide_conc: Optional[float] = None,
) -> Tuple[List[FKmer], List[RKmer], List[str]]:
    """Digester(msa_path, ncores, remap).digest(...), with thermo_check on
    by default."""
def do_seqs_interact(seq1: bytes, seq2: bytes, t: float) -> bool: ...
def do_pool_interact(seqs1: List[bytes], seqs2: List[bytes], t: float) -> bool: ...
def calc_at_offset_py(seq1: str, seq2: str, offset: int) -> float: ...
def calc_annealing(
    sequence: str,
    dna_nm: float,
    k_mm: float,
    divalent_conc: float,
    dntp_conc: float,
    dmso_conc: float,
    dmso_fact: float,
    formamide_conc: float,
    annealing_temp_c: float,
) -> float:
    """% of sequence bound to its complement. DMSO/formamide are accepted
    but don't affect binding. ValueError for non-ACGT bases."""
def generate_primerpairs_py(
    fkmers: List[FKmer],
    rkmers: List[RKmer],
    t: float,
    amplicon_size_min: int,
    amplicon_size_max: int,
) -> List[Tuple[FKmer, RKmer]]: ...

# Thermodynamics via primer3-thal-rs: drop-in for primer3-py 2.3.0's
# calc_hairpin / calc_homodimer / calc_heterodimer.

class ThermoResult:
    @property
    def structure_found(self) -> bool: ...
    @property
    def tm(self) -> float: ...
    @property
    def dg(self) -> float: ...
    @property
    def dh(self) -> float: ...
    @property
    def ds(self) -> float: ...
    @property
    def ascii_structure(self) -> Optional[str]: ...
    @property
    def ascii_structure_lines(self) -> Optional[List[str]]: ...
    def todict(self) -> Dict[str, Any]: ...

def calc_hairpin(
    seq: Union[str, bytes],
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dna_conc: float = 50.0,
    temp_c: float = 37.0,
    max_loop: int = 30,
    output_structure: bool = False,
    dmso_conc: float = 0.0,
    dmso_fact: float = 0.6,
    formamide_conc: float = 0.0,
) -> ThermoResult: ...
def calc_homodimer(
    seq: Union[str, bytes],
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dna_conc: float = 50.0,
    temp_c: float = 37.0,
    max_loop: int = 30,
    output_structure: bool = False,
    dmso_conc: float = 0.0,
    dmso_fact: float = 0.6,
    formamide_conc: float = 0.0,
) -> ThermoResult: ...
def calc_heterodimer(
    seq1: Union[str, bytes],
    seq2: Union[str, bytes],
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dna_conc: float = 50.0,
    temp_c: float = 37.0,
    max_loop: int = 30,
    output_structure: bool = False,
    dmso_conc: float = 0.0,
    dmso_fact: float = 0.6,
    formamide_conc: float = 0.0,
) -> ThermoResult: ...
def calc_hairpin_batch(
    seqs: List[Union[str, bytes]],
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dna_conc: float = 50.0,
    temp_c: float = 37.0,
    max_loop: int = 30,
    output_structure: bool = False,
    dmso_conc: float = 0.0,
    dmso_fact: float = 0.6,
    formamide_conc: float = 0.0,
) -> List[ThermoResult]:
    """calc_hairpin of each sequence, in parallel without the GIL."""
def calc_heterodimer_batch(
    pairs: List[Tuple[Union[str, bytes], Union[str, bytes]]],
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dna_conc: float = 50.0,
    temp_c: float = 37.0,
    max_loop: int = 30,
    output_structure: bool = False,
    dmso_conc: float = 0.0,
    dmso_fact: float = 0.6,
    formamide_conc: float = 0.0,
) -> List[ThermoResult]:
    """calc_heterodimer of each (seq1, seq2), in parallel without the GIL."""
def calc_tm(
    seq: Union[str, bytes],
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dna_conc: float = 50.0,
    dmso_conc: float = 0.0,
    dmso_fact: float = 0.6,
    formamide_conc: float = 0.0,
) -> float:
    """Tm (°C) of seq with its complement, on thal's model (the digest's).
    ValueError for an empty sequence or any base other than A/C/G/T."""
def calc_annealing_profile(
    seq: Union[str, bytes],
    temps: List[float],
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dna_conc: float = 50.0,
) -> List[float]:
    """% bound at each temperature. ValueError for an empty sequence or any
    base other than A/C/G/T."""
def calc_annealing_hetero(
    seq1: Union[str, bytes],
    seq2: Union[str, bytes],
    temp_c: float,
    mv_conc: float = 50.0,
    dv_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dna_conc: float = 50.0,
    max_loop: int = 30,
) -> float:
    """% of seq1 bound to seq2 (both 5'->3') via thal's heterodimer; 0 if no
    structure forms. Unlike calc_tm/calc_annealing, N and other non-ACGT
    bases are accepted (as unpairable, as primer3)."""
