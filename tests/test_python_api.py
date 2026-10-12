"""Python API tests for the primer3-thal-rs integration.

Run after building the extension (e.g. `maturin develop --release`):

    pip install pytest "primer3-py==2.3.0"
    pytest tests/

Tests that compare against primer3-py are skipped if it isn't installed.
"""

import random
import re

import pytest

import primalschemers as ps

# primalscheme3's conditions (matching digest.rs) and primer3-py's defaults.
PS3_CONDITIONS = dict(mv_conc=100, dv_conc=2, dntp_conc=0.8, dna_conc=15)
CONDITIONS = [PS3_CONDITIONS, {}]
FIELDS = ("tm", "dg", "dh", "ds")


def rseq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def rc(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def assert_same(ours, ref, label):
    assert ours.structure_found == ref.structure_found, label
    for f in FIELDS:
        assert getattr(ours, f) == pytest.approx(getattr(ref, f), abs=1e-6), (label, f)
    assert ours.ascii_structure == ref.ascii_structure, label
    assert ours.ascii_structure_lines == ref.ascii_structure_lines, label


# --- calc_* against primer3-py ---------------------------------------------


def primer_cases(seed, n=150):
    """Random primers, designed hairpins and primer pairs with 3' overlaps."""
    rng = random.Random(seed)
    cases = []
    for _ in range(n):
        a = rseq(rng, rng.randint(18, 34))
        stem = rseq(rng, rng.randint(5, 9))
        hairpin = rseq(rng, rng.randint(0, 8)) + stem + rseq(rng, rng.randint(3, 8)) + rc(stem)
        k = rng.randint(4, 10)
        partner = rseq(rng, rng.randint(10, 24)) + rc(a[-k:]) if rng.random() < 0.4 else rseq(rng, rng.randint(18, 34))
        cases.append((a, hairpin, partner))
    return cases


@pytest.mark.parametrize("conditions", CONDITIONS, ids=["ps3", "defaults"])
@pytest.mark.parametrize("output_structure", [False, True])
def test_calc_functions_match_primer3(conditions, output_structure):
    primer3 = pytest.importorskip("primer3")
    for a, hairpin, partner in primer_cases(seed=11):
        kw = dict(conditions, output_structure=output_structure)
        for seq in (a, hairpin):
            assert_same(ps.calc_hairpin(seq, **kw), primer3.calc_hairpin(seq, **kw), ("hairpin", seq))
        assert_same(ps.calc_homodimer(a, **kw), primer3.calc_homodimer(a, **kw), ("homodimer", a))
        assert_same(
            ps.calc_heterodimer(a, partner, **kw),
            primer3.calc_heterodimer(a, partner, **kw),
            ("heterodimer", a, partner),
        )


def test_calc_conditions_change_results():
    seq = "CCCCCATCCGATCAGGGGG"
    default = ps.calc_hairpin(seq)
    salty = ps.calc_hairpin(seq, mv_conc=500, dv_conc=10)
    hot = ps.calc_hairpin(seq, temp_c=70)
    assert salty.tm != default.tm
    assert hot.dg != default.dg and hot.tm == pytest.approx(default.tm)
    assert ps.calc_hairpin(seq, max_loop=0).structure_found


# --- ThermoResult behaviour --------------------------------------------------


def test_thermo_result_attributes():
    res = ps.calc_hairpin("CCCCCATCCGATCAGGGGG", output_structure=True)
    assert isinstance(res, ps.ThermoResult)
    assert res.structure_found
    assert res.ascii_structure == "SEQ\t/////---------\\\\\\\\\\\nSTR\tCCCCCATCCGATCAGGGGG\n"
    assert res.ascii_structure_lines == ["SEQ\t/////---------\\\\\\\\\\", "STR\tCCCCCATCCGATCAGGGGG"]
    assert res.todict() == {
        "structure_found": True,
        "ascii_structure": res.ascii_structure,
        "tm": res.tm,
        "dg": res.dg,
        "dh": res.dh,
        "ds": res.ds,
    }
    assert re.fullmatch(
        r"ThermoResult\(structure_found=True, tm=-?\d+\.\d\d, dg=-?\d+\.\d\d, dh=-?\d+\.\d\d, ds=-?\d+\.\d\d\)",
        repr(res),
    )
    assert str(res) == repr(res)
    with pytest.raises(AttributeError):
        res.tm = 1.0  # frozen


def test_ascii_structure_none_vs_empty():
    # As primer3-py: None unless requested; "" if requested but nothing forms.
    assert ps.calc_hairpin("ACGT").ascii_structure is None
    none = ps.calc_hairpin("ACGT", output_structure=True)
    assert not none.structure_found
    assert none.ascii_structure == "" and none.ascii_structure_lines is None
    assert none.tm == none.dg == none.dh == none.ds == 0.0
    assert ps.calc_hairpin("CCCCCATCCGATCAGGGGG").ascii_structure_lines is None


def test_repr_matches_primer3():
    primer3 = pytest.importorskip("primer3")
    for seq in ("CCCCCATCCGATCAGGGGG", "ACGT"):
        assert repr(ps.calc_hairpin(seq)) == repr(primer3.calc_hairpin(seq))


def test_str_and_bytes_input():
    a, b = "GTAAAACGACGGCCAGT", "ACTGGCCGTCGTTTTAC"
    assert ps.calc_heterodimer(a, b).todict() == ps.calc_heterodimer(a.encode(), b.encode()).todict()
    assert ps.calc_heterodimer(a, b.encode()).todict() == ps.calc_heterodimer(a.encode(), b).todict()
    assert ps.calc_hairpin(a).todict() == ps.calc_hairpin(a.encode()).todict()
    assert ps.calc_homodimer(a).todict() == ps.calc_homodimer(a.encode()).todict()
    # Case-insensitive, as primer3
    assert ps.calc_hairpin(a.lower()).todict() == ps.calc_hairpin(a).todict()


# --- errors, as primer3-py raises them -------------------------------------


@pytest.mark.parametrize(
    "call, exc",
    [
        (lambda: ps.calc_hairpin(""), RuntimeError),
        (lambda: ps.calc_hairpin("A" * 61), RuntimeError),
        (lambda: ps.calc_homodimer("A" * 61), RuntimeError),
        (lambda: ps.calc_heterodimer("A" * 61, "A" * 61), RuntimeError),
        (lambda: ps.calc_heterodimer("", "ACGT"), RuntimeError),
        (lambda: ps.calc_hairpin("ACGTACGT", max_loop=31), ValueError),
        (lambda: ps.calc_hairpin("ACGTACGT", max_loop=-1), ValueError),
        (lambda: ps.calc_hairpin(123), TypeError),
        (lambda: ps.calc_heterodimer("ACGT", None), TypeError),
    ],
)
def test_errors(call, exc):
    with pytest.raises(exc):
        call()


def test_errors_match_primer3():
    primer3 = pytest.importorskip("primer3")
    for args, kw in [(("",), {}), (("A" * 61,), {}), (("ACGTACGT",), {"max_loop": 31})]:
        with pytest.raises(Exception) as ours:
            ps.calc_hairpin(*args, **kw)
        with pytest.raises(Exception) as ref:
            primer3.calc_hairpin(*args, **kw)
        assert type(ours.value) is type(ref.value)


# --- one thermodynamic model ---------------------------------------------


# The digest's conditions (digest.rs `conditions()`)
DIGEST_CONDITIONS = dict(PS3_CONDITIONS, formamide_conc=0.8)
MODEL_CONDITIONS = [
    {},
    PS3_CONDITIONS,
    DIGEST_CONDITIONS,
    dict(dmso_conc=5.0, dmso_fact=0.6, formamide_conc=1.5),
]


def gc_fraction(*seqs):
    joined = "".join(seqs).upper()
    return sum(c in "GC" for c in joined) / len(joined)


def tm_shift(gc, dmso_conc=0.0, dmso_fact=0.6, formamide_conc=0.0, **_):
    return -dmso_conc * dmso_fact + (0.453 * gc - 2.88) * formamide_conc


@pytest.mark.parametrize("conditions", MODEL_CONDITIONS)
def test_calc_tm_equals_thal_duplex_tm(conditions):
    rng = random.Random(21)
    seqs = [rseq(rng, rng.randint(15, 34)) for _ in range(200)]
    seqs += ["ACGTACGTACGTACGTACGT", "GAATTCGAATTC", "AAAAAAAAAAAAAAAAAAAA"]
    for s in seqs:
        thal = ps.calc_heterodimer(s, rc(s), **conditions).tm
        assert ps.calc_tm(s, **conditions) == pytest.approx(thal, abs=1e-9), s


def test_calc_tm_input():
    s = "ACCAACGATGGTGTGTCCAT"
    assert ps.calc_tm(s) == ps.calc_tm(s.encode()) == ps.calc_tm(s.lower())
    for bad in ["", "ACGTNACGT", "ACGR"]:
        with pytest.raises(ValueError):
            ps.calc_tm(bad)
    with pytest.raises(TypeError):
        ps.calc_tm(123)


def test_dmso_and_formamide_shift_only_tm():
    hp, a, b = "CCCCCATCCGATCAGGGGG", "GTAAAACGACGGCCAGTGAAT", "ACACCTGTGCCTGTTAAACCAT"
    extra = dict(dmso_conc=2.0, formamide_conc=0.8)
    cases = [
        (ps.calc_hairpin, (hp,), gc_fraction(hp)),
        (ps.calc_homodimer, (a,), gc_fraction(a)),
        (ps.calc_heterodimer, (a, b), gc_fraction(a, b)),
    ]
    for f, seqs, gc in cases:
        plain = f(*seqs, **PS3_CONDITIONS, output_structure=True)
        corrected = f(*seqs, **PS3_CONDITIONS, **extra, output_structure=True)
        assert plain.structure_found
        assert corrected.tm == pytest.approx(plain.tm + tm_shift(gc, **extra), abs=1e-9)
        assert (corrected.dg, corrected.dh, corrected.ds) == (plain.dg, plain.dh, plain.ds)
        assert corrected.ascii_structure == plain.ascii_structure


def test_calc_annealing_is_half_bound_at_calc_tm():
    rng = random.Random(4)
    for _ in range(50):
        s = rseq(rng, rng.randint(18, 30))
        tm = ps.calc_tm(s, **PS3_CONDITIONS)
        assert ps.calc_annealing(s, 15.0, 100.0, 2.0, 0.8, 0.0, 0.6, 0.0, tm) == pytest.approx(50.0, abs=1e-6)
    with pytest.raises(ValueError):
        ps.calc_annealing("ACGTN", 15.0, 100.0, 2.0, 0.8, 0.0, 0.6, 0.0, 60.0)


def test_tm_mode_digest_uses_calc_tm(msa_path):
    f, r, _ = ps.Digester(msa_path, 1, False).digest(
        thermo_check=True, primer_tm_min=59.5, primer_tm_max=62.5
    )
    assert f + r
    for k in f + r:
        for s in k.seqs():
            assert 59.5 <= ps.calc_tm(s, **DIGEST_CONDITIONS) < 62.5, s


def test_one_long_heterodimer_sequence_allowed():
    # As primer3: only both > 60 nt is an error
    long_target = "ACGTTGCA" * 20
    assert ps.calc_heterodimer("ACGTACGTACGTACGT", long_target).structure_found in (True, False)
    assert ps.calc_heterodimer(long_target, "ACGTACGTACGTACGT") is not None


# --- digest: hairpin threshold and return_all -------------------------------


@pytest.fixture(scope="module")
def msa_path(tmp_path_factory):
    """A 3 kb alignment of three variants with inverted repeats, so some
    kmers form hairpins, plus an N run and an ambiguity code."""
    rng = random.Random(3)
    parts = []
    for _ in range(60):
        parts.append(rseq(rng, rng.randint(20, 40)))
        stem = rseq(rng, rng.randint(6, 9))
        parts.append(stem + rseq(rng, rng.randint(3, 6)) + rc(stem))
    genome = "".join(parts)
    seqs = [genome]
    for _ in range(4):
        g = list(genome)
        for i in rng.sample(range(len(g)), 15):
            g[i] = rng.choice("ACGT")
        seqs.append("".join(g))
    seqs[2] = seqs[2][:1500] + "NNNNNN" + seqs[2][1506:]
    seqs[3] = seqs[3][:2000] + "R" + seqs[3][2001:]
    path = tmp_path_factory.mktemp("msa") / "msa.fasta"
    path.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(seqs)))
    return str(path)


def kmer_key(k):
    return ("F", k.end) if isinstance(k, ps.FKmer) else ("R", k.start)


def by_key(kmers):
    return {kmer_key(k): (k.seqs(), k.counts()) for k in kmers}


def test_hairpin_threshold_matches_python_check(msa_path):
    primer3 = pytest.importorskip("primer3")
    d = ps.Digester(msa_path, 1, False)

    def hairpin_tm(s):
        # primer3's hairpin Tm, shifted by the digest's formamide when a
        # structure forms
        res = primer3.calc_hairpin(s, **PS3_CONDITIONS)
        shift = tm_shift(gc_fraction(s), **DIGEST_CONDITIONS) if res.structure_found else 0.0
        return res.tm + shift

    def forms_hairpin(seqs, th):  # primalscheme3's core/thermo.py check
        return any(hairpin_tm(s) > th for s in seqs)

    f, r, _ = d.digest(thermo_check=True, primer_hairpin_th_max=1000.0)
    loose = f + r
    for th in (55.0, 45.0):
        f, r, _ = d.digest(thermo_check=True, primer_hairpin_th_max=th)
        expected = {kmer_key(k) for k in loose if not forms_hairpin(k.seqs(), th)}
        assert {kmer_key(k) for k in f + r} == expected
        assert len(expected) < len(loose)
    # The default threshold is 55 °C
    f, r, _ = d.digest(thermo_check=True)
    f55, r55, _ = d.digest(thermo_check=True, primer_hairpin_th_max=55.0)
    assert by_key(f + r) == by_key(f55 + r55)


def annealing(seq, temp_c):
    # The digest's conditions (digest.rs `conditions()`)
    return ps.calc_annealing(seq, 15.0, 100.0, 2.0, 0.8, 0.0, 0.6, 0.8, temp_c)


@pytest.mark.parametrize("entry", ["digest", "digest_seq"])
@pytest.mark.parametrize("temp_c,prop,diff", [(60.0, 50.0, 5.0), (65.0, 30.0, 10.0)])
def test_annealing_mode_filters_on_annealing(msa_path, entry, temp_c, prop, diff):
    def run(**kw):
        kw = dict(thermo_check=True, primer_annealing_prop=prop, **kw)
        if entry == "digest":
            f, r, _ = ps.Digester(msa_path, 1, False).digest(**kw)
        else:
            f, r, _ = ps.digest_seq(msa_path, 1, False, **kw)
        return f + r

    kmers = run(primer_annealing_diff=diff, annealing_temp_c=temp_c)
    assert kmers
    for k in kmers:
        for s in k.seqs():
            a = annealing(s.decode() if isinstance(s, bytes) else s, temp_c)
            assert prop - diff <= a <= prop + diff, (s, a)
    # The filter does work: an open window lets some outside it through
    wide = run(primer_annealing_diff=100.0, annealing_temp_c=temp_c)
    assert any(
        not prop - diff <= annealing(s.decode() if isinstance(s, bytes) else s, temp_c) <= prop + diff
        for k in wide
        for s in k.seqs()
    )
    # and the temperature matters
    other = run(primer_annealing_diff=diff, annealing_temp_c=temp_c - 10.0)
    assert {kmer_key(k) for k in other} != {kmer_key(k) for k in kmers}


def test_hairpin_not_applied_without_thermo_check(msa_path):
    d = ps.Digester(msa_path, 1, False)
    f, r, _ = d.digest(thermo_check=False, primer_hairpin_th_max=0.0)
    f2, r2, _ = d.digest(thermo_check=False)
    assert by_key(f + r) == by_key(f2 + r2) and len(f + r) > 0


@pytest.mark.parametrize("min_freq", [0.0, 0.3])
@pytest.mark.parametrize("ignore_n", [False, True])
def test_return_all_matches_normal_digest(msa_path, min_freq, ignore_n):
    d = ps.Digester(msa_path, 1, False)
    kw = dict(thermo_check=True, min_freq=min_freq, ignore_n=ignore_n)
    f0, r0, _ = d.digest(**kw)
    f1, r1, _ = d.digest(return_all=True, **kw)
    normal = by_key(f0 + r0)
    assert all(k.statuses() is None and k.passed() and k.dimers() == [] for k in f0 + r0)

    returned = f1 + r1
    assert len(returned) > len(normal)
    multi = 0
    for k in returned:
        statuses = k.statuses()
        assert statuses is not None and len(statuses) == len(k.seqs()) == len(k.counts())
        multi += any(len(t) > 1 for t in statuses)
        if k.passed():
            keep = [i for i, t in enumerate(statuses) if "LowFreq" not in t]
            assert normal[kmer_key(k)] == ([k.seqs()[i] for i in keep], [k.counts()[i] for i in keep])
        else:
            assert kmer_key(k) not in normal
    assert sum(k.passed() for k in returned) == len(normal)
    assert multi > 0
    if min_freq > 0:
        assert any("LowFreq" in t for k in returned for t in k.statuses())


def test_return_all_status_vocabulary(msa_path):
    d = ps.Digester(msa_path, 1, False)
    f, r, _ = d.digest(thermo_check=True, return_all=True)
    seq_tags = {t for k in f + r for tags in k.statuses() for t in tags}
    kmer_tags = {t for k in f + r for t in k.kmer_statuses()}
    error_tags = {t for k in f + r for tags, _ in k.errors() for t in tags}
    known = {"LowGC", "HighGC", "Homopolymer", "HighTm", "LowTm", "Hairpin", "Fail", "LowFreq"}
    assert seq_tags and seq_tags <= known
    assert "Hairpin" in seq_tags
    assert kmer_tags <= {"PrimerDimer"}
    assert error_tags <= {
        "InvalidBase", "WalkedOutLeft", "WalkedOutRight", "GapOnSetBase", "NoValidPrimer",
        "MaxWalk", "ToLong", "ContainsN", "LowFreq", "Ignored",
    }
    assert all(count > 0 for k in f + r for _, count in k.errors())


def n_errors(kmers):
    return [tags for k in kmers for tags, _ in k.errors() if "ContainsN" in tags]


def test_n_reporting_with_and_without_ignore_n(msa_path):
    d = ps.Digester(msa_path, 1, False)
    # Without ignore_n: ContainsN fails the kmer
    f, r, _ = d.digest(thermo_check=True, return_all=True, ignore_n=False)
    strict = n_errors(f + r)
    assert strict and all("Ignored" not in t for t in strict)
    assert all(not k.passed() for k in f + r if any("ContainsN" in t for t, _ in k.errors()))
    # With ignore_n: still reported, tagged Ignored, and doesn't fail the kmer
    f, r, _ = d.digest(thermo_check=True, return_all=True, ignore_n=True)
    ignored = n_errors(f + r)
    assert ignored and all("Ignored" in t for t in ignored)
    assert any(k.passed() for k in f + r if any("ContainsN" in t for t, _ in k.errors()))
    # The same N positions are reported either way
    assert len(ignored) == len(strict)


def test_lowercase_n_in_fasta(tmp_path, msa_path):
    # The same alignment with its N run written in lower case
    text = open(msa_path).read()
    assert "NNNNNN" in text
    lower = tmp_path / "lower.fasta"
    lower.write_text(text.replace("NNNNNN", "nnnnnn"))
    for ignore_n in (False, True):
        kw = dict(thermo_check=True, return_all=True, ignore_n=ignore_n)
        f1, r1, logs1 = ps.Digester(msa_path, 1, False).digest(**kw)
        f2, r2, logs2 = ps.Digester(str(lower), 1, False).digest(**kw)
        assert by_key(f1 + r1) == by_key(f2 + r2)
        assert [k.errors() for k in f1 + r1] == [k.errors() for k in f2 + r2]
        assert logs1 == logs2
        assert n_errors(f2 + r2)
        assert "InvalidBase" not in str(logs2)


def test_return_all_dimers(msa_path):
    d = ps.Digester(msa_path, 1, False)
    f, r, _ = d.digest(thermo_check=True, return_all=True)
    with_dimers = 0
    for k in f + r:
        pairs = k.dimers()
        assert bool(pairs) == ("PrimerDimer" in k.kmer_statuses())
        assert pairs == sorted(set(pairs))
        for i, j in pairs:
            assert 0 <= i <= j < len(k.seqs())
            assert "LowFreq" not in k.statuses()[i] and "LowFreq" not in k.statuses()[j]
            assert ps.do_seqs_interact(k.seqs()[i].encode(), k.seqs()[j].encode(), -26.0)
        with_dimers += bool(pairs)
    assert with_dimers > 0


def test_return_all_log_counts(msa_path):
    d = ps.Digester(msa_path, 1, False)
    f, r, logs = d.digest(thermo_check=True, return_all=True)
    fail_count = sum(not k.passed() for k in f)
    pass_count = sum(k.passed() for k in f)
    assert f"(ThermoResult(Fail), {fail_count})" in logs[0]
    assert f"(Pass, {pass_count})" in logs[0]


def test_return_all_with_remap_and_digest_seq(msa_path):
    d = ps.Digester(msa_path, 1, True)
    f, r, _ = d.digest(thermo_check=True, return_all=True)
    assert all(k.statuses() is not None for k in f + r)
    f2, r2, _ = ps.digest_seq(msa_path, 1, True, return_all=True, primer_hairpin_th_max=55.0)
    assert by_key(f + r) == by_key(f2 + r2)
    assert [k.statuses() for k in f] == [k.statuses() for k in f2]
