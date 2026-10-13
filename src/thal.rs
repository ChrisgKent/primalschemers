//! Python bindings for primer3-thal-rs: drop-in replacements for
//! primer3-py's `calc_hairpin`, `calc_homodimer` and `calc_heterodimer`.
//!
//! Signatures, defaults, results and errors follow primer3-py 2.3.0, which
//! primer3-thal-rs reproduces, so primalscheme3 can switch with only an
//! import change. They also take optional DMSO/formamide (not in primer3's
//! thal; off by default), and `calc_tm` gives a primer's Tm on the same
//! model, as the digest uses.

use primer3_thal_rs::{ThalArgs, ThalError};
use pyo3::exceptions::{PyRuntimeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyString};

/// Result of a thal calculation, mirroring primer3-py's `ThermoResult`.
#[pyclass(name = "ThermoResult", frozen, get_all, skip_from_py_object)]
#[derive(Clone, Debug)]
pub struct PyThermoResult {
    /// Whether any structure was found.
    structure_found: bool,
    /// Melting temperature, °C.
    tm: f64,
    /// Gibbs free energy at `temp_c`, cal/mol.
    dg: f64,
    /// Enthalpy, cal/mol.
    dh: f64,
    /// Entropy, cal/K/mol.
    ds: f64,
    /// primer3's ASCII drawing: `None` unless `output_structure` was set,
    /// and `""` if it was set but no structure formed (as primer3-py).
    ascii_structure: Option<String>,
}

#[pymethods]
impl PyThermoResult {
    /// The drawing split into lines, or `None` if there is none (as primer3-py).
    #[getter]
    fn ascii_structure_lines(&self) -> Option<Vec<String>> {
        self.ascii_structure
            .as_deref()
            .filter(|s| !s.is_empty())
            .map(|s| s.trim_matches('\n').split('\n').map(String::from).collect())
    }

    /// The result as a dict, as primer3-py's `ThermoResult.todict()`.
    fn todict(&self, py: Python<'_>) -> PyResult<Py<pyo3::types::PyDict>> {
        let d = pyo3::types::PyDict::new(py);
        d.set_item("structure_found", self.structure_found)?;
        d.set_item("ascii_structure", self.ascii_structure.clone())?;
        d.set_item("tm", self.tm)?;
        d.set_item("dg", self.dg)?;
        d.set_item("dh", self.dh)?;
        d.set_item("ds", self.ds)?;
        Ok(d.unbind())
    }

    fn __repr__(&self) -> String {
        format!(
            "ThermoResult(structure_found={}, tm={:.2}, dg={:.2}, dh={:.2}, ds={:.2})",
            if self.structure_found {
                "True"
            } else {
                "False"
            },
            self.tm,
            self.dg,
            self.dh,
            self.ds
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

impl PyThermoResult {
    fn new(res: primer3_thal_rs::ThermoResult, output_structure: bool) -> Self {
        Self {
            structure_found: res.structure_found,
            tm: res.tm,
            dg: res.dg,
            dh: res.dh,
            ds: res.ds,
            // primer3-py reports "" when a drawing was requested but nothing formed.
            ascii_structure: output_structure.then(|| res.ascii_structure.unwrap_or_default()),
        }
    }
}

/// Errors as primer3-py raises them: bad `max_loop` is a ValueError, sequence
/// problems a RuntimeError.
fn to_py_err(e: ThalError) -> PyErr {
    match e {
        ThalError::InvalidMaxLoop(_) => PyValueError::new_err(e.to_string()),
        _ => PyRuntimeError::new_err(e.to_string()),
    }
}

/// A sequence given as `str` or `bytes`, as primer3-py accepts.
fn seq_bytes(seq: &Bound<'_, PyAny>) -> PyResult<Vec<u8>> {
    if let Ok(s) = seq.cast::<PyString>() {
        Ok(s.to_str()?.as_bytes().to_vec())
    } else if let Ok(b) = seq.cast::<PyBytes>() {
        Ok(b.as_bytes().to_vec())
    } else {
        Err(PyTypeError::new_err("sequence must be str or bytes"))
    }
}

#[allow(clippy::too_many_arguments)]
fn thal_args(
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    temp_c: f64,
    max_loop: i32,
    dmso_conc: f64,
    dmso_fact: f64,
    formamide_conc: f64,
) -> ThalArgs {
    ThalArgs {
        max_loop,
        mv: mv_conc,
        dv: dv_conc,
        dntp: dntp_conc,
        dna_conc,
        dmso_conc,
        dmso_fact,
        formamide_conc,
        ..ThalArgs::default()
    }
    .with_temp_c(temp_c)
}

/// Hairpin thermodynamics of `seq`, as primer3-py's `calc_hairpin`.
///
/// Concentrations are mM (`mv_conc`, `dv_conc`, `dntp_conc`) and nM
/// (`dna_conc`); `temp_c` is in °C. `dmso_conc` (%), `dmso_fact` (°C per %)
/// and `formamide_conc` (mol/l) shift Tm only, as primer3's oligotm; they
/// default to no correction. Raises RuntimeError for an empty or >60 nt
/// sequence, ValueError for `max_loop` outside 0..=30.
#[pyfunction]
#[pyo3(signature = (seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30, output_structure=false, dmso_conc=0.0, dmso_fact=0.6, formamide_conc=0.0))]
#[allow(clippy::too_many_arguments)]
pub fn calc_hairpin(
    py: Python<'_>,
    seq: &Bound<'_, PyAny>,
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    temp_c: f64,
    max_loop: i32,
    output_structure: bool,
    dmso_conc: f64,
    dmso_fact: f64,
    formamide_conc: f64,
) -> PyResult<PyThermoResult> {
    let seq = seq_bytes(seq)?;
    let args = thal_args(
        mv_conc,
        dv_conc,
        dntp_conc,
        dna_conc,
        temp_c,
        max_loop,
        dmso_conc,
        dmso_fact,
        formamide_conc,
    );
    let res = py
        .detach(|| primer3_thal_rs::calc_hairpin(&seq, &args, output_structure))
        .map_err(to_py_err)?;
    Ok(PyThermoResult::new(res, output_structure))
}

/// Homodimer thermodynamics of `seq`, as primer3-py's `calc_homodimer`.
/// Arguments and errors as `calc_hairpin`.
#[pyfunction]
#[pyo3(signature = (seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30, output_structure=false, dmso_conc=0.0, dmso_fact=0.6, formamide_conc=0.0))]
#[allow(clippy::too_many_arguments)]
pub fn calc_homodimer(
    py: Python<'_>,
    seq: &Bound<'_, PyAny>,
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    temp_c: f64,
    max_loop: i32,
    output_structure: bool,
    dmso_conc: f64,
    dmso_fact: f64,
    formamide_conc: f64,
) -> PyResult<PyThermoResult> {
    let seq = seq_bytes(seq)?;
    let args = thal_args(
        mv_conc,
        dv_conc,
        dntp_conc,
        dna_conc,
        temp_c,
        max_loop,
        dmso_conc,
        dmso_fact,
        formamide_conc,
    );
    let res = py
        .detach(|| primer3_thal_rs::calc_homodimer(&seq, &args, output_structure))
        .map_err(to_py_err)?;
    Ok(PyThermoResult::new(res, output_structure))
}

/// Heterodimer thermodynamics of `seq1` binding `seq2` (both 5'->3'), as
/// primer3-py's `calc_heterodimer`. Arguments as `calc_hairpin`; only raises
/// for length if both sequences exceed 60 nt.
#[pyfunction]
#[pyo3(signature = (seq1, seq2, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30, output_structure=false, dmso_conc=0.0, dmso_fact=0.6, formamide_conc=0.0))]
#[allow(clippy::too_many_arguments)]
pub fn calc_heterodimer(
    py: Python<'_>,
    seq1: &Bound<'_, PyAny>,
    seq2: &Bound<'_, PyAny>,
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    temp_c: f64,
    max_loop: i32,
    output_structure: bool,
    dmso_conc: f64,
    dmso_fact: f64,
    formamide_conc: f64,
) -> PyResult<PyThermoResult> {
    let (seq1, seq2) = (seq_bytes(seq1)?, seq_bytes(seq2)?);
    let args = thal_args(
        mv_conc,
        dv_conc,
        dntp_conc,
        dna_conc,
        temp_c,
        max_loop,
        dmso_conc,
        dmso_fact,
        formamide_conc,
    );
    let res = py
        .detach(|| primer3_thal_rs::calc_heterodimer(&seq1, &seq2, &args, output_structure))
        .map_err(to_py_err)?;
    Ok(PyThermoResult::new(res, output_structure))
}

/// `calc_hairpin` of each of `seqs`, computed in parallel without the GIL.
/// Arguments as `calc_hairpin`; raises on the first sequence that fails.
#[pyfunction]
#[pyo3(signature = (seqs, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30, output_structure=false, dmso_conc=0.0, dmso_fact=0.6, formamide_conc=0.0))]
#[allow(clippy::too_many_arguments)]
pub fn calc_hairpin_batch(
    py: Python<'_>,
    seqs: Vec<Bound<'_, PyAny>>,
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    temp_c: f64,
    max_loop: i32,
    output_structure: bool,
    dmso_conc: f64,
    dmso_fact: f64,
    formamide_conc: f64,
) -> PyResult<Vec<PyThermoResult>> {
    let seqs = seqs.iter().map(seq_bytes).collect::<PyResult<Vec<_>>>()?;
    let args = thal_args(
        mv_conc,
        dv_conc,
        dntp_conc,
        dna_conc,
        temp_c,
        max_loop,
        dmso_conc,
        dmso_fact,
        formamide_conc,
    );
    py.detach(|| primer3_thal_rs::calc_hairpin_batch(&seqs, &args, output_structure))
        .into_iter()
        .map(|r| r.map(|res| PyThermoResult::new(res, output_structure)).map_err(to_py_err))
        .collect()
}

/// `calc_heterodimer` of each `(seq1, seq2)` in `pairs`, computed in
/// parallel without the GIL. Arguments as `calc_heterodimer`; raises on the
/// first pair that fails.
#[pyfunction]
#[pyo3(signature = (pairs, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30, output_structure=false, dmso_conc=0.0, dmso_fact=0.6, formamide_conc=0.0))]
#[allow(clippy::too_many_arguments)]
pub fn calc_heterodimer_batch(
    py: Python<'_>,
    pairs: Vec<(Bound<'_, PyAny>, Bound<'_, PyAny>)>,
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    temp_c: f64,
    max_loop: i32,
    output_structure: bool,
    dmso_conc: f64,
    dmso_fact: f64,
    formamide_conc: f64,
) -> PyResult<Vec<PyThermoResult>> {
    let pairs = pairs
        .iter()
        .map(|(a, b)| Ok((seq_bytes(a)?, seq_bytes(b)?)))
        .collect::<PyResult<Vec<_>>>()?;
    let args = thal_args(
        mv_conc,
        dv_conc,
        dntp_conc,
        dna_conc,
        temp_c,
        max_loop,
        dmso_conc,
        dmso_fact,
        formamide_conc,
    );
    py.detach(|| primer3_thal_rs::calc_heterodimer_batch(&pairs, &args, output_structure))
        .into_iter()
        .map(|r| r.map(|res| PyThermoResult::new(res, output_structure)).map_err(to_py_err))
        .collect()
}

/// Tm (°C) of `seq` paired with its complement: equal to
/// `calc_heterodimer(seq, reverse_complement(seq)).tm` at the same
/// conditions, and the Tm the digest checks (thal's nearest-neighbour
/// model, not primer3's `calc_tm`). Arguments as `calc_hairpin`. Raises
/// ValueError for an empty sequence or a base other than A/C/G/T.
#[pyfunction]
#[pyo3(signature = (seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, dmso_conc=0.0, dmso_fact=0.6, formamide_conc=0.0))]
#[allow(clippy::too_many_arguments)]
pub fn calc_tm(
    seq: &Bound<'_, PyAny>,
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    dmso_conc: f64,
    dmso_fact: f64,
    formamide_conc: f64,
) -> PyResult<f64> {
    let seq = seq_bytes(seq)?;
    let args = ThalArgs {
        mv: mv_conc,
        dv: dv_conc,
        dntp: dntp_conc,
        dna_conc,
        dmso_conc,
        dmso_fact,
        formamide_conc,
        ..ThalArgs::default()
    };
    crate::tm::oligo_tm(&seq, &args).map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Percentage of `seq` bound to its complement at each of `temps` (°C), in
/// order: `calc_annealing` at each temperature, on the digest's model.
/// Raises ValueError for an empty sequence or a base other than A/C/G/T.
#[pyfunction]
#[pyo3(signature = (seq, temps, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0))]
pub fn calc_annealing_profile(
    seq: &Bound<'_, PyAny>,
    temps: Vec<f64>,
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
) -> PyResult<Vec<f64>> {
    let seq = seq_bytes(seq)?;
    let args = ThalArgs {
        mv: mv_conc,
        dv: dv_conc,
        dntp: dntp_conc,
        dna_conc,
        ..ThalArgs::default()
    };
    crate::tm::oligo_annealing_profile(&seq, &args, &temps)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Percentage of `seq1` bound to `seq2` (both 5'->3', mismatches allowed) at
/// `temp_c` °C, from thal's heterodimer: equal to `calc_annealing` when
/// `seq2` is `seq1`'s reverse complement. For a primer on a template site,
/// pass the site's reverse complement as `seq2`. 0 if no structure forms.
/// Errors as `calc_heterodimer`.
#[pyfunction]
#[pyo3(signature = (seq1, seq2, temp_c, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, max_loop=30))]
#[allow(clippy::too_many_arguments)]
pub fn calc_annealing_hetero(
    py: Python<'_>,
    seq1: &Bound<'_, PyAny>,
    seq2: &Bound<'_, PyAny>,
    temp_c: f64,
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    max_loop: i32,
) -> PyResult<f64> {
    let (seq1, seq2) = (seq_bytes(seq1)?, seq_bytes(seq2)?);
    let args = thal_args(
        mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_loop, 0.0, 0.6, 0.0,
    );
    py.detach(|| crate::tm::hetero_annealing(&seq1, &seq2, &args, temp_c))
        .map_err(to_py_err)
}
