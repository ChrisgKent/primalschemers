//! Python bindings for primer3-thal-rs: drop-in replacements for
//! primer3-py's `calc_hairpin`, `calc_homodimer` and `calc_heterodimer`.
//!
//! Signatures, defaults, results and errors follow primer3-py 2.3.0, which
//! primer3-thal-rs reproduces, so primalscheme3 can switch with only an
//! import change.

use primer3_thal_rs::{ThalArgs, ThalError};
use pyo3::exceptions::{PyRuntimeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyString};

/// Result of a thal calculation, mirroring primer3-py's `ThermoResult`.
#[pyclass(name = "ThermoResult", frozen, get_all)]
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
            if self.structure_found { "True" } else { "False" },
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

fn thal_args(
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    temp_c: f64,
    max_loop: i32,
) -> ThalArgs {
    ThalArgs {
        max_loop,
        mv: mv_conc,
        dv: dv_conc,
        dntp: dntp_conc,
        dna_conc,
        ..ThalArgs::default()
    }
    .with_temp_c(temp_c)
}

/// Hairpin thermodynamics of `seq`, as primer3-py's `calc_hairpin`.
///
/// Concentrations are mM (`mv_conc`, `dv_conc`, `dntp_conc`) and nM
/// (`dna_conc`); `temp_c` is in °C. Raises RuntimeError for an empty or
/// >60 nt sequence, ValueError for `max_loop` outside 0..=30.
#[pyfunction]
#[pyo3(signature = (seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30, output_structure=false))]
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
) -> PyResult<PyThermoResult> {
    let seq = seq_bytes(seq)?;
    let args = thal_args(mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_loop);
    let res = py
        .detach(|| primer3_thal_rs::calc_hairpin(&seq, &args, output_structure))
        .map_err(to_py_err)?;
    Ok(PyThermoResult::new(res, output_structure))
}

/// Homodimer thermodynamics of `seq`, as primer3-py's `calc_homodimer`.
/// Arguments and errors as `calc_hairpin`.
#[pyfunction]
#[pyo3(signature = (seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30, output_structure=false))]
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
) -> PyResult<PyThermoResult> {
    let seq = seq_bytes(seq)?;
    let args = thal_args(mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_loop);
    let res = py
        .detach(|| primer3_thal_rs::calc_homodimer(&seq, &args, output_structure))
        .map_err(to_py_err)?;
    Ok(PyThermoResult::new(res, output_structure))
}

/// Heterodimer thermodynamics of `seq1` binding `seq2` (both 5'->3'), as
/// primer3-py's `calc_heterodimer`. Arguments as `calc_hairpin`; only raises
/// for length if both sequences exceed 60 nt.
#[pyfunction]
#[pyo3(signature = (seq1, seq2, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30, output_structure=false))]
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
) -> PyResult<PyThermoResult> {
    let (seq1, seq2) = (seq_bytes(seq1)?, seq_bytes(seq2)?);
    let args = thal_args(mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_loop);
    let res = py
        .detach(|| primer3_thal_rs::calc_heterodimer(&seq1, &seq2, &args, output_structure))
        .map_err(to_py_err)?;
    Ok(PyThermoResult::new(res, output_structure))
}
