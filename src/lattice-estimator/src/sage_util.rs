use pyo3::{PyResult, Python};
use pyo3::types::{IntoPyDict, PyList};

pub(crate) fn sage_run<F, T>(f: F) -> PyResult<T>
    where F: Fn(Python) -> PyResult<T> {
    Python::with_gil(|py| {
        let locals = [("sys", py.import("sys")?)].into_py_dict(py);
        let syspath: &PyList = py.eval("sys.path", None, Some(&locals))?.extract()?;

        // Add paths manually to ensure relevant python/sagemath files can be found
        syspath.insert(0, "../lattice-estimator/lattice-estimator").unwrap();
        syspath.insert(0, "../lattice-estimator/src").unwrap();

        f(py)
    })
}