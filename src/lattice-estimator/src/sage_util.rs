use std::fmt;
use std::fmt::{Debug, Display};
use std::path::Path;
use std::process::Command;

// TODO: unfortunately, using pyo3 segfaults a lot when calling even simple python/sagemath scripts, use a dumb shell call for now
// use pyo3::{PyResult, Python};
// use pyo3::types::{IntoPyDict, PyList};
//
// pub(crate) fn sage_run<F, T>(f: F) -> PyResult<T>
//     where F: Fn(Python) -> PyResult<T> {
//     Python::with_gil(|py| {
//         let locals = [("sys", py.import("sys")?)].into_py_dict(py);
//         let syspath: &PyList = py.eval("sys.path", None, Some(&locals))?.extract()?;
//
//         // Add paths manually to ensure relevant python/sagemath files can be found
//         let root = Path::new(env!("CARGO_MANIFEST_DIR"));
//         // Ensure `import estimator' works
//         syspath.insert(0, root.join("lattice-estimator")).unwrap();
//         // Ensure `import sis' works
//         syspath.insert(0, root.join("src")).unwrap();
//
//         f(py)
//     })
// }

#[derive(Clone, Debug)]
pub struct SageMathError(String);

impl Display for SageMathError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SageMathError({})", self.0)
    }
}

// impl std::error::Error for SageMathError {
//     fn description(&self) -> &str {
//         &self.0
//     }
// }

impl<E: std::error::Error> From<E> for SageMathError {
    fn from(err: E) -> Self {
        SageMathError(err.to_string())
    }
}

pub(crate) fn sagemath_eval<F, T, E>(eval: String, parse: F) -> Result<T, SageMathError>
    where F: Fn(String) -> Result<T, E>,
          E: std::error::Error
{
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let output = Command::new("sage")
        .arg("--python")
        .arg("-c")
        .arg(
            format!("import sys;sys.path.insert(0, '{}');sys.path.insert(0, '{}');from sis import *;print({})",
                    root.join("lattice-estimator").to_str().ok_or(SageMathError("could not construct path".to_string()))?,
                    root.join("src").to_str().ok_or(SageMathError("could not construct path".to_string()))?,
                    eval)).output()?;
    if !output.status.success() {
        return Err(SageMathError(format!("Command {} terminated with exit code {}:\n {}", eval, output.status, String::from_utf8(output.stderr)?)));
    }
    let stdout = String::from_utf8(output.stdout).map_err(SageMathError::from)?;
    parse(stdout).map_err(SageMathError::from)
}