use pyo3::prelude::*;
use pyo3::types::{IntoPyDict, PyList};

enum Norm {
    L2,
    Linf,
}

fn sage_run<F, T>(f: F) -> PyResult<T>
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

fn sis_security_level(n: usize, q: u64, length_bound: f64, m: usize, norm: Norm) -> f64 {
    sage_run(|py| {
        let sis = py.import("sis")?;
        let security_level = match norm {
            Norm::L2 => sis.getattr("security_level_l2")?,
            Norm::Linf => sis.getattr("security_level_linf")?
        };
        security_level.call1((n, q, length_bound, m))?.extract()
    }).unwrap()
}

#[cfg(test)]
mod test {
    use super::*;

    // from lattice-estimator/schemes
    const FALCON512_UNF: (usize, u64, f64, usize, Norm) = (512, 12289, 5833.9072, 2426, Norm::L2);
    const DILITHIUM2_MSIS_WK_UNF: (usize, u64, f64, usize, Norm) = (1024, 8380417, 350209., 2304, Norm::Linf);


    // TODO: this segfaults on some runs, probably because something at the pyo3 <> sagemath boundary breaks

    #[test]
    fn test_sis_security_level_l2()
    {
        let (n, q, length_bound, m, norm) = FALCON512_UNF;
        let lambda = sis_security_level(n, q, length_bound, m, norm);
        println!("lambda: {:?}", lambda);
    }

    #[test]
    fn test_sis_security_level_linf()
    {
        let (n, q, length_bound, m, norm) = DILITHIUM2_MSIS_WK_UNF;
        let lambda = sis_security_level(n, q, length_bound, m, norm);
        println!("lambda: {:?}", lambda);
    }
}