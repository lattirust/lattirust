use std::fmt;
use std::fmt::{Debug, Display};
use crate::norms::Norm;
use crate::errors::LatticeEstimatorError;
use crate::sage_util::sage_run;

pub struct SIS {
    n: usize,
    q: u64,
    length_bound: f64,
    m: usize,
    norm: Norm,
}

impl Display for SIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "SIS_{{n={}, q={}, length_bound={}, m={}, norm={}}}", self.n, self.q, self.length_bound, self.m, self.norm)
    }
}

impl Debug for SIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "SIS_{{n={}, q={}, length_bound={}, m={}, norm={}}}", self.n, self.q, self.length_bound, self.m, self.norm)
    }
}

impl SIS {
    pub const fn new(n: usize, q: u64, length_bound: f64, m: usize, norm: Norm) -> Self {
        SIS { n, q, length_bound, m, norm }
    }

    pub const fn with_n(&self, n: usize) -> Self {
        SIS { n, q: self.q, length_bound: self.length_bound, m: self.m, norm: self.norm }
    }

    /// Return lambda such that SIS_{n, q, length_bound, m} is 2^lambda-hard (for a given norm).
    /// Internally, this calls out to the lattice-estimator via a wrapper Python script.
    pub fn security_level(&self) -> f64 {
        sage_run(|py| {
            let sis = py.import("sis")?;
            let security_level = match self.norm {
                Norm::L2 => sis.getattr("security_level_l2")?,
                Norm::Linf => sis.getattr("security_level_linf")?
            };
            security_level.call1((self.n, self.q, self.length_bound, self.m))?.extract()
        }).unwrap()
    }

    /// Return the smallest m such that SIS_{n, q, length_bound, m} is 2^lambda-hard (for a given norm).
    pub fn find_optimal_n(&self, lambda: usize) -> Result<usize, LatticeEstimatorError> {
        let mut hi: usize = self.m; // (m as f64 / (q as f64).log2()).floor() as usize;
        let mut lo: usize = 1;

        let sis = self.with_n(hi);
        let lambda_hi = sis.security_level();
        debug_assert!(lambda_hi >= lambda as f64, "{sis} has sec. param. {lambda_hi}  < target lambda = {lambda}");
        // Loop invariant: SIS_{hi, q, length_bound, m} is 2^lambda_hi-hard with lambda_hi >= lambda
        while hi > lo {
            let mid = lo + (hi - lo) / 2;
            let sis = self.with_n(mid);
            let mid_lambda = sis.security_level();
            if mid_lambda >= lambda as f64 { // Search for smaller n in [lo, mid]
                hi = mid;
            } else { // Search for smaller n in [mid+1, hi]
                lo = mid + 1;
            }
        }
        assert_eq!(hi, lo);
        Ok(hi)
    }
}


#[cfg(test)]
mod test {
    use crate::norms::Norm;
    use crate::sis::SIS;

    // from lattice-estimator/schemes
    const FALCON512_UNF: SIS = SIS::new(512, 12289, 5833.9072, 1024, Norm::L2);
    const DILITHIUM2_MSIS_WK_UNF: SIS = SIS::new(1024, 8380417, 350209., 2304, Norm::Linf);


    // TODO: this segfaults on some runs, probably because something at the pyo3 <> sagemath boundary breaks

    #[test]
    fn test_sis_security_level_l2()
    {
        let lambda = FALCON512_UNF.security_level();
        println!("{FALCON512_UNF} -> lambda: {lambda}");
    }

    #[test]
    fn test_sis_security_level_linf()
    {
        let lambda = DILITHIUM2_MSIS_WK_UNF.security_level();
        println!("{DILITHIUM2_MSIS_WK_UNF} -> lambda: {lambda}");
    }

    #[test]
    fn test_find_optimal_n_l2()
    {
        let n_opt = FALCON512_UNF.find_optimal_n(128).unwrap();
        let sis = FALCON512_UNF.with_n(n_opt);
        let lambda = sis.security_level();
        assert!(lambda >= 128.0);
        println!("{FALCON512_UNF} -> lambda: {lambda}");
        println!("{sis} -> lambda: {lambda}");
    }

    #[test]
    fn test_find_optimal_n_linf()
    {
        let n_opt = DILITHIUM2_MSIS_WK_UNF.find_optimal_n(128).unwrap();
        let sis = DILITHIUM2_MSIS_WK_UNF.with_n(n_opt);
        let lambda = sis.security_level();
        assert!(lambda >= 128.0);
        println!("{DILITHIUM2_MSIS_WK_UNF} -> lambda: {lambda}");
        println!("{sis} -> lambda: {lambda}");
    }
}
