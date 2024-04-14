use std::fmt;
use std::fmt::{Debug, Display};
use log::debug;

use crate::errors::LatticeEstimatorError;
use crate::norms::Norm;
use crate::sage_util::sagemath_eval;
use crate::sis::SIS;

/// An MSIS instance is a matrix A in R_q^{m x n}, R_q = Z_q[X]/(X^d+1) such that A * s = 0 for some s in R_q^n with ||s||_norm <= length_bound.
pub struct MSIS {
    pub n: usize,
    pub d: usize,
    pub q: u64,
    pub length_bound: f64,
    pub m: usize,
    pub norm: Norm,
}

impl Display for MSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "MSIS_{{n={}, d={}, q={}, length_bound={}, m={}, norm={}}}", self.n, self.d, self.q, self.length_bound, self.m, self.norm)
    }
}

impl Debug for MSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "MSIS_{{n={}, d={}, q={}, length_bound={}, m={}, norm={}}}", self.n, self.d, self.q, self.length_bound, self.m, self.norm)
    }
}

impl MSIS {
    pub fn with_n(&self, n: usize) -> Self {
        MSIS { n, d: self.d, q: self.q, length_bound: self.length_bound, m: self.m, norm: self.norm }
    }

    pub const fn with_length_bound(&self, length_bound: f64) -> Self {
        MSIS { n: self.n, d: self.d, q: self.q, length_bound, m: self.m, norm: self.norm }
    }

    pub fn to_sis(&self) -> SIS {
        SIS::new(self.n * self.d, self.q, self.length_bound, self.m * self.d, self.norm)
    }

    pub fn upper_bound_n(&self) -> usize {
        self.to_sis().upper_bound_n().div_floor(self.d)
    }

    /// Return lambda such that MSIS_{n, d, q, length_bound, m} is 2^lambda-hard (for a given norm).
    /// We estimate the security by reducing to SIS_{n*d, q, length_bound, m*d} and calling the SIS security estimator.
    pub fn security_level(&self) -> f64 {
        let func = match self.norm {
            Norm::L2 => "msis_security_level_l2",
            Norm::Linf => "msis_security_level_linf"
        };
        sagemath_eval(format!("{}(n={}, d={}, q={}, length_bound={}, m={})", func, self.n, self.d, self.q, self.length_bound, self.m), SIS::parse_f64).unwrap()
    }

    /// Return the smallest n such that MSIS_{n, d, q, length_bound, m} is 2^lambda-hard (for a given norm).
    pub fn find_optimal_n(&self, lambda: usize) -> Result<usize, LatticeEstimatorError> {
        let lo: usize = 1;
        let hi: usize = self.upper_bound_n();

        debug!("Searching optimal n in [{lo}, {hi}] for {self} with target lambda={lambda}...");
        // Start from n=1 and exhaustively test upwards, since the upper bound hi is often not very large, and computing the security level for larger n is expensive.
        for n in lo..=hi {
            let curr = self.with_n(n);
            let lambda_curr = curr.security_level();
            debug!("\t{curr} -> {lambda_curr} bits of security");
            if lambda_curr >= lambda as f64 {
                return Ok(n);
            }
        }
        Err(LatticeEstimatorError::from(format!("no suitable n found in range [{lo}, {hi}] for {self}")))
    }

    /// Return the smallest n such that MSIS_{n, d, q, length_bound(n), m} is 2^lambda-hard (for a given norm), where length_bound is a function of n.
    pub fn find_optimal_n_dynamic<F>(&self, length_bound: F, lambda: usize) -> Result<usize, LatticeEstimatorError>
        where F: Fn(usize) -> f64
    {
        let lo: usize = 1;
        let hi: usize = self.upper_bound_n();

        // Start from n=1 and exhaustively test upwards, since the upper bound hi is often not very large, and computing the security level for larger n is expensive.
        debug!("Searching optimal n in [{lo}, {hi}] for {self} with target lambda={lambda}...");
        for n in lo..=hi {
            let curr = self.with_n(n).with_length_bound(length_bound(n));
            let lambda_curr = curr.security_level();
            debug!("\t{curr} -> {lambda_curr} bits of security");
            if lambda_curr >= lambda as f64 {
                return Ok(n);
            }
        }
        Err(LatticeEstimatorError::from(format!("no suitable n found in range [{lo}, {hi}] for {self}")))
    }
}

#[cfg(test)]
mod test {
    use crate::norms::Norm;

    use super::*;

    const Q: u64 = 2147483649;
    const SQRT_Q: f64 = 46340.95;
    const TEST_L2: MSIS = MSIS { n: 5, d: 64, q: Q, length_bound: SQRT_Q, m: 512, norm: Norm::L2 };
    const TEST_LINF: MSIS = MSIS { n: 2, d: 64, q: Q, length_bound: 1., m: 512, norm: Norm::Linf };

    #[test]
    fn test_msis_security_level_l2()
    {
        let lambda = TEST_L2.security_level();
        println!("{TEST_L2} -> lambda: {lambda}");
    }

    #[test]
    fn test_msis_security_level_linf()
    {
        let lambda = TEST_LINF.security_level();
        println!("{TEST_LINF} -> lambda: {lambda}");
    }

    #[test]
    fn test_find_optimal_n_l2()
    {
        let n_opt = TEST_L2.with_n(0).find_optimal_n(128).unwrap();
        let msis = TEST_L2.with_n(n_opt);
        println!("{TEST_L2} -> lambda: {}", TEST_L2.security_level());
        println!("{msis} -> lambda: {}", msis.security_level());
    }

    #[test]
    fn test_find_optimal_n_linf()
    {
        let n_opt = TEST_LINF.with_n(0).find_optimal_n(128).unwrap();
        let msis = TEST_LINF.with_n(n_opt);
        println!("{TEST_L2} -> lambda: {}", TEST_LINF.security_level());
        println!("{msis} -> lambda: {}", msis.security_level());
    }
}
