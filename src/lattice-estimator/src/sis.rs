use std::fmt;
use std::fmt::{Debug, Display};
use std::num::ParseFloatError;
use std::str::FromStr;

use crate::errors::LatticeEstimatorError;
use crate::norms::Norm;
use crate::sage_util::sagemath_eval;

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

    pub const fn with_length_bound(&self, length_bound: f64) -> Self {
        SIS { n: self.n, q: self.q, length_bound, m: self.m, norm: self.norm }
    }

    pub fn parse_f64(s: String) -> Result<f64, ParseFloatError> {
        // Both lattice-estimator and security-estimator logs some additional info, we only care about the last line of stdout
        f64::from_str(&s.lines().last().unwrap())
    }

    /// Return lambda such that SIS_{n, q, length_bound, m} is 2^lambda-hard (for a given norm).
    /// Internally, this calls out to the lattice-estimator via a wrapper Python script.
    pub fn security_level(&self) -> f64 {
        let func = match self.norm {
            Norm::L2 => "sis_security_level_l2",
            Norm::Linf => "sis_security_level_linf"
        };
        sagemath_eval(format!("{}({}, {}, {}, {})", func, self.n, self.q, self.length_bound, self.m), SIS::parse_f64).unwrap()
    }

    pub fn upper_bound_n(&self) -> usize {
        // let log_q = match self.norm {
        //     Norm::L2 => (self.q as f64).log2(),
        //     Norm::Linf => (self.q as f64).log(2. * self.length_bound + 1.)
        // };
        // (self.m as f64 / (2. * log_q)).ceil() as usize
        self.m
    }

    /// Return the smallest n such that SIS_{n, q, length_bound, m} is 2^lambda-hard (for a given norm).
    pub fn find_optimal_n(&self, lambda: usize) -> Result<usize, LatticeEstimatorError> {
        let mut hi: usize = self.upper_bound_n();
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

    /// Return the smallest n such that SIS_{n, q, length_bound(n), m} is 2^lambda-hard (for a given norm), where length_bound is a function of n.
    pub fn find_optimal_n_dynamic<F>(&self, length_bound: F, lambda: usize) -> Result<usize, LatticeEstimatorError>
        where F: Fn(usize) -> f64
    {
        // We can't assume that length_bound is monotonic, so we can't use binary search.
        // Instead, exhaustively search powers of 2 until we find a suitable n.
        // TODO: use a better search algorithm / return a more fine-grained result
        let log2_m = self.m.ilog2();
        let candidates = (1..log2_m).map(|i| 2usize.pow(i));
        candidates.map(|n| self.with_n(n).with_length_bound(length_bound(n)))
            .find(|sis| sis.security_level() >= lambda as f64)
            .map(|sis| sis.n).ok_or(LatticeEstimatorError::from("no suitable n found".to_string()))
    }
}


#[cfg(test)]
mod test {
    use crate::norms::Norm;
    use crate::sis::SIS;

    // from lattice-estimator/schemes
    const FALCON512_UNF: SIS = SIS::new(512, 12289, 5833.9072, 1024, Norm::L2);
    const DILITHIUM2_MSIS_WK_UNF: SIS = SIS::new(1024, 8380417, 350209., 2304, Norm::Linf);

    #[test]
    fn test_sis_security_level_l2()
    {
        let lambda = FALCON512_UNF.security_level();
        assert!(lambda >= 128.);
        println!("{FALCON512_UNF} -> lambda: {lambda}");
    }

    #[test]
    fn test_sis_security_level_linf()
    {
        let lambda = DILITHIUM2_MSIS_WK_UNF.security_level();
        assert!(lambda >= 128.);
        println!("{DILITHIUM2_MSIS_WK_UNF} -> lambda: {lambda}");
    }

    #[test]
    fn test_find_optimal_n_l2()
    {
        let n_opt = FALCON512_UNF.find_optimal_n(128).unwrap();
        let sis = FALCON512_UNF.with_n(n_opt);
        let lambda = sis.security_level();
        assert!(lambda >= 128.0);
        println!("{FALCON512_UNF} -> lambda: {}", FALCON512_UNF.security_level());
        println!("{sis} -> lambda: {lambda}");
    }

    #[test]
    fn test_find_optimal_n_linf()
    {
        let n_opt = DILITHIUM2_MSIS_WK_UNF.find_optimal_n(128).unwrap();
        let sis = DILITHIUM2_MSIS_WK_UNF.with_n(n_opt);
        let lambda = sis.security_level();
        assert!(lambda >= 128.0);
        println!("{DILITHIUM2_MSIS_WK_UNF} -> lambda: {}", DILITHIUM2_MSIS_WK_UNF.security_level());
        println!("{sis} -> lambda: {lambda}");
    }
}
