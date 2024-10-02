use std::f64::consts::{E, PI};
use std::fmt;
use std::fmt::{Debug, Display};
use std::num::ParseFloatError;
use std::str::FromStr;

use num_bigint::BigUint;
use num_traits::ToPrimitive;
use roots::SimpleConvergency;

use crate::errors::LatticeEstimatorError;
use crate::norms::Norm;
use crate::sage_util::sagemath_eval;
use crate::reduction::EstimatesParams;
use crate::reduction::Estimates;
use crate::reduction;

pub struct SIS {
    h: usize,
    w: usize,
    q: BigUint,
    length_bound: f64,
    norm: Norm,
}

impl Display for SIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "SIS[h={}, w={}, q={}, length_bound={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.norm
        )
    }
}

impl Debug for SIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "SIS[h={}, w={}, q={}, length_bound={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.norm
        )
    }
}

impl SIS {

    pub const fn new(h: usize, q: BigUint, length_bound: f64, w: usize, norm: Norm) -> Self {
        SIS {
            h,
            w,
            q,
            length_bound,
            norm,
        }
    }

    pub fn with_h(&self, h: usize) -> Self {
        SIS {
            h,
            q: self.q.clone(),
            length_bound: self.length_bound,
            w: self.w,
            norm: self.norm,
        }
    }

    pub fn with_length_bound(&self, length_bound: f64) -> Self {
        SIS {
            h: self.h,
            q: self.q.clone(),
            length_bound,
            w: self.w,
            norm: self.norm,
        }
    }

    fn optimal_w(&self) -> f64 {
        let q: f64 = self.q.to_f64().unwrap();
        let log_delta = self.length_bound.log2().powi(2) / (4.0 * self.h as f64 * q.log2());
        (self.h as f64 * q.log2() / log_delta).sqrt() 
    }

    fn optimal_block_size(&self, delta: f64) -> Result<f64, roots::SearchError> {
        let f = |beta| {beta / (2.0 * PI * E) * (PI * beta as f64).powf(1.0 / beta).powf(1.0 / (2.0 * (beta - 1.0)))- delta};
        let mut conv: SimpleConvergency<f64> = SimpleConvergency {eps: 1e-15, max_iter: 100};
        roots::find_root_secant(f64::powf(2.0, 16.0), 40f64, f, &mut conv)
    }

    pub fn security_level_internal(&self, estimate_type: Estimates) -> f64 {
        let q = self.q.to_f64().unwrap();
        //Check for trivial case and impossible case
        if self.length_bound < (self.w as f64 * q.log2()).sqrt() {
            panic!("The length bound is too small for the given parameters");
        } else if self.length_bound >= q {
            panic!("The solution is trivial for the given parameters");
        } else {

            //find optimal lattice dimension for reduction
            let w_opt = f64::min(self.optimal_w(), self.w as f64);
            println!("Optimal w: {}", w_opt);

            //find the root hermite factor required for the optimal lattice dimension
            let log_delta: f64 = (1.0 / (w_opt - 1.0)) * (self.length_bound.log2() - (self.h as f64 / w_opt) * q.log2());
            let delta: f64 = log_delta.exp2();

            //For now use the rule of thumb in "On the concrete hardness with learning with errors"
            let block_size_thumb: usize = (1.0 / log_delta).round() as usize; 
            println!("Block size thumb: {}", block_size_thumb);
            
            //find optimal block size
            let block_size: usize = match self.optimal_block_size(delta) {
                Ok(bl_size) => bl_size as usize,
                Err(e) => panic!("Error finding optimal block size: {}", e)
            }; 

            println!("Block size: {}", block_size);  
            
            if delta >= 1.0 && block_size <= w_opt.round() as usize {
                let estimate_params: EstimatesParams = match estimate_type {
                    Estimates::LLL => EstimatesParams::LLL(reduction::LLL::new(self.h, q.to_bits() as usize)),
                    Estimates::CheNgue12=> EstimatesParams::CheNgue12(reduction::BKZ::new(block_size, self.h, q.to_bits() as usize))
                };
                reduction::cost(estimate_params).unwrap().log2()
            }
            else {
                panic!("The reduction is not doable with the given parameters");
            }
        }
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
            Norm::Linf => "sis_security_level_linf",
        };
        sagemath_eval(
            format!(
                "{}({}, {}, {}, {})",
                func, self.h, self.q, self.length_bound, self.w
            ),
            SIS::parse_f64,
        )
        .unwrap()
    }

    pub fn upper_bound_h(&self) -> usize {
        let log_q = match self.norm {
            Norm::L2 => (self.q.to_f64().unwrap()).log2(),
            Norm::Linf => (self.q.to_f64().unwrap()).log(2. * self.length_bound + 1.),
        };
        let mut h = (self.w as f64 / log_q).floor() as usize;
        // Deal with the case where e.g. w and q are powers of 2, to ensure that w > h * log_q still holds without having to set h = w / (2*log_q)
        if self.w as f64 == (h as f64) * log_q {
            h -= 1;
        }
        assert!(self.w as f64 > (h as f64) * log_q);
        h
    }

    /// Return the smallest `h` such that `SIS\[h, w, q, length_bound\]` is $2^\lambda$-hard (for a given norm).
    pub fn find_optimal_h(&self, lambda: usize) -> Result<usize, LatticeEstimatorError> {
        let mut hi: usize = self.upper_bound_h();
        let mut lo: usize = 1;

        let sis = self.with_h(hi);
        let lambda_hi = sis.security_level();
        println!(
            "{sis} has sec. param. {lambda_hi}  < target lambda = {lambda}"
        );
        // Loop invariant: SIS[hi, w, q, length_bound] is 2^lambda_hi-hard with lambda_hi >= lambda
        while hi > lo {
            let mid = lo + (hi - lo) / 2;
            let sis = self.with_h(mid);
            let mid_lambda = sis.security_level();
            if mid_lambda >= lambda as f64 {
                // Search for smaller n in [lo, mid]
                hi = mid;
            } else {
                // Search for smaller n in [mid+1, hi]
                lo = mid + 1;
            }
        }
        assert_eq!(hi, lo);
        Ok(hi)
    }

    /// Return the smallest `h` such that `SIS\[h, w, q, length_bound(h)\]` is $2^\lambda$-hard (for a given norm), where `length_bound` is a function of `h`.
    pub fn find_optimal_h_dynamic<F>(
        &self,
        length_bound: F,
        lambda: usize,
    ) -> Result<usize, LatticeEstimatorError>
    where
        F: Fn(usize) -> f64,
    {
        // We can't assume that length_bound is monotonic, so we can't use binary search.
        // Instead, exhaustively search powers of 2 until we find a suitable n.
        // TODO: use a better search algorithm / return a more fine-grained result
        let log2_m = self.w.ilog2();
        let candidates = (1..log2_m).map(|i| 2usize.pow(i));
        candidates
            .map(|n| self.with_h(n).with_length_bound(length_bound(n)))
            .find(|sis| sis.security_level() >= lambda as f64)
            .map(|sis| sis.h)
            .ok_or(LatticeEstimatorError::from(
                "no suitable h found".to_string(),
            ))
    }
}

#[cfg(test)]
mod test {
    use crate::norms::Norm;
    use crate::sis::SIS;
    use crate::reduction::Estimates;

    #[test]
    fn test_sis_security_level_l2() {
        let falcon512_unf: SIS = SIS::new(512, 12289u64.into(), 5833.9072, 1024, Norm::L2);
        let lambda = falcon512_unf.security_level();
        assert!(lambda >= 128.);
        println!("External : {falcon512_unf} -> lambda: {lambda}");


        let cost_internal_lll: f64 = falcon512_unf.security_level_internal(Estimates::LLL);
        println!("Internal LLL: {falcon512_unf} -> lambda: {cost_internal_lll}");

        let cost_internal_chengue12: f64 = falcon512_unf.security_level_internal(Estimates::CheNgue12);
        println!("Internal CheNgue12: {falcon512_unf} -> lambda: {cost_internal_chengue12}");
    }

    #[test]
    fn test_sis_security_level_linf() {
        let dilithium2_msis_wk_unf: SIS =
            SIS::new(1024, 8380417u64.into(), 350209., 2304, Norm::Linf);

        let lambda = dilithium2_msis_wk_unf.security_level();
        assert!(lambda >= 128.);
        println!("{dilithium2_msis_wk_unf} -> lambda: {lambda}");
    }

    #[test]
    fn test_find_optimal_h_l2() {
        let falcon512_unf: SIS = SIS::new(512, 12289u64.into(), 5833.9072, 1024, Norm::L2);
        let h_opt = falcon512_unf.find_optimal_h(128).unwrap();
        let sis = falcon512_unf.with_h(h_opt);
        println!("{} ", h_opt);
        let lambda = sis.security_level();
        //assert!(lambda >= 128.0);
        println!(
            "{falcon512_unf} -> lambda: {}",
            falcon512_unf.security_level()
        );
        println!("{sis} -> lambda: {lambda}");
    }

    #[test]
    fn test_find_optimal_h_linf() {
        let dilithium2_msis_wk_unf: SIS =
            SIS::new(1024, 8380417u64.into(), 350209., 2304, Norm::Linf);

        let h_opt = dilithium2_msis_wk_unf.find_optimal_h(128).unwrap();
        let sis = dilithium2_msis_wk_unf.with_h(h_opt);
        let lambda = sis.security_level();
        assert!(lambda >= 128.0);
        println!(
            "{dilithium2_msis_wk_unf} -> lambda: {}",
            dilithium2_msis_wk_unf.security_level()
        );
        println!("{sis} -> lambda: {lambda}");
    }
}
