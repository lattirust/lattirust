use core::{f64, panic};
use std::fmt;
use std::fmt::{Debug, Display};
use std::num::ParseFloatError;
use std::str::FromStr;

use num_bigint::BigUint;
use num_traits::ToPrimitive;

use crate::utils::CostParameters;
use crate::errors::LatticeEstimatorError;
use crate::norms::Norm;
use crate::sage_util::sagemath_eval;
use crate::reduction::matzov_short_vectors;
use crate::reduction::Estimates;
use crate::{reduction, utils};
use statrs::distribution::{ContinuousCDF, Normal};

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

    fn cost_infinity(&self, zeta: usize, block_size: usize, classical: bool, estimate_type: Estimates) -> CostParameters{
        
        let q: f64 = self.q.to_f64().unwrap();
        let dim: usize = self.w - zeta;

        if dim < block_size || dim < self.h {
            return CostParameters::with_values(f64::INFINITY, f64::INFINITY, 0.0, block_size, 0,  zeta, dim);
        }

        // Step 1: Simulate Gram-Schmidt vector lengths using reduction simulator
        let gs_vectors_sizes =  reduction::gsa_simulator(dim, dim - self.h, q as usize, block_size, None);

        let sampling_cost: (f64, f64, usize, usize) = matzov_short_vectors(block_size, dim, q as usize, Some(dim), classical);
        let bkz_cost: f64 = reduction::bkz_cost(estimate_type, block_size, dim, q as usize);
        
        let mut log_proba: f64;
        if (self.w as f64).sqrt() * self.length_bound <= q {

            let vec_len = sampling_cost.0 * gs_vectors_sizes[0].sqrt();
            let std_dev: f64 = vec_len / (dim as f64).sqrt();
            let cdf = Normal::new(0.0, std_dev).unwrap().cdf(-self.length_bound);
            log_proba = (1.0 - 2.0 * cdf).log2() * dim as f64;
        
        } else {
            
            let idx_start = if (gs_vectors_sizes[0].sqrt() - q).abs() < 1e-8 {
                gs_vectors_sizes.iter().position(|&r| r < gs_vectors_sizes[0]).unwrap_or(0)
            } else {
                0
            };
    
            let idx_end = if (gs_vectors_sizes[gs_vectors_sizes.len() - 1] - 1.0).abs() < 1e-8 {
                gs_vectors_sizes.iter().position(|&r| r.sqrt() <= 1.0 + 1e-8).unwrap_or(dim - 1)
            } else {
                dim - 1
            };
    
            let vector_length = gs_vectors_sizes[idx_start].sqrt();
            let gaussian_coords = (idx_end - idx_start + 1).max(block_size);
            let std_dev = vector_length / (gaussian_coords as f64).sqrt();
            let cdf = Normal::new(0.0, std_dev).unwrap().cdf(-self.length_bound);
            log_proba = (1.0 - 2.0 * cdf).log2() * gaussian_coords as f64;
            log_proba += ((2.0 * self.length_bound + 1.0) / q).log2() * idx_start as f64;
        }
        
        let proba = (2.0 as f64).powf(f64::min(0.0, log_proba + (sampling_cost.2 as f64).log2()));

        let amplificator: f64 = utils::amplify_via_trials(0.99, proba);
        if amplificator == f64::INFINITY {
            return CostParameters::with_values(f64::INFINITY, f64::INFINITY, 0.0, block_size, 0, zeta, dim);
        }
        return CostParameters::with_values(
            (sampling_cost.1 * amplificator).log2(),
            (bkz_cost * amplificator).log2(),
            (f64::max(sampling_cost.1 - bkz_cost, 1e-100) * amplificator).log2(), 
            block_size, 
            sampling_cost.3, 
            zeta, 
            dim);

    }


    fn optimal_w(&self) -> f64 {

        //see how to deal with BigUint and f64 conversion
        let q: f64 = self.q.to_f64().unwrap();
        let log_delta = self.length_bound.log2().powf(2.0) / (4.0 * self.h as f64 * q.log2());
        let optimal_w: f64 = (self.h as f64 * q.log2() / log_delta).sqrt();
        optimal_w 
    }

    fn optimal_block_size(&self, delta: f64) -> usize {
        let mut beta: usize = 40; 

        while reduction::bkz_delta(2 * beta) > delta {
            beta *= 2;
        }
        while reduction::bkz_delta(beta + 10) > delta {
            beta += 10;
        }
        while reduction::bkz_delta(beta) >= delta {
            beta += 1;
        }
        beta
    }


    pub fn security_level_internal(&self, estimate_type: Estimates) -> f64 {
        let q = self.q.to_f64().unwrap();
        

        //Check for trivial case and impossible case
        //0.1 is arbitrary see if it can be changed
        if self.length_bound < (self.w as f64 * q.log2()).sqrt() || q < self.length_bound * (self.h as f64).powf(0.1) {
            panic!("The length bound is too small for the given parameters, no short vector is expected to be found.");
        } else if self.length_bound >= q {
            panic!("The solution is trivial for the given parameters. Please set the norm < q");
        } else {

            match self.norm {
                Norm::L2 => {
                    
                    //find optimal lattice dimension for reduction
                    let w_opt: f64 = f64::min(self.optimal_w(), self.w as f64);

                    //find the root hermite factor required for the optimal lattice dimension
                    let log_delta: f64 = (1.0 / (w_opt - 1.0)) * (self.length_bound.log2() - (self.h as f64 / w_opt) * q.log2());
                    let delta: f64 = log_delta.exp2();

                    //Find optimal BKZ block size
                    let block_size: usize = self.optimal_block_size(delta);  
                    
                    //Actual reduction
                    if delta >= 1.0 && block_size <= w_opt.round() as usize {
                        return reduction::bkz_cost(estimate_type, block_size, self.h, q as usize).log2();
                    }
                    else {
                        panic!("The reduction is not doable with the given parameters");
                    }
                }
                
                Norm::Linf => {

                    // Define the range for zeta 
                    let zeta_range = (1, self.w); 
                    
                    //Get the L2 upperbound on blocksize
                    let w_opt: f64 = f64::min(self.optimal_w(), self.w as f64);
                    let log_delta: f64 = (1.0 / (w_opt - 1.0)) * (self.length_bound.log2() - (self.h as f64 / w_opt) * q.log2());
                    let delta: f64 = log_delta.exp2();
                    let l2_block_bound: usize = self.optimal_block_size(delta);  
                   
                    // Define the range for block size
                    let block_size_range = (40, l2_block_bound);

                    // Run grid search to find optimal zeta and block_size that minimize cost
                    let (_best_zeta, _best_beta, best_cost) = utils::adaptive_grid_search(
                        |zeta, block_size| self.cost_infinity(zeta, block_size, true, estimate_type),
                        zeta_range,
                        block_size_range,
                        10
                    );
                    
                    // Update best cost value based on found parameters
                    let final_cost = best_cost.rop;

                    // Return the security level corresponding to the minimal found cost
                    return final_cost;
                }
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
        let log2_m: u32 = self.w.ilog2();
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
    use statrs::assert_almost_eq;

    use crate::norms::Norm;
    use crate::sis::SIS;
    use crate::reduction::Estimates;

    #[test]
    fn test_sis_security_level_l2() {
        let falcon512_unf: SIS = SIS::new(512, 12289u64.into(), 5833.9072, 1024, Norm::L2);
        let lambda = falcon512_unf.security_level();
        assert!(lambda >= 128.);
        println!("External : {falcon512_unf} -> lambda: {lambda}");

        let cost_internal: f64 = falcon512_unf.security_level_internal(Estimates::Kyber(true));
        println!("Internal : {falcon512_unf} -> lambda: {cost_internal}");
    }

    #[test]
    fn test_cost_models_l2() {
        let falcon512_unf: SIS = SIS::new(512, 12289u64.into(), 5833.9072, 1024, Norm::L2);
        
        let cost_chengue12: f64 = falcon512_unf.security_level_internal(Estimates::CheNgueEnum);
        println!("Cost CheNgue12: {falcon512_unf} -> lambda: {cost_chengue12}");

        let cost_bdgl16: f64 = falcon512_unf.security_level_internal(Estimates::BdglSieve);
        println!("Cost BDGL16: {falcon512_unf} -> lambda: {cost_bdgl16}");

        let cost_laamospol14: f64 = falcon512_unf.security_level_internal(Estimates::QSieve);
        println!("Cost LaaMosPol14: {falcon512_unf} -> lambda: {cost_laamospol14}");

        let cost_abfksw20: f64 = falcon512_unf.security_level_internal(Estimates::AbfEnum(true));
        println!("Cost ABFKSW20: {falcon512_unf} -> lambda: {cost_abfksw20}");

        let cost_ablr21: f64 = falcon512_unf.security_level_internal(Estimates::AblrEnum);
        println!("Cost ABLR21: {falcon512_unf} -> lambda: {cost_ablr21}");

        let cost_adps16: f64 = falcon512_unf.security_level_internal(Estimates::AdpsSieve(true));
        println!("Cost ADPS16: {falcon512_unf} -> lambda: {cost_adps16}");

        let cost_chaloy21: f64 = falcon512_unf.security_level_internal(Estimates::ChaLoySieve);
        println!("Cost ChaLoy21: {falcon512_unf} -> lambda: {cost_chaloy21}");

        let cost_kyber: f64 = falcon512_unf.security_level_internal(Estimates::Kyber(true));
        println!("Cost Kyber: {falcon512_unf} -> lambda: {cost_kyber}");
    }

    #[test]
    fn test_sis_security_level_linf() {
        let dilithium2_msis_wk_unf: SIS =
            SIS::new(1024, 8380417u64.into(), 350209., 2304, Norm::Linf);

        let lambda: f64 = dilithium2_msis_wk_unf.security_level();
        assert!(lambda >= 128.);
        println!("External {dilithium2_msis_wk_unf} -> lambda: {lambda}");

        let cost_internal: f64 = dilithium2_msis_wk_unf.security_level_internal(Estimates::AbfEnum(true));
        println!("Internal Linf {dilithium2_msis_wk_unf} -> lambda: {cost_internal}");
    }

    #[test]
    fn test_find_optimal_h_l2() {
        let falcon512_unf: SIS = SIS::new(512, 12289u64.into(), 5833.9072, 1024, Norm::L2);
        let h_opt = falcon512_unf.find_optimal_h(128).unwrap();
        let sis = falcon512_unf.with_h(h_opt);
        println!("{} ", h_opt);
        let lambda = sis.security_level();
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

    #[test]
    fn test_cost_infinity() {
        let dilithium2_msis_wk_unf: SIS =
            SIS::new(1024, 8380417u64.into(), 350209., 2304, Norm::Linf);
        
        let cost = dilithium2_msis_wk_unf.cost_infinity(575, 776, true, Estimates::Matzov(true));
        assert_almost_eq!(cost.rop, 248.0, 1e1);

        let cost = dilithium2_msis_wk_unf.cost_infinity(2303, 224, true, Estimates::Matzov(true));
        assert_almost_eq!(cost.rop, f64::INFINITY, 1e1);

        let cost = dilithium2_msis_wk_unf.cost_infinity(1152, 42, true, Estimates::Matzov(true));
        assert_almost_eq!(cost.rop, f64::INFINITY, 1e1);
    }
}
