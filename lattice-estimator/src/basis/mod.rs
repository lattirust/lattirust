use std::fmt;
use std::fmt::{Debug, Display};

use num_bigint::BigUint;
use crate::errors::LatticeEstimatorError;
use crate::norms::Norm;
use crate::sis::SIS;
use crate::reduction::Estimates;

pub struct BASISrand {
    h: usize,
    w: usize,
    q: BigUint,
    length_bound: f64,
    sigma: f64,
    k: usize,
    norm: Norm,
}

impl Display for BASISrand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "BASIS_rand[h={}, w={}, q={}, length_bound={}, sigma={}, k={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.sigma, self.k, self.norm
        )
    }
}

impl Debug for BASISrand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "BASIS_rand[h={}, w={}, q={}, length_bound={}, sigma={}, k={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.sigma, self.k, self.norm
        )
    }
}

impl BASISrand {


    pub fn new(h: usize, w: usize, q: BigUint, length_bound: f64, sigma: f64, k: usize, norm: Norm) -> Self {
        BASISrand {
            h,
            w,
            q,
            length_bound,
            sigma,
            k,
            norm,
        }
    }

    pub fn with_h(&self, h: usize) -> Self {
        BASISrand {
            h,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_w(&self, w: usize) -> Self {
        BASISrand {
            h: self.h,
            w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_length_bound(&self, length_bound: f64) -> Self {
        BASISrand {
            h: self.h,
            w: self.w,
            q: self.q.clone(),
            length_bound,
            sigma: self.sigma,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_k(&self, k: usize) -> Self {
        BASISrand {
            h: self.h,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            k,
            norm: self.norm,
        }
    }


    pub fn to_sis(&self) -> SIS {
        SIS::new(
            self.h,
            self.w,
            self.q.clone(),
            self.length_bound,
            self.norm)
    }

    pub fn security_level_internal(&self, estimate_type: Estimates) -> Result<f64, LatticeEstimatorError> { 

        if (self.sigma as f64) < (self.k as f64 * self.w as f64 * (self.h as f64 * self.k as f64).log2()) {
            return Err(LatticeEstimatorError::InvalidParameter {
                param_name: "sigma".to_string(),
                reason: "sigma < k*w*h*log2(h*k)".to_string(),
            });

        } else {   
            self.to_sis().security_level_internal(estimate_type)
        }
    }

    pub fn upper_bound_h(&self) -> usize {
        self.to_sis().upper_bound_h()
    }


    pub fn find_optimal_h_annealing(&self, lambda: usize, estimate_type: Estimates) -> usize {
        match self.to_sis().find_optimal_h_annealing(lambda, estimate_type) {
            Ok(h) => h,
            //print the lattice estimator error
            Err(e) => panic!("Error: {:?}", e)
        }
    }

    pub fn find_optimal_length_bound_annealing(&self, lambda: usize, estimate_type: Estimates) -> f64 {
        match self.to_sis().find_optimal_length_bound_annealing(lambda, estimate_type) {
            Ok(length_bound) => length_bound,
            Err(e) => panic!("Error: {:?}", e)
        }
    }
}



#[cfg(test)]

mod test{

    use crate::norms::Norm;
    use crate::reduction::Estimates;

    use super::*;

    #[test]
    fn test_basis_security_level_l2() {
        let falcon512_unf: BASISrand = BASISrand::new(512, 1024, 12289u64.into(), 5833.9072, 235332.0, 1, Norm::L2);

        let security_level = falcon512_unf.security_level_internal(Estimates::Matzov(true)).unwrap();
        print!("Security level: {}", security_level);

    }

    #[test]
    fn test_find_optimal_h(){
        let falcon512_unf: BASISrand = BASISrand::new(512, 1024, 12289u64.into(), 5833.9072, 235332.0, 1, Norm::L2);
        let h_optimal = falcon512_unf.find_optimal_h_annealing(122, Estimates::Matzov(true));
        println!("Optimal h: {}", h_optimal);

        let security_level = falcon512_unf.with_h(h_optimal).security_level_internal(Estimates::Matzov(true)).unwrap();
        println!("Security level: {}", security_level);
    }

    #[test]
    fn test_find_optimal_lenbound(){

        let falcon512_unf: BASISrand = BASISrand::new(512, 1024, 12289u64.into(), 5833.9072, 235332.0, 1, Norm::L2);
        let lenbound_optimal = falcon512_unf.find_optimal_length_bound_annealing(122, Estimates::Matzov(true));
        println!("Optimal length bound: {}", lenbound_optimal);

        let security_level = falcon512_unf.with_length_bound(lenbound_optimal).security_level_internal(Estimates::Matzov(true)).unwrap();
        println!("Security level: {}", security_level);
    }
}