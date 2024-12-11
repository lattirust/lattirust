use std::fmt;
use std::fmt::{Debug, Display};

use num_bigint::BigUint;

use crate::norms::Norm;
use crate::rsis::RSIS;
use crate::sis::SIS;
use crate::errors::LatticeEstimatorError;
use crate::reduction::Estimates;
use crate::msis::{self, MSIS};

pub struct KMSIS {
    pub h: usize,
    pub w: usize,
    pub d: usize,
    pub q: BigUint,
    pub length_bound: f64,
    pub sigma: f64,
    pub k: usize,
    pub norm: Norm,
}

impl Display for KMSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-M-SIS[h={}, w={}, d={}, q={}, length_bound={}, sigma={}, k={}, norm={}]",
            self.h, self.w, self.d, self.q, self.length_bound, self.sigma, self.k, self.norm
        )
    }
}

impl Debug for KMSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-M-SIS[h={}, w={}, d={}, q={}, length_bound={}, sigma={}, k={}, norm={}]",
            self.h, self.w, self.d, self.q, self.length_bound, self.sigma, self.k, self.norm
        )
    }
}

impl KMSIS {
    pub fn new(h: usize, w: usize, d: usize, q: BigUint, length_bound: f64, sigma: f64, k: usize, norm: Norm) -> Self {
        KMSIS {
            h,
            w,
            d,
            q,
            length_bound,
            sigma,
            k,
            norm,
        }
    }

    pub fn with_h(&self, h: usize) -> Self {
        KMSIS {
            h,
            w: self.w,
            d: self.d,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_w(&self, w: usize) -> Self {
        KMSIS {
            h: self.h,
            w,
            d: self.d,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_length_bound(&self, length_bound: f64) -> Self {
        KMSIS {
            h: self.h,
            w: self.w,
            d: self.d,
            q: self.q.clone(),
            length_bound,
            sigma: self.sigma,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_sigma(&self, sigma: f64) -> Self {
        KMSIS {
            h: self.h,
            w: self.w,
            d: self.d,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_k(&self, k: usize) -> Self {
        KMSIS {
            h: self.h,
            w: self.w,
            d: self.d,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            k,
            norm: self.norm,
        }
    }
    
    pub fn upper_bound_h(&self) -> usize {
        let msis: MSIS = MSIS::new(self.h, self.w - self.k, self.d, self.q.clone(), self.length_bound, self.norm);
        msis.upper_bound_h()
    }

    pub fn to_sis(&self) -> SIS {        
        let msis: MSIS = MSIS::new(self.h, self.w - self.k, self.d, self.q.clone(), self.length_bound, self.norm);
        msis.to_sis()
    }

    pub fn security_level_internal(&self, estimate_type: Estimates) -> Result<f64, LatticeEstimatorError> {
        if self.k >= self.w {
            return Err(LatticeEstimatorError::InvalidParameter { param_name : self.k.to_string() , reason : "k should be smaller than w, otherwise reduction does not work".to_string() });
        }
        self.to_sis().security_level_internal(estimate_type)
    }

    pub fn find_optimal_h_annealing(&self, lambda: usize, estimate_type: Estimates) -> usize {
        let msis: MSIS = MSIS::new(self.h, self.w - self.k, self.d, self.q.clone(), self.length_bound, self.norm);
        return msis.find_optimal_h_annealing(lambda, estimate_type);
    }

    pub fn find_optimal_length_bound_annealing(&self, lambda: usize, estimate_type: Estimates) -> f64 {
        match self.to_sis().find_optimal_length_bound_annealing(lambda, estimate_type) {
            Ok(length_bound) => length_bound,
            Err(e) => panic!("Error: {:?}", e)
        }
    }
}

pub struct KRSIS {
    pub h: usize,
    pub w: usize,
    pub q: BigUint,
    pub length_bound: f64,
    pub sigma: f64,
    pub k: usize,
    pub norm: Norm,
}

impl Display for KRSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-R-SIS[h={}, w={}, q={}, length_bound={}, sigma={}, k={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.sigma, self.k, self.norm
        )
    }
}

impl Debug for KRSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-R-SIS[h={}, w={}, q={}, length_bound={}, sigma={}, k={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.sigma, self.k, self.norm
        )
    }
}


impl KRSIS {

    pub fn new(h: usize, w: usize, q: BigUint, length_bound: f64, sigma: f64, k: usize, norm: Norm) -> Self {
        KRSIS {
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
        KRSIS {
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
        KRSIS {
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
        KRSIS {
            h: self.h,
            w: self.w,
            q: self.q.clone(),
            length_bound,
            sigma: self.sigma,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_sigma(&self, sigma: f64) -> Self {
        KRSIS {
            h: self.h,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_k(&self, k: usize) -> Self {
        KRSIS {
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

        RSIS::new(
            self.h,
            self.w - self.k,
            self.q.clone(),
            self.length_bound,
            self.norm,
        ).to_sis()

    }

    pub fn upper_bound_h(&self) -> usize {
        self.to_sis().upper_bound_h()
    }

    pub fn security_level_internal(&self, estimate_type: Estimates) -> Result<f64, LatticeEstimatorError> {
        self.to_sis().security_level_internal(estimate_type)
    }

    pub fn find_optimal_h_annealing(&self, lambda: usize, estimate_type: Estimates) -> usize {
        match self.to_sis().find_optimal_h_annealing(lambda, estimate_type) {
            Ok(h) => h,
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
    use crate::{kmsis::KMSIS, norms::Norm, reduction::Estimates};

    
    #[test]
    fn test_sis_security_level_l2() {
        let falcon512_unf: KMSIS = KMSIS::new(24, 512,4, 122898899u64.into(), 5833.9072, 4.0, 2,  Norm::L2);


        let cost_internal:Result<f64, crate::errors::LatticeEstimatorError> = falcon512_unf.security_level_internal(Estimates::Kyber(true));

        match cost_internal {
            Ok(cost) => {
                println!("Internal : {falcon512_unf} -> lambda: {cost}");
            },
            Err(e) => {
                panic!("Error: {:?}", e);
            }
        }
    
    }

    #[test]
    fn test_sis_security_level_linf() {
        let falcon512_unf: KMSIS = KMSIS::new(24, 512,4, 122898899u64.into(), 5833.9072, 4.0, 2,  Norm::Linf);

        let cost_internal:Result<f64, crate::errors::LatticeEstimatorError> = falcon512_unf.security_level_internal(Estimates::Kyber(true));

        match cost_internal {
            Ok(cost) => {
                println!("Internal : {falcon512_unf} -> lambda: {cost}");
            },
            Err(e) => {
                panic!("Error: {:?}", e);
            }
            
        }
    
    }


    #[test]
    fn test_sis_h_annealing() {
        let falcon512_unf: KMSIS = KMSIS::new(24, 512,4, 122898899u64.into(), 5833.9072, 4.0, 2,  Norm::L2);
        let h_opt = falcon512_unf.with_h(falcon512_unf.upper_bound_h()).find_optimal_h_annealing(64, Estimates::Kyber(true));
        let lambda_opt = falcon512_unf.with_h(h_opt).security_level_internal(Estimates::Kyber(true)).unwrap();

        println!("Optimal h: {h_opt} -> lambda: {lambda_opt}");

    }

    #[test]
    fn test_sis_lenbound_annealing() {

        let falcon512_unf: KMSIS = KMSIS::new(24, 512,4, 122898899u64.into(), 5833.9072, 4.0, 2,  Norm::L2);


        let len_opt = falcon512_unf.find_optimal_length_bound_annealing(64, Estimates::Kyber(true));
        let lambda_opt = falcon512_unf.with_length_bound(len_opt).security_level_internal(Estimates::Kyber(true)).unwrap();

        println!("Optimal length bound: {len_opt} -> lambda: {:?}", lambda_opt);

        let lambda_test = falcon512_unf.with_length_bound(40000.0).security_level_internal(Estimates::Kyber(true)).unwrap();
        println!("Test length bound: 40000 -> lambda: {:?}", lambda_test);
    }
}