use std::fmt;
use std::fmt::{Debug, Display};

use core::{f64, panic};


use num_bigint::BigUint;
use num_traits::ToPrimitive;

use crate::errors::LatticeEstimatorError;
use crate::norms::Norm;
use crate::reduction::Estimates;
use statrs::function::factorial::factorial;
use crate::sis::SIS;

pub struct KSIS {
    pub h: usize,
    pub k: usize,
    pub q: BigUint,
    pub length_bound: f64,
    pub sigma: f64,
    pub w: usize,
    pub norm: Norm
}

impl Display for KSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-SIS[h={}, w={}, k={}, q={}, length_bound={}, sigma={}, norm={}]",
            self.h, self.w, self.k, self.q, self.length_bound, self.sigma, self.norm
        )
    }
}

impl Debug for KSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-SIS[h={}, w={}, k={}, q={}, length_bound={}, sigma={}, norm={}]",
            self.h, self.w, self.k, self.q, self.length_bound, self.sigma, self.norm
        )
    }
}

impl KSIS {
    pub fn with_h(&self, h: usize) -> Self {
        KSIS {
            h,
            k: self.k,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            norm: self.norm,
        }
    }

    pub fn with_length_bound(&self, length_bound: f64) -> Self {
        KSIS {
            h: self.h,
            k: self.k,
            w: self.w,
            q: self.q.clone(),
            length_bound,
            sigma: self.sigma,
            norm: self.norm,
        }
    }

    pub fn with_k(&self, k: usize) -> Self {
        KSIS {
            h: self.h,
            k,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            norm: self.norm,
        }
    }

    pub fn with_sigma(&self, sigma: f64) -> Self {
        KSIS {
            h: self.h,
            k: self.k,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma,
            norm: self.norm,
        }
    }

    pub fn to_sis(&self) -> SIS {
        let t = (self.h as f64).log2().sqrt();
        let new_bound: f64 = self.length_bound * (self.k as f64).powf(3.0 / 2.0) * factorial(self.k as u64) as f64 * (self.sigma * t * (self.h as f64).log2()).powf(self.k as f64);

        SIS::new(
            self.h,
            self.q.clone(),
            new_bound,
            self.w - self.k,
            self.norm,
        )
    }

    pub fn security_level(&self) -> f64 {

        self.to_sis().security_level()
    }  

    pub fn security_level_internal(&self, estimate_type: Estimates) -> Result<f64, LatticeEstimatorError> {

        if (self.w as f64) < (2.0 as f64) * self.h as f64 * self.q.to_f64().unwrap().log2(){
            return Err(LatticeEstimatorError::InvalidParameter { param_name : self.w.to_string() , reason : "w should be bigger than 2hlog_2(q)".to_string() });
        }

        if self.w / self.k <= self.h {
            return Err(LatticeEstimatorError::InvalidParameter { param_name : self.k.to_string() , reason : "h should be smaller than w/k".to_string() });
        }

        if (self.sigma as f64) < ((self.w as f64).log2()).sqrt() {
            return Err(LatticeEstimatorError::InvalidParameter { param_name : self.sigma.to_string() , reason : "sigma should be bigger than sqrt(log(w))".to_string() });
        }

        if self.q.to_f64().unwrap() < self.sigma as f64 * ((self.w as f64).log2()).sqrt() {
            return Err(LatticeEstimatorError::InvalidParameter { param_name : self.q.to_string() , reason : "q should be bigger than sigma*sqrt(log(w))".to_string() });
        }

        self.to_sis().security_level_internal(estimate_type)
    }
}



mod test{
    
    use num_bigint::{ToBigInt, ToBigUint};
    use num_traits::FromPrimitive;

    use crate::norms::Norm;

    use super::*;

    const Q: u128 = 2147483649;
    const SQRT_Q: f64 = 46340.95;

    //TODO find applicable tests
    
    #[test]
    fn test_rsis_security_level_l2() {
        let test_l2: KSIS = KSIS {
            h: 512,
            q: BigUint::from_f64((2.0 as f64).powf(170.0)).unwrap(),
            k: 10,
            length_bound: 5833.9072,
            w: 1024,
            sigma: 8.0,
            norm: Norm::L2,
        };

        print!("{test_l2}");

        let lambda: f64 = test_l2.security_level_internal(Estimates::CheNgueEnum).unwrap();
        println!("Internal  CheNgue12 {test_l2} -> lambda: {lambda}");
    }
}