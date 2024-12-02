use std::fmt;
use std::fmt::{Debug, Display};

use num_bigint::BigUint;

use crate::errors::LatticeEstimatorError;
use crate::norms::Norm;
use crate::sis::SIS;
use crate::reduction::Estimates;

pub struct KMISIS {
    pub h: usize,
    pub k: usize,
    pub d: usize,
    pub q: BigUint,
    pub length_bound: f64,
    pub sigma: f64,
    pub w: usize,
    pub norm: Norm,
}

impl Display for KMISIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-M-ISIS[h={}, w={}, k={}, d={}, q={}, length_bound={}, sigma={}, norm={}]",
            self.h, self.w, self.k, self.d, self.q, self.length_bound, self.sigma, self.norm
        )
    }
}

impl Debug for KMISIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-M-ISIS[h={}, w={}, k={}, d={}, q={}, length_bound={}, sigma={}, norm={}]",
            self.h, self.w, self.k, self.d, self.q, self.length_bound, self.sigma, self.norm
        )
    }
}

impl KMISIS{
    pub fn with_h(&self, h: usize) -> Self {
        KMISIS {
            h,
            k: self.k,
            w: self.w,
            d: self.d,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            norm: self.norm,
        }
    }

    pub fn with_length_bound(&self, length_bound: f64) -> Self {
        KMISIS {
            h: self.h,
            k: self.k,
            w: self.w,
            d: self.d,
            q: self.q.clone(),
            length_bound,
            sigma: self.sigma,
            norm: self.norm,
        }
    }

    pub fn with_k(&self, k: usize) -> Self {
        KMISIS {
            h: self.h,
            k,
            w: self.w,
            d: self.d,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            norm: self.norm,
        }
    }

    pub fn with_sigma(&self, sigma: f64) -> Self {
        KMISIS {
            h: self.h,
            k: self.k,
            d: self.d,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma,
            norm: self.norm,
        }
    }

    pub fn with_d(&self, d: usize) -> Self {
        KMISIS {
            h: self.h,
            k: self.k,
            w: self.w,
            d,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            norm: self.norm,
        }
    }

    pub fn to_sis(&self) -> SIS {        
        
        SIS::new(
            self.h * self.d,
            self.q.clone(),
            self.length_bound,
            (self.w + 1) * self.h * self.d,
            self.norm,
        )

    }

    pub fn security_level(&self) -> f64 {
        self.to_sis().security_level()
    }

    pub fn security_level_internal(&self, estimate_type: Estimates) -> Result<f64, LatticeEstimatorError> {
        self.to_sis().security_level_internal(estimate_type)
    }

    pub fn upper_bound_h(&self) -> usize {
        self.to_sis().upper_bound_h()
    }
}


pub struct KRISIS {
    pub h: usize,
    pub k: usize,
    pub q: BigUint,
    pub length_bound: f64,
    pub sigma: f64,
    pub w: usize,
    pub norm: Norm,
}

impl Display for KRISIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-R-ISIS[h={}, w={}, k={}, q={}, length_bound={}, sigma={}, norm={}]",
            self.h, self.w, self.k, self.q, self.length_bound, self.sigma, self.norm
        )
    }
}

impl Debug for KRISIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-R-ISIS[h={}, w={}, k={} q={}, length_bound={}, sigma={}, norm={}]",
            self.h, self.w, self.k, self.q, self.length_bound, self.sigma, self.norm
        )
    }
}

impl KRISIS{
    pub fn with_h(&self, h: usize) -> Self {
        KRISIS {
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
        KRISIS {
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
        KRISIS {
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
        KRISIS {
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
        SIS::new(
            self.h,
            self.q.clone(),
            self.length_bound,
            (self.w + 1) * self.h,
            self.norm,
        )
    }

    pub fn security_level(&self) -> f64 {
        self.to_sis().security_level()
    }

    pub fn security_level_internal(&self, estimate_type: Estimates) -> Result<f64, LatticeEstimatorError> {
        self.to_sis().security_level_internal(estimate_type)
    }

    pub fn upper_bound_h(&self) -> usize {
        self.to_sis().upper_bound_h()
    }
}


#[cfg(test)]

mod test{

    use crate::norms::Norm;

    use super::*;

    const Q: u128 = 2147483649;
    const SQRT_Q: f64 = 46340.95;

    //TODO find applicable tests

    #[test]
    fn test_rsis_security_level_l2() {
    }
}