use std::fmt;
use std::fmt::{Debug, Display};

use num_bigint::BigUint;

use crate::norms::Norm;
use crate::sis::SIS;
use crate::errors::LatticeEstimatorError;
use crate::reduction::Estimates;

pub struct PRISIS {
    pub h: usize,
    pub w: usize,
    pub d: usize,
    pub q: BigUint,
    pub length_bound: f64,
    pub sigma: f64,
    pub k: usize,
    pub norm: Norm,
}

impl Display for PRISIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "PRISIS[h={}, w={}, d={}, q={}, length_bound={}, sigma={}, k={}, norm={}]",
            self.h, self.w, self.d, self.q, self.length_bound, self.sigma, self.k, self.norm
        )
    }
}

impl Debug for PRISIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "PRISIS[h={}, w={}, d={}, q={}, length_bound={}, sigma={}, k={}, norm={}]",
            self.h, self.w, self.d, self.q, self.length_bound, self.sigma, self.k, self.norm
        )
    }
}

impl PRISIS {
    pub fn new(h: usize, w: usize, d: usize, q: BigUint, length_bound: f64, sigma: f64, k: usize, norm: Norm) -> Self {
        PRISIS {
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
        PRISIS {
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
        PRISIS {
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

    pub fn with_d(&self, d: usize) -> Self {
        PRISIS {
            h: self.h,
            w: self.w,
            d,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_length_bound(&self, length_bound: f64) -> Self {
        PRISIS {
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
        PRISIS {
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
        PRISIS {
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

    pub fn to_sis(&self) -> SIS {        
        
        //reduces to M-SIS
        SIS::new(
            self.h * self.d,
            self.w * self.h,
            self.q.clone(),
            self.length_bound,
            self.norm,
        )

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
            Ok(l) => l,
            Err(e) => panic!("Error: {:?}", e)
        }
    }
}