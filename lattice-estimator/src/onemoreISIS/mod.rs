use std::fmt;
use std::fmt::{Debug, Display};

use num_bigint::BigUint;

use crate::norms::Norm;
use crate::sis::SIS;
use crate::errors::LatticeEstimatorError;
use crate::reduction::Estimates;

pub struct onemoreISIS {
    pub h: usize,
    pub w: usize,
    pub q: BigUint,
    pub length_bound: f64,
    pub sigma: f64,
    pub norm: Norm,
}

impl Display for onemoreISIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "onemoreISIS[h={}, w={}, q={}, length_bound={}, sigma={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.sigma, self.norm
        )
    }
}

impl Debug for onemoreISIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "onemoreISIS[h={}, w={}, q={}, length_bound={}, sigma={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.sigma, self.norm
        )
    }
}
impl onemoreISIS {
    pub fn new(h: usize, w: usize, q: BigUint, length_bound: f64, sigma: f64, norm: Norm) -> Self {
        onemoreISIS {
            h,
            w,
            q,
            length_bound,
            sigma,
            norm,
        }
    }

    pub fn with_h(&self, h: usize) -> Self {
        onemoreISIS {
            h,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            norm: self.norm,
        }
    }

    pub fn with_w(&self, w: usize) -> Self {
        onemoreISIS {
            h: self.h,
            w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            norm: self.norm,
        }
    }

    pub fn with_length_bound(&self, length_bound: f64) -> Self {
        onemoreISIS {
            h: self.h,
            w: self.w,
            q: self.q.clone(),
            length_bound,
            sigma: self.sigma,
            norm: self.norm,
        }
    }

    pub fn with_sigma(&self, sigma: f64) -> Self {
        onemoreISIS {
            h: self.h,
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
            self.w,
            self.q.clone(),
            self.length_bound,
            self.norm,
        )
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
            Ok(l) => l,
            Err(e) => panic!("Error: {:?}", e)
        }
    }
}
