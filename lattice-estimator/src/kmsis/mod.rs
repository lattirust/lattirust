use std::fmt;
use std::fmt::{Debug, Display};

use num_bigint::BigUint;

use crate::norms::Norm;
use crate::sis::SIS;
use crate::errors::LatticeEstimatorError;
use crate::reduction::Estimates;

pub struct KMSIS {
    pub h: usize,
    pub k: usize,
    pub d: usize,
    pub q: BigUint,
    pub length_bound: f64,
    pub sigma: f64,
    pub w: usize,
    pub norm: Norm,
}

impl Display for KMSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-M-SIS[h={}, w={}, k={}, d={}, q={}, length_bound={}, sigma={}, norm={}]",
            self.h, self.w, self.k, self.d, self.q, self.length_bound, self.sigma, self.norm
        )
    }
}

impl Debug for KMSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-M-SIS[h={}, w={}, k={}, d={}, q={}, length_bound={}, sigma={}, norm={}]",
            self.h, self.w, self.k, self.d, self.q, self.length_bound, self.sigma, self.norm
        )
    }
}

impl KMSIS {
    pub fn with_h(&self, h: usize) -> Self {
        KMSIS {
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
        KMSIS {
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
        KMSIS {
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
        KMSIS {
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
        KMSIS {
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
            self.h,
            self.q.clone(),
            self.length_bound,
            self.w - self.k,
            self.norm,
        )

    }

    pub fn security_level(&self) -> f64 {
        self.to_sis().security_level()
    }

    pub fn security_level_internal(&self, estimate_type: Estimates) -> Result<f64, LatticeEstimatorError> {

        if self.k >= self.w {
            return Err(LatticeEstimatorError::InvalidParameter { param_name : self.k.to_string() , reason : "k should be smaller than w, otherwise reduction does not work".to_string() });
        }

        self.to_sis().security_level_internal(estimate_type)
    }
}

pub struct KRSIS {
    pub h: usize,
    pub k: usize,
    pub q: BigUint,
    pub length_bound: f64,
    pub sigma: f64,
    pub w: usize,
    pub norm: Norm,
}

impl Display for KRSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-R-SIS[h={}, w={}, k={}, q={}, length_bound={}, sigma={}, norm={}]",
            self.h, self.w, self.k, self.q, self.length_bound, self.sigma, self.norm
        )
    }
}

impl Debug for KRSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-R-SIS[h={}, w={}, k={}, q={}, length_bound={}, sigma={}, norm={}]",
            self.h, self.w, self.k, self.q, self.length_bound, self.sigma, self.norm
        )
    }
}


impl KRSIS {
    pub fn with_h(&self, h: usize) -> Self {
        KRSIS {
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
        KRSIS {
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
        KRSIS {
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
        KRSIS {
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
            self.w - self.k,
            self.norm,
        )

    }

    pub fn security_level(&self) -> f64 {
        self.to_sis().security_level()
    }

    pub fn security_level_internal(&self, estimate_type: Estimates) -> Result<f64, LatticeEstimatorError> {
        self.to_sis().security_level_internal(estimate_type)
    }
}