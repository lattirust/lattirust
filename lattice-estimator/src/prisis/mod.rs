use std::fmt;
use std::fmt::{Debug, Display};

use num_bigint::BigUint;
use num_traits::ToPrimitive;

use crate::norms::Norm;
use crate::sis::SIS;
use crate::errors::LatticeEstimatorError;
use crate::reduction::Estimates;

pub struct PRISIS {
    pub h: usize,
    pub k: usize,
    pub d: usize,
    pub q: BigUint,
    pub length_bound: f64,
    pub sigma: f64,
    pub w: usize,
    pub norm: Norm,
}

impl Display for PRISIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "PRISIS[h={}, w={}, k={}, d={}, q={}, length_bound={}, sigma={}, norm={}]",
            self.h, self.w, self.k, self.d, self.q, self.length_bound, self.sigma, self.norm
        )
    }
}

impl Debug for PRISIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "PRISIS[h={}, w={}, k={}, d={}, q={}, length_bound={}, sigma={}, norm={}]",
            self.h, self.w, self.k, self.d, self.q, self.length_bound, self.sigma, self.norm
        )
    }
}

impl PRISIS {
    pub fn with_h(&self, h: usize) -> Self {
        PRISIS {
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
        PRISIS {
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
        PRISIS {
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
        PRISIS {
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
        PRISIS {
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
        
        //reduces to M-SIS
        SIS::new(
            self.h * self.d,
            self.q.clone(),
            self.length_bound,
            self.w * self.h,
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
