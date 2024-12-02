use std::fmt;
use std::fmt::{Debug, Display};

use num_bigint::BigUint;

use crate::errors::LatticeEstimatorError;
use crate::norms::Norm;
use crate::sis::SIS;
use crate::reduction::Estimates;

pub struct MSIS {
    pub h: usize,
    pub d: usize,
    pub q: BigUint,
    pub length_bound: f64,
    pub w: usize,
    pub norm: Norm,
}

impl Display for MSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "MSIS[h={}, w={}, d={}, q={}, length_bound={}, norm={}]",
            self.h, self.w, self.d, self.q, self.length_bound, self.norm
        )
    }
}

impl Debug for MSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "MSIS[h={}, w={}, d={}, q={}, length_bound={}, norm={}]",
            self.h, self.w, self.d, self.q, self.length_bound, self.norm
        )
    }
}

impl MSIS {
    pub fn with_h(&self, h: usize) -> Self {
        MSIS {
            h,
            w: self.w,
            d: self.d,
            q: self.q.clone(),
            length_bound: self.length_bound,
            norm: self.norm,
        }
    }

    pub fn with_length_bound(&self, length_bound: f64) -> Self {
        MSIS {
            h: self.h,
            w: self.w,
            d: self.d,
            q: self.q.clone(),
            length_bound,
            norm: self.norm,
        }
    }


    //should it be dxd or w*h
    pub fn to_sis(&self) -> SIS {
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

    pub fn upper_bound_h(&self) -> usize {
        self.to_sis().upper_bound_h().div_floor(self.d)
    }
}
