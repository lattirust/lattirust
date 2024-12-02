use std::fmt;
use std::fmt::{Debug, Display};

use num_bigint::BigUint;
use crate::errors::LatticeEstimatorError;
use crate::norms::Norm;
use crate::sis::SIS;
use crate::reduction::Estimates;

pub struct BASISrand {
    h: usize,
    q: BigUint,
    length_bound: f64,
    w: usize,
    sigma: usize,
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
    pub fn with_h(&self, h: usize) -> Self {
        BASISrand {
            h,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            k: self.k,
            norm: self.norm
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
            norm: self.norm
        }
    }

    pub fn to_sis(&self) -> SIS {
        SIS::new(
            self.h,
            self.q.clone(),
            self.length_bound,
            self.w,
            self.norm,
        )
    }

    pub fn security_level(&self) -> f64 {
        self.to_sis().security_level()
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
}



#[cfg(test)]

mod test{
    use crate::norms::Norm;

    use super::*;

    const Q: u128 = 2147483649;
    const SQRT_Q: f64 = 46340.95;

    //TODO find applicable tests

    //TODO BASISstruct relies on k-M-SIS, which to my knowledge does not have a reduction yet
    #[test]
    fn test_rsis_security_level_l2() {
        panic!("No tests implemented for BASISrand yet");
    }
}