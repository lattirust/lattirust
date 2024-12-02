use std::fmt;
use std::fmt::{Debug, Display};

use num_bigint::BigUint;

use crate::errors::LatticeEstimatorError;
use crate::norms::Norm;
use crate::sis::SIS;
use crate::reduction::Estimates;

pub struct VSIS {
    h: usize,
    d: usize,
    q: BigUint,
    length_bound: f64,
    w: usize,
    norm: Norm,
}

impl Display for VSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "RSIS[h={}, w={}, q={}, length_bound={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.norm
        )
    }
}

impl Debug for VSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "RSIS[h={}, w={}, q={}, length_bound={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.norm
        )
    }
}

impl VSIS {
    pub fn with_h(&self, h: usize) -> Self {
        VSIS {
            h,
            d: self.d,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            norm: self.norm,
        }
    }

    pub fn with_length_bound(&self, length_bound: f64) -> Self {
        VSIS {
            h: self.h,
            d: self.d,
            w: self.w,
            q: self.q.clone(),
            length_bound,
            norm: self.norm,
        }
    }

    pub fn with_d(&self, d: usize) -> Self {
        VSIS {
            h: self.h,
            d,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            norm: self.norm,
        }
    }

    pub fn to_sis(&self) -> SIS {
        SIS::new(
            self.h,
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
        let test_l2: VSIS = VSIS {
            h: 50,
            d: 50,
            q: Q.into(),
            length_bound: SQRT_Q,
            w: 512,
            norm: Norm::L2,
        };
        let lambda = test_l2.security_level();
        println!("External {test_l2} -> lambda: {lambda}");

        let lambda = test_l2.security_level_internal(Estimates::CheNgueEnum).unwrap();
        println!("Internal  CheNgue12 {test_l2} -> lambda: {lambda}");
    }
}