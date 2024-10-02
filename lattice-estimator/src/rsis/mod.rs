use std::fmt;
use std::fmt::{Debug, Display};

use num_bigint::BigUint;

use crate::norms::Norm;
use crate::sis::SIS;
use crate::reduction::Estimates;


/// MSIS parameters for instances $A \in R\_q^{\texttt{h}\times\texttt{w}}$ where $R\_q = \mathbb{Z}\_\texttt{q}\[X\]/(X^\texttt{d}+1)$ such that $A s = 0$ for some $s \in R\_\texttt{q}^\texttt{w}$ with ${\lVert s \rVert\}_\texttt{norm} \leq \texttt{length\\_bound}$.
pub struct RSIS {
    pub h: usize,
    pub q: BigUint,
    pub length_bound: f64,
    pub w: usize,
    pub norm: Norm,
}

impl Display for RSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "RSIS[h={}, w={}, q={}, length_bound={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.norm
        )
    }
}

impl Debug for RSIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "RSIS[h={}, w={}, q={}, length_bound={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.norm
        )
    }
}

impl RSIS {
    pub fn with_h(&self, h: usize) -> Self {
        RSIS {
            h,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            norm: self.norm,
        }
    }

    pub fn with_length_bound(&self, length_bound: f64) -> Self {
        RSIS {
            h: self.h,
            w: self.w,
            q: self.q.clone(),
            length_bound,
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

    /// Return $\lambda$ such that `MSIS\[h, w, d, q, length_bound\]` is $2^\lambda$-hard (for a given norm).
    /// We estimate the security by reducing to `SIS\[h\*d, w\*d, q, length_bound\]` and calling the SIS security estimator.
    pub fn security_level(&self) -> f64 {
        self.to_sis().security_level()
    }

    pub fn security_level_internal(&self, estimate_type: Estimates) -> f64 {
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


    #[test]
    fn test_rsis_security_level_l2() {
        let test_l2: RSIS = RSIS {
            h: 5,
            q: Q.into(),
            length_bound: SQRT_Q,
            w: 512,
            norm: Norm::L2,
        };
        let lambda = test_l2.security_level();
        println!("External {test_l2} -> lambda: {lambda}");

        let lambda = test_l2.security_level_internal(Estimates::LLL);
        println!("Internal  LLL {test_l2} -> lambda: {lambda}");

        let lambda = test_l2.security_level_internal(Estimates::CheNgue12);
        println!("Internal  CheNgue12 {test_l2} -> lambda: {lambda}");
    }
}