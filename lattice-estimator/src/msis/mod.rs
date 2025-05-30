use std::fmt;
use std::fmt::{Debug, Display};

use num_bigint::BigUint;
use num_traits::ToPrimitive;

use crate::norms::Norm;
use crate::sis::SIS;

pub mod lattice_estimator;
pub mod security_estimates;

/// MSIS parameters for instances $A \in R\_q^{\texttt{h}\times\texttt{w}}$ where $R\_q = \mathbb{Z}\_\texttt{q}\[X\]/(X^\texttt{d}+1)$ such that $A s = 0$ for some $s \in R\_\texttt{q}^\texttt{w}$ with ${\lVert s \rVert\}_\texttt{norm} \leq \texttt{length\\_bound}$.
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

    pub fn to_sis(&self) -> SIS {
        SIS::new(
            self.h * self.d,
            self.q.clone(),
            self.length_bound,
            self.w * self.d,
            self.norm,
        )
    }

    /// Return $\lambda$ such that `MSIS\[h, w, d, q, length_bound\]` is $2^\lambda$-hard (for a given norm).
    /// We estimate the security by reducing to `SIS\[h\*d, w\*d, q, length_bound\]` and calling the SIS security estimator.
    pub fn security_level(&self) -> f64 {
        self.to_sis().security_level()
    }

    pub fn upper_bound_h(&self) -> usize {
        self.to_sis().upper_bound_h().div_floor(self.d)
    }
}

const ROOT_HERMITE_FACTOR: f64 = 1.0045;

pub fn msis_h_128_l2(msis: &MSIS) -> Option<usize> {
    let log_q = msis.q.to_f64().unwrap().log2();

    if msis.length_bound >= msis.q.to_f64().unwrap() / 2. {
        return None;
    }

    let min_h = (f64::log2(2. * msis.length_bound) / 2.).powi(2)
        / (msis.d as f64 * log_q * f64::log2(ROOT_HERMITE_FACTOR));
    Some(min_h.ceil() as usize)
}

pub fn msis_h_128_linf(msis: &MSIS) -> Option<usize> {
    let msis_l2 = MSIS {
        h: msis.h,
        d: msis.d,
        q: msis.q.clone(),
        length_bound: msis.length_bound * ((msis.d * msis.w) as f64).sqrt(),
        w: msis.w,
        norm: Norm::L2,
    };

    msis_h_128_l2(&msis_l2)
}
