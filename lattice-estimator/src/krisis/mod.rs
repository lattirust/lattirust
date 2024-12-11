use std::fmt;
use std::fmt::{Debug, Display};

use num_bigint::BigUint;

use crate::errors::LatticeEstimatorError;
use crate::norms::Norm;
use crate::sis::SIS;
use crate::reduction::Estimates;

pub struct KMISIS {
    pub h: usize,
    pub w: usize,
    pub d: usize,
    pub q: BigUint,
    pub length_bound: f64,
    pub sigma: f64,
    pub k: usize,
    pub norm: Norm,
}

impl Display for KMISIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-M-ISIS[h={}, w={}, d={}, q={}, length_bound={}, sigma={}, k={}, norm={}]",
            self.h, self.w, self.d, self.q, self.length_bound, self.sigma, self.k, self.norm
        )
    }
}

impl Debug for KMISIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-M-ISIS[h={}, w={}, d={}, q={}, length_bound={}, sigma={}, k={}, norm={}]",
            self.h, self.w, self.d, self.q, self.length_bound, self.sigma, self.k, self.norm
        )
    }
}

impl KMISIS{

    pub fn new(h: usize, w: usize, d: usize, q: BigUint, length_bound: f64, sigma: f64, k: usize, norm: Norm) -> Self {
        KMISIS {
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
        KMISIS {
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
        KMISIS {
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

    pub fn with_length_bound(&self, length_bound: f64) -> Self {
        KMISIS {
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
        KMISIS {
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
        KMISIS {
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
        
        SIS::new(
            self.h * self.d,
            (self.w + 1) * self.h * self.d,
            self.q.clone(),
            self.length_bound,
            self.norm,
        )

    }

    pub fn upper_bound_h(&self) -> usize {
        self.to_sis().upper_bound_h().div_floor(self.d)
    }
    pub fn security_level_internal(&self, estimate_type: Estimates) -> Result<f64, LatticeEstimatorError> {
        self.to_sis().security_level_internal(estimate_type)
    }

    pub fn find_optimal_h_annealing(&self, lambda: usize, estimate_type: Estimates) -> usize {
        match self.to_sis().find_optimal_h_annealing(lambda, estimate_type) {
            Ok(h) => h.div_floor(self.d),
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


pub struct KRISIS {
    pub h: usize,
    pub w: usize,
    pub q: BigUint,
    pub length_bound: f64,
    pub sigma: f64,
    pub k: usize,
    pub norm: Norm,
}

impl Display for KRISIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-R-ISIS[h={}, w={}, q={}, length_bound={}, sigma={}, k={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.sigma, self.k, self.norm
        )
    }
}

impl Debug for KRISIS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "K-R-ISIS[h={}, w={}, q={}, length_bound={}, sigma={}, k={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.sigma, self.k, self.norm
        )
    }
}

impl KRISIS{

    pub fn new(h: usize, w: usize, q: BigUint, length_bound: f64, sigma: f64, k: usize, norm: Norm) -> Self {
        KRISIS {
            h,
            w,
            q,
            length_bound,
            sigma,
            k,
            norm,
        }
    }

    pub fn with_h(&self, h: usize) -> Self {
        KRISIS {
            h,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_w(&self, w: usize) -> Self {
        KRISIS {
            h: self.h,
            w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_length_bound(&self, length_bound: f64) -> Self {
        KRISIS {
            h: self.h,
            w: self.w,
            q: self.q.clone(),
            length_bound,
            sigma: self.sigma,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_sigma(&self, sigma: f64) -> Self {
        KRISIS {
            h: self.h,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_k(&self, k: usize) -> Self {
        KRISIS {
            h: self.h,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            sigma: self.sigma,
            k,
            norm: self.norm,
        }
    }

    pub fn to_sis(&self) -> SIS {
        SIS::new(
            self.h,
            (self.w + 1) * self.h,
            self.q.clone(),
            self.length_bound,
            self.norm,
        )
    }

    pub fn security_level_internal(&self, estimate_type: Estimates) -> Result<f64, LatticeEstimatorError> {
        self.to_sis().security_level_internal(estimate_type)
    }

    pub fn upper_bound_h(&self) -> usize {
        self.to_sis().upper_bound_h()
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


#[cfg(test)]

mod test{

    use crate::{krisis::KMISIS, norms::Norm, reduction::Estimates};

    
    #[test]
    fn test_sis_security_level_l2() {
        let falcon512_unf: KMISIS = KMISIS::new(24, 512,4, 122898899u64.into(), 5833.9072, 4.0, 2, Norm::L2);

        let cost_internal:Result<f64, crate::errors::LatticeEstimatorError> = falcon512_unf.security_level_internal(Estimates::Kyber(true));

        match cost_internal {
            Ok(cost) => {
                println!("Internal : {falcon512_unf} -> lambda: {cost}");
            },
            Err(e) => {
                panic!("Error: {:?}", e);
            }
            
        }
    
    }

    #[test]
    fn test_sis_security_level_linf() {

        let falcon512_unf: KMISIS = KMISIS::new(24, 512,4, 122898899u64.into(), 5833.9072, 4.0, 2, Norm::Linf);

        let cost_internal:Result<f64, crate::errors::LatticeEstimatorError> = falcon512_unf.security_level_internal(Estimates::Kyber(true));

        match cost_internal {
            Ok(cost) => {
                println!("Internal : {falcon512_unf} -> lambda: {cost}");
            },
            Err(e) => {
                panic!("Error: {:?}", e);
            }
            
        }
    
    }


    #[test]
    fn test_sis_h_annealing() {
        let falcon512_unf: KMISIS = KMISIS::new(24, 512,4, 122898899u64.into(), 5833.9072, 4.0, 2, Norm::L2);

        falcon512_unf.with_h(falcon512_unf.upper_bound_h());
        let h_opt = falcon512_unf.find_optimal_h_annealing(64, Estimates::Kyber(true));
        let lambda_opt = falcon512_unf.with_h(h_opt).security_level_internal(Estimates::Kyber(true)).unwrap();

        println!("Optimal h: {h_opt} -> lambda: {lambda_opt}");

        let test = falcon512_unf.with_h(29).security_level_internal(Estimates::Kyber(true)).unwrap();
        println!("Test h: 29 -> lambda: {test}");


        let test = falcon512_unf.with_h(30).security_level_internal(Estimates::Kyber(true)).unwrap();
        println!("Test h: 30 -> lambda: {test}");

        let test = falcon512_unf.with_h(31).security_level_internal(Estimates::Kyber(true)).unwrap();
        println!("Test h: 31 -> lambda: {test}");
    }

    #[test]
    fn test_sis_lenbound_annealing() {
        let falcon512_unf: KMISIS = KMISIS::new(24, 512,4, 122898899u64.into(), 5833.9072, 4.0, 2, Norm::L2);

        let len_opt = falcon512_unf.find_optimal_length_bound_annealing(64, Estimates::Kyber(true));
        let lambda_opt = falcon512_unf.with_length_bound(len_opt).security_level_internal(Estimates::Kyber(true)).unwrap();

        println!("Optimal length bound: {len_opt} -> lambda: {:?}", lambda_opt);

        let lambda_test = falcon512_unf.with_length_bound(40000.0).security_level_internal(Estimates::Kyber(true)).unwrap();
        println!("Test length bound: 40000 -> lambda: {:?}", lambda_test);
    }
}