use std::{fmt, usize};
use std::fmt::{Debug, Display};

use num_bigint::BigUint;

use crate::errors::LatticeEstimatorError;
use crate::norms::Norm;
use crate::sis::SIS;
use crate::reduction::Estimates;

pub struct MSIS {
    pub h: usize,
    pub w: usize,
    pub d: usize,
    pub q: BigUint,
    pub length_bound: f64,
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
    pub fn new(h: usize, w: usize, d: usize, q: BigUint, length_bound: f64, norm: Norm) -> Self {
        MSIS {
            h,
            w,
            d,
            q,
            length_bound,
            norm,
        }
    }

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

    pub fn with_w(&self, w: usize) -> Self {
        MSIS {
            h: self.h,
            w,
            d: self.d,
            q: self.q.clone(),
            length_bound: self.length_bound,
            norm: self.norm,
        }
    }

    pub fn with_d(&self, d: usize) -> Self {
        MSIS {
            h: self.h,
            w: self.w,
            d,
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
            self.w * self.h,
            self.q.clone(),
            self.length_bound,
            self.norm,
        )
    }

    pub fn security_level_internal(&self, estimate_type: Estimates) -> Result<f64, LatticeEstimatorError> {
        self.to_sis().security_level_internal(estimate_type)
    }

    pub fn upper_bound_h(&self) -> usize {
        self.to_sis().upper_bound_h().div_floor(self.d)
    }

    pub fn find_optimal_h_annealing(&self, lambda: usize, estimate_type: Estimates) -> usize {
        match self.to_sis().find_optimal_h_annealing(lambda, estimate_type) {
            Ok(h) => h.div_floor(self.d),
            Err(e) => panic!("Error: {:?}", e),
        }
    }

    pub fn find_optimal_length_bound_annealing(&self, lambda: usize, estimate_type: Estimates) -> f64 {
        match self.to_sis().find_optimal_length_bound_annealing(lambda, estimate_type) {
            Ok(length_bound) => length_bound,
            Err(e) => panic!("Error: {:?}", e),
        }
    }

    pub fn simultaneous_optimize(&self, lambda: usize, estimate_type: Estimates) ->(f64, usize) {
        match self.to_sis().simultaneous_optimize(lambda, estimate_type){
            Ok((len_b, h)) => (len_b, h.div_floor(self.d)),
            Err(e) => panic!("Error: {:?}", e),
        }
    }


}

#[cfg(test)]
mod test {

    use crate::msis::MSIS;
    use crate::norms::Norm;
    use crate::reduction::Estimates;

    #[test]
    fn test_sis_security_level_l2() {
        let falcon512_unf: MSIS = MSIS::new(24, 512, 4, 122898899u64.into(), 5833.9072, Norm::L2);
        let lambda = falcon512_unf.security_level_internal(Estimates::Matzov(true)).unwrap();
        println!("External : {falcon512_unf} -> lambda: {lambda}");

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
        let falcon512_unf: MSIS = MSIS::new(24, 512, 4, 122898899u64.into(), 5833.9072, Norm::Linf);


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
        let falcon512_unf: MSIS = MSIS::new(24, 512, 4, 122898899u64.into(), 5833.9072, Norm::L2);

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
        let falcon512_unf: MSIS = MSIS::new(24, 512, 4, 122898899u64.into(), 5833.9072, Norm::L2);

        let len_opt = falcon512_unf.find_optimal_length_bound_annealing(64, Estimates::Kyber(true));
        let lambda_opt = falcon512_unf.with_length_bound(len_opt).security_level_internal(Estimates::Kyber(true)).unwrap();

        println!("Optimal length bound: {len_opt} -> lambda: {:?}", lambda_opt);

        let lambda_test = falcon512_unf.with_length_bound(40000.0).security_level_internal(Estimates::Kyber(true)).unwrap();
        println!("Test length bound: 40000 -> lambda: {:?}", lambda_test);
    }


    #[test]
    fn test_simultaneous_pareto(){
        let falcon512_unf: MSIS = MSIS::new(24, 512, 4, 122898899u64.into(), 5833.9072, Norm::L2);


        let sol: (f64, usize) = falcon512_unf.simultaneous_optimize(64, Estimates::Kyber(true));

        println!("Simultaneous optimization: {:?}", sol);
        let lambda_test = falcon512_unf.with_h(sol.1).with_length_bound(sol.0).security_level_internal(Estimates::Kyber(true)).unwrap();
        println!("Test h: 31, length bound: 75556.41376748933 -> lambda: {:?}", lambda_test);
    }
}
