use std::fmt;
use std::fmt::{Debug, Display};

use num_bigint::BigUint;
use crate::errors::LatticeEstimatorError;
use crate::norms::Norm;
use crate::sis::SIS;
use crate::reduction::Estimates;

pub struct ISISf {
    h: usize,
    w: usize,
    q: BigUint,
    length_bound: f64,
    k: usize,
    norm: Norm,
}

impl Display for ISISf {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "ISISf[h={}, w={}, q={}, length_bound={}, k={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.k, self.norm
        )
    }
}

impl Debug for ISISf {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "ISISf[h={}, w={}, q={}, length_bound={}, k={}, norm={}]",
            self.h, self.w, self.q, self.length_bound, self.k, self.norm
        )
    }
}

impl ISISf {
    pub fn new(h: usize, w: usize, q: BigUint, length_bound: f64, k: usize, norm: Norm) -> Self {
        ISISf {
            h,
            w,
            q,
            length_bound,
            k,
            norm,
        }
    }

    pub fn with_h(&self, h: usize) -> Self {
        ISISf {
            h,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_w(&self, w: usize) -> Self {
        ISISf {
            h: self.h,
            w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_length_bound(&self, length_bound: f64) -> Self {
        ISISf {
            h: self.h,
            w: self.w,
            q: self.q.clone(),
            length_bound,
            k: self.k,
            norm: self.norm,
        }
    }

    pub fn with_k(&self, k: usize) -> Self {
        ISISf {
            h: self.h,
            w: self.w,
            q: self.q.clone(),
            length_bound: self.length_bound,
            k,
            norm: self.norm,
        }
    }

    pub fn to_sis(&self) -> SIS {

        //The reduction can be to SIS/ISIS if f modelled with ROM
        //The reduction can be to M-SIS/M-ISIS if we do an attack regardless of f
        //We take the ROM road here

        SIS::new(
            self.h,
            self.w - self.k,
            self.q.clone(),
            self.length_bound, 
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


    pub fn find_optimal_h_annealing(&self, lambda: usize, estimate_type: Estimates) -> usize {
        match self.to_sis().find_optimal_h_annealing(lambda, estimate_type) {
            Ok(h) => h,
            Err(e) => panic!("Error: {:?}", e)
        }
    }

    pub fn find_optimal_length_bound_annealing(&self, lambda: usize, estimate_type: Estimates) -> f64 {
        match self.to_sis().find_optimal_length_bound_annealing(lambda, estimate_type) {
            Ok(length_bound) => length_bound,
            Err(e) => panic!("Error: {:?}", e)
        }
    }
}



#[cfg(test)]

mod test{
    use crate::isisf::ISISf;
    use crate::norms::Norm;
    use crate::reduction::Estimates;

    
    #[test]
    fn test_basis_security_level_l2() {
        let falcon512_unf: ISISf = ISISf::new(512, 1024, 12289u64.into(), 5833.9072, 1, Norm::L2);

        let security_level = falcon512_unf.security_level_internal(Estimates::Matzov(true)).unwrap();
        print!("Security level: {}", security_level);

    }

    #[test]
    fn test_find_optimal_h(){
        let falcon512_unf: ISISf = ISISf::new(512, 1024, 12289u64.into(), 5833.9072, 1, Norm::L2);

        let h_optimal = falcon512_unf.find_optimal_h_annealing(122, Estimates::Matzov(true));
        println!("Optimal h: {}", h_optimal);

        let security_level = falcon512_unf.with_h(h_optimal).security_level_internal(Estimates::Matzov(true)).unwrap();
        println!("Security level: {}", security_level);
    }

    #[test]
    fn test_find_optimal_lenbound(){

        let falcon512_unf: ISISf = ISISf::new(512, 1024, 12289u64.into(), 5833.9072, 1, Norm::L2);
        let lenbound_optimal = falcon512_unf.find_optimal_length_bound_annealing(122, Estimates::Matzov(true));
        println!("Optimal length bound: {}", lenbound_optimal);

        let security_level = falcon512_unf.with_length_bound(lenbound_optimal).security_level_internal(Estimates::Matzov(true)).unwrap();
        println!("Security level: {}", security_level);
    }
}