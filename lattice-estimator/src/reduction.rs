use num_traits::{Float, Pow};

use crate::errors::LatticeEstimatorError as EstError;
use core::panic;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::f64::consts::{PI, E};
use std::thread::panicking;

//parameter holder for LLL
#[derive(Copy, Clone)]
pub struct LLL{
    d: usize,
    entry_bitsize: usize
}

//cost estimate for LLL reduction, only one possibility
impl LLL {
    fn cost(&self) -> Result<f64, EstError> {
        if self.entry_bitsize == 0 {
            Ok(f64::powi(self.d as f64, 3))
        } else {
            Ok(f64::powi(self.d as f64, 3) * f64::powi(self.entry_bitsize as f64, 2))
        }
    }
}

struct ShortVectorCost(f64, f64, usize, usize);

//parameter holder for BKZ_like reduction
#[derive(Copy, Clone)]
pub struct BKZ{
    block_size: usize,
    d: usize,
    entry_bitsize: usize
}

//Enum of available cost estimates for lattice reduction
pub enum Estimates{
    CheNgue12,
    BDGL16,
    LaaMosPol14,
    ABFKSW20,
    ABLR21,
    ADPS16,
    ChaLoy21,
    Kyber(bool)
}

pub enum EstimatesParams{
    CheNgue12(BKZ),
    BDGL16(BKZ),
    LaaMosPol14(BKZ),
    ABFKSW20(BKZ),
    ABLR21(BKZ),
    ADPS16(BKZ),
    ChaLoy21(BKZ),
    Kyber(BKZ, bool),
}

pub fn cost(reduction_type: EstimatesParams) -> Result<f64, EstError> {
    match reduction_type {
        EstimatesParams::CheNgue12(bkz) => chengue12_cost(bkz),
        EstimatesParams::BDGL16(bkz) => bdgl16_cost(bkz),
        EstimatesParams::LaaMosPol14(bkz) => laamospol14_cost(bkz),
        EstimatesParams::ABFKSW20(bkz) => abfksw20_cost(bkz),
        EstimatesParams::ABLR21(bkz) => ablr21_cost(bkz),
        EstimatesParams::ADPS16(bkz) => adbs16_cost(bkz),
        EstimatesParams::ChaLoy21(bkz) => chaloy21_cost(bkz),
        EstimatesParams::Kyber(bkz, classical) => kyber_cost(bkz, classical),
    }
}

//LLL builder
impl LLL {
    pub fn new(d: usize, entry_bitsize: usize) -> Self {
        LLL { d, entry_bitsize }
    }
}

impl BKZ {
    pub fn new(block_size: usize, d: usize, entry_bitsize: usize) -> Self {
        BKZ { block_size, d, entry_bitsize }
    }
}

fn chengue12_cost(bkz: BKZ)-> Result<f64, EstError> {
    let svp_nb: usize;
    if bkz.block_size < bkz.d {
        svp_nb = 8 * bkz.d
    } else {
        svp_nb = 1
    }
    let cost: f64 = 0.270188776350190 as f64 * bkz.block_size as f64 * (bkz.block_size as f64).ln()
                    - 1.0192050451318417 * bkz.block_size as f64 + 16.10253135200765   + (100.0 as f64).log2();


    
    let lll_cost: f64 = LLL::new(bkz.d, bkz.entry_bitsize).cost().unwrap();
    let new_cost = lll_cost + (svp_nb as f64) * (2.0 as f64).powf(cost);
    Ok(new_cost)
}

//This is only the classical cost estimate, the quantum would be 0.2650
fn adbs16_cost(bkz: BKZ) -> Result<f64, EstError> {
    Ok((2.0 as f64).powf(0.292 * bkz.block_size as f64))
}

fn chaloy21_cost(bkz: BKZ) -> Result<f64, EstError> {
    Ok((2.0 as f64).powf(0.2570 * bkz.block_size as f64))
}

fn bdgl16_cost(bkz: BKZ) -> Result<f64, EstError> {
    let tours: usize;
    if bkz.block_size < bkz.d {
        tours = 8 * bkz.d
    } else {
        tours = 1
    }
    
    let lll_cost: f64 = LLL::new(bkz.d, bkz.entry_bitsize).cost().unwrap();
    if bkz.block_size < 90 {
        let new_cost = lll_cost + (2.0 as f64).powf(0.387 * bkz.block_size as f64 + 16.4 + (tours as f64).log2());
        Ok(new_cost)
    } else {
        let new_cost = lll_cost + (2.0 as f64).powf(0.292 * bkz.block_size as f64 + 16.4 + (tours as f64).log2());
        Ok(new_cost)
    }
}

//todo find the difference
fn laamospol14_cost(bkz: BKZ) -> Result<f64, EstError> {
    let tours: usize;
    if bkz.block_size < bkz.d {
        tours = 8 * bkz.d
    } else {
        tours = 1
    }
    let lll_cost: f64 = LLL::new(bkz.d, bkz.entry_bitsize).cost().unwrap();
    let new_cost = lll_cost + (2.0 as f64).powf(0.265 * bkz.block_size as f64 + 16.4 + (tours as f64).log2());
    Ok(new_cost)
}

fn abfksw20_cost(bkz: BKZ) -> Result<f64, EstError> {
    let tours: usize;
    if bkz.block_size < bkz.d {
        tours = 8 * bkz.d
    } else {
        tours = 1
    }
    let lll_cost: f64 = LLL::new(bkz.d, bkz.entry_bitsize).cost().unwrap();
    let power:f64;
    if bkz.block_size <= 92 {
        power = 0.1839 * bkz.block_size as f64 * (bkz.block_size as f64).log2() - 0.995 * bkz.block_size as f64 + 16.25 + (64.0 as f64).log2();
    } else {
        power = 0.125 * bkz.block_size as f64 * (bkz.block_size as f64).log2() - 0.547 * bkz.block_size as f64 + 10.4 + (64.0 as f64).log2();
    }
    Ok(tours as f64 * (2.0 as f64).powf(power) + lll_cost)
}

fn ablr21_cost(bkz: BKZ) -> Result<f64, EstError> {
    let tours: usize;
    if bkz.block_size < bkz.d {
        tours = 8 * bkz.d
    } else {
        tours = 1
    }
    let lll_cost: f64 = LLL::new(bkz.d, bkz.entry_bitsize).cost().unwrap();
    let power:f64;
    if bkz.block_size <= 97 {
        power = 0.1839 * bkz.block_size as f64 * (bkz.block_size as f64).log2() - 1.077 * bkz.block_size as f64 + 29.12 + (64.0 as f64).log2();
    } else {
        power = 0.125 * bkz.block_size as f64 * (bkz.block_size as f64).log2() - 0.654 * bkz.block_size as f64 + 25.84 + (64.0 as f64).log2();
    }
    Ok(tours as f64 * (2.0 as f64).powf(power) + lll_cost)
}


fn matzov_cost(bkz: BKZ) {
    unimplemented!()
}

fn gj21_cost(bkz: BKZ) -> Result<f64, EstError> {
    unimplemented!()
}

//classical vs quantum decoding
fn kyber_cost(bkz: BKZ, classical: bool) -> Result<f64, EstError> {
    let t: (f64, f64);
    if classical {
        t = (0.2988026130564745, 26.011121212891872);
    } else {
        t = (0.26944796385592995, 28.97237346443934);
    }
    let sieve_free_dim: f64 = f64::max(bkz.block_size as f64 * (4.0/3.0 as f64).log2() / (bkz.block_size as f64 / (2.0 * PI * E)).log2(), 0.0);
    
    if bkz.block_size < 20 {
        return chengue12_cost(bkz);
    } else {
        let svp_calls: f64 = 5.46 * max(bkz.d- bkz.block_size , 1) as f64;
        let beta = bkz.block_size as f64 - sieve_free_dim;
        let gates: f64 = 5.46 * (2.0 as f64).powf(t.0 * beta + t.1);
        let lll_cost: f64 = LLL::new(bkz.d, bkz.entry_bitsize).cost().unwrap();
        Ok(lll_cost + svp_calls * gates)
    }
}

fn kyber_short_vectors(bkz: BKZ, nb_vec_out: Option<usize>, sampling_cost: Option<usize>, classical: bool) -> ShortVectorCost {
    let new_block_size: usize = (bkz.block_size as f64 - f64::max(bkz.block_size as f64 * (4.0/3.0 as f64).log2() / (bkz.block_size as f64 / (2.0 * PI * E)).log2(), 0.0)).floor() as usize;

    match nb_vec_out{
        None => {
            let new_nb_vec_out: usize = ((2.0).powf(0.2075 *  new_block_size as f64)).floor() as usize;
            let cost:f64 = new_nb_vec_out as f64 / (2.0).powf(0.2075 *  new_block_size as f64).floor();
            return ShortVectorCost(1.1547, cost.ceil() * kyber_cost(bkz, classical), cost.ceil() * ((2.0).powf(0.2075 *  new_block_size as f64)).floor(), new_block_size);
        }
        Some(nb) => {
            if nb == 1 {
                return ShortVectorCost(1.0, kyber_cost(bkz, true), bkz.block_size, 1);
            } else {
                let cost:f64 = nb as f64 / (2.0).powf(0.2075 *  new_block_size as f64).floor();
                return ShortVectorCost(1.1547, cost.ceil() * kyber_cost(bkz, classical), cost.ceil() * ((2.0).powf(0.2075 *  new_block_size as f64)).floor(), new_block_size);
            }
        }
    }

}

fn gj21_short_vectors(bkz: BKZ, nb_vec_out: Option<usize>, sampling_cost: Option<usize>, classical: bool) -> ShortVectorCost {
    let new_block_size: usize = (bkz.block_size as f64 - f64::max(bkz.block_size as f64 * (4.0/3.0 as f64).log2() / (bkz.block_size as f64 / (2.0 * PI * E)).log2(), 0.0)).floor() as usize;
    let mut sieve_dim: usize = new_block_size;
    
    if bkz.block_size < bkz.d {
        sieve_dim = min(bkz.d, (new_block_size + ((bkz.d - bkz.block_size) as f64 * 5.46 as f64).log2()).floor() as usize);
    }

    let scaling_fact_rho: f64 = (4.0/3.0).sqrt() * bkz_delta(sieve_dim).powf(sieve_dim - 1 as f64) * bkz_delta(bkz.block_size).powf(1 - sieve_dim);
    let mut c : Option<f64> = None;
    match nb_vec_out{
        None => {
            let c = Some(1.0);
        }
        Some(nb) => {
            if nb == 1 {
                return ShortVectorCost(1.0, kyber_cost(bkz, true), bkz.block_size, 1);
            } else {
                let cost:f64 = nb as f64 / (2.0).powf(0.2075 *  new_block_size as f64).floor();
                return ShortVectorCost(1.1547, cost.ceil() * kyber_cost(bkz, classical), cost.ceil() * ((2.0).powf(0.2075 *  new_block_size as f64)).floor(), new_block_size);
            }
        }
    }


}

pub fn bkz_delta(block_size: usize) -> f64 {

    let small_approximations = [
        (2, 1.02190),
        (5, 1.01862),
        (10, 1.01616),
        (15, 1.01485),
        (20, 1.01420),
        (25, 1.01342),
        (28, 1.01331),
        (40, 1.01295),
    ];

    // Collect the small approximations into a HashMap
    let approx_map: HashMap<usize, f64> = small_approximations.iter().cloned().collect();
    
    // Case for block_size <= 2
    if block_size <= 2 {
        return *approx_map.get(&2).unwrap();
    } 
    // Case for 2 < block_size < 40: find the closest smaller value
    else if block_size < 40 {
        // Sort the keys
        let mut keys: Vec<usize> = approx_map.keys().cloned().collect();
        keys.sort();
        
        // Find the largest key smaller than or equal to block_size
        let mut closest = 2;
        for &key in &keys {
            if key > block_size {
                break;
            }
            closest = key;
        }
        return *approx_map.get(&closest).unwrap();
    } 
    // Case for block_size == 40
    else if block_size == 40 {
        return *approx_map.get(&40).unwrap();
    } 
    // Case for block_size > 40: apply the formula
    else {
        return (block_size as f64 / (2.0 * PI * E) * (PI * block_size as f64).powf(1.0 / block_size as f64))
            .powf(1.0 / (2.0 * (block_size as f64 - 1.0)));
    }

}

//approx is the Kannan factor for the GSA simulator
pub fn gsa(d: usize, n: usize, q: usize, block_size: usize, approx: Option<f64>) -> Vec<f64> {
    let log_volume: f64 = match approx {
        None => (q as f64).log2() * (d - n) as f64 + (1.0 as f64).log2() * n as f64,
        Some(approx) => (q as f64).log2() * (d - n - 1) as f64 + (1.0 as f64).log2() * n as f64 + approx.log2(),
    };

    let delta = bkz_delta(block_size);
    let mut r: Vec<f64> = Vec::with_capacity(d);
    for i in 0..d {
        r.push((2.0_f64).powf((2.0 * (d - 1 - 2 * i) as f64) * delta.log2() + log_volume / d as f64));
    }
    r
}


