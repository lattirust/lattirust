use num_traits::Float;
use crate::errors::LatticeEstimatorError as EstError;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::f64::consts::{PI, E};

//Cost estimates for LLL (heuristic cost)
fn lll_cost(lattice_dim: usize, bit_size: Option<usize>) -> f64 {

    match bit_size {
        None => f64::powi(lattice_dim as f64, 3),
        Some(bit_size) => f64::powi(lattice_dim as f64, 3) * f64::powi(bit_size as f64, 2)
    }
}

//parameter holder for BKZ_like reduction
#[derive(Copy, Clone)]
pub struct BKZ{
    block_size: usize,
    d: usize,
    entry_bitsize: Option<usize>
}

impl BKZ {
    pub fn new(block_size: usize, d: usize, entry_bitsize: Option<usize>) -> Self {
        BKZ { block_size, d, entry_bitsize }
    }
}

//Enum of available cost estimates for lattice reduction
#[derive(Debug, Clone, Copy)]
pub enum Estimates{
    CheNgue12,
    BDGL16,
    LaaMosPol14,
    ABFKSW20,
    ABLR21,
    ADPS16,
    ChaLoy21,
    Kyber(bool),
    Matzov(bool)
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
    Matzov(BKZ, bool)
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
        EstimatesParams::Matzov(bkz, classical) => matzov_cost(bkz, classical)
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


    
    let lll_cost: f64 = lll_cost(bkz.d, bkz.entry_bitsize);
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
    
    let lll_cost: f64 = lll_cost(bkz.d, bkz.entry_bitsize);
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
    let lll_cost: f64 = lll_cost(bkz.d, bkz.entry_bitsize);
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
    let lll_cost: f64 = lll_cost(bkz.d, bkz.entry_bitsize);
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
    let lll_cost: f64 = lll_cost(bkz.d, bkz.entry_bitsize);
    let power:f64;
    if bkz.block_size <= 97 {
        power = 0.1839 * bkz.block_size as f64 * (bkz.block_size as f64).log2() - 1.077 * bkz.block_size as f64 + 29.12 + (64.0 as f64).log2();
    } else {
        power = 0.125 * bkz.block_size as f64 * (bkz.block_size as f64).log2() - 0.654 * bkz.block_size as f64 + 25.84 + (64.0 as f64).log2();
    }
    Ok(tours as f64 * (2.0 as f64).powf(power) + lll_cost)
}

fn reduce_dimension(block_size: usize) -> f64{
    f64::max(block_size as f64 * f64::ln(4.0 / 3.0) / (block_size as f64 / (2.0 * PI * E)).ln(), 0.0)
}

//classical vs quantum decoding
fn kyber_cost(bkz: BKZ, classical: bool) -> Result<f64, EstError> {
    let t: (f64, f64);
    if classical {
        t = (0.2988026130564745, 26.011121212891872);
    } else {
        t = (0.26944796385592995, 28.97237346443934);
    }

    if bkz.block_size < 20 {
        return chengue12_cost(bkz);
    } else {
        let svp_calls: f64 = 5.46 * max(bkz.d - bkz.block_size , 1) as f64;
        let beta = bkz.block_size as f64 - reduce_dimension(bkz.block_size);
        let gates: f64 = 5.46 * (2.0 as f64).powf(t.0 * beta + t.1);
        let lll_cost: f64 = lll_cost(bkz.d, bkz.entry_bitsize);
        Ok(lll_cost + svp_calls * gates)
    }
}

fn matzov_cost(bkz: BKZ, classical: bool) -> Result<f64, EstError> {
    let t: (f64, f64);
    if classical {
        t = (0.29613500308205365, 20.387885985467914);
    } else {
        t = (0.2663676536352464, 25.299541499216627);
    }

    if bkz.block_size < 20 {
        return chengue12_cost(bkz);
    } else {
        let svp_calls: f64 = 5.46 * max(bkz.d - bkz.block_size , 1) as f64;
        let beta = bkz.block_size as f64 - reduce_dimension(bkz.block_size);
        let gates: f64 = 5.46 * (2.0 as f64).powf(t.0 * beta + t.1);
        let lll_cost: f64 = lll_cost(bkz.d, bkz.entry_bitsize);
        Ok(lll_cost + svp_calls * gates)
    }
}

fn kyber_short_vectors(bkz: BKZ, nb_vec_out: Option<usize>, classical: bool) -> (f64, f64, usize, usize) {

    let beta_: f64 = bkz.block_size as f64 - reduce_dimension(bkz.block_size).floor() as f64;

    let nb_vec_out: usize = match nb_vec_out {
        Some(1) => return (1.0, kyber_cost(bkz, classical).unwrap(), bkz.block_size, 1),
        Some(n) => n,
        None => return(1.1547,
        kyber_cost(bkz, classical).unwrap(),
        ((2.0_f64).powf(0.2075 * beta_)).floor() as usize,
        beta_ as usize,
        )
    };
  
    let c = nb_vec_out as f64 / (2.0_f64).powf(0.2075 * beta_);
    (
        1.1547,
        c.ceil() * kyber_cost(bkz, classical).unwrap(),
        (c.ceil() * (2.0_f64).powf(0.2075 * beta_)).floor() as usize,
        beta_ as usize,
    )
}

//As a convention let's use usize_max as a flag for the cost being impossible to compute
pub fn matzov_short_vectors(bkz: BKZ, nb_vec_out: Option<usize>, classical: bool) -> (f64, f64, usize, usize) {
    
    let _beta: usize = bkz.block_size - reduce_dimension(bkz.block_size).floor() as usize;

    let t: (f64, f64);
    if classical {
        t = (0.29613500308205365, 20.387885985467914);
    } else {
        t = (0.2663676536352464, 25.299541499216627);
    }

    let sieve_dim: usize;
    if bkz.block_size < bkz.d {
        sieve_dim = min(bkz.d, (_beta as f64 + ((bkz.d - bkz.block_size) as f64 * 5.46 as f64).log2() / t.0).floor() as usize);
    } else {
        sieve_dim = _beta;
    }

    let scaling_fact_rho: f64 = (4.0/3.0).sqrt() * bkz_delta(sieve_dim).powf(sieve_dim as f64 - 1.0) * bkz_delta(bkz.block_size).powf(1.0 - sieve_dim as f64);
    
    let new_nb_vec_out: usize = match nb_vec_out {
        None => ((2.0).powf(0.2075 *  sieve_dim as f64)).floor() as usize,
        Some(1) => return (1.0, kyber_cost(bkz, true).unwrap(), bkz.block_size, 1), 
        Some(_) => nb_vec_out.unwrap()
    };

    let c: f64 = new_nb_vec_out as f64 / (2.0).powf(0.2075 *  sieve_dim as f64).floor();
    let sieve_cost: f64 = 5.46 * (2.0).powf(t.0 *  sieve_dim as f64 + t.1);
    
    if c > (2.0).powf(10000.0 as f64) {
        return (scaling_fact_rho, f64::INFINITY, usize::MAX, sieve_dim);
    } else {
        let final_cost = c.ceil() * (matzov_cost(bkz, classical).unwrap() + sieve_cost);
        return (
            scaling_fact_rho,
            final_cost,
            (c.ceil() * (2.0).powf(0.2075 *  sieve_dim as f64).floor()) as usize,
            sieve_dim);
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
pub fn gsa_simulator(d: usize, n: usize, q: usize, block_size: usize, approx: Option<f64>) -> Vec<f64> {
    let log_volume: f64 = match approx {
        None => (q as f64).log2() * (d - n) as f64 + (1.0_f64).log2() * n as f64,
        Some(approx) => (q as f64).log2() * (d - n - 1) as f64 + (1.0_f64).log2() * n as f64 + approx.log2(),
    };

    let delta = bkz_delta(block_size);
    let mut r_log: Vec<f64> = Vec::with_capacity(d);

    // First loop to calculate `next` values and store them in `r_log`
    for i in 0..d {
        let next: f64 = (d as isize - 1 - 2 * i as isize) as f64 * delta.log2() + log_volume / d as f64;
        r_log.push(next);
    }

    // Second loop to apply the power transformation
    for i in 0..d {
        r_log[i] = (2.0_f64).powf(2.0 * r_log[i]);
    }

    r_log
}


#[cfg(test)]
mod test {
    use statrs::assert_almost_eq;
    use crate::{norms::Norm, sis::SIS};

    use super::*;

    #[test]
    fn test_chengue12_cost() {
        let bkz: BKZ = BKZ::new(500, 1024, None);
        let cost = chengue12_cost(bkz).unwrap().log2();
        assert_eq!(cost.round(), 366.0);
    }

    #[test]
    fn test_abfskw20_cost() {
        let bkz: BKZ = BKZ::new(500, 1024, None);
        let cost = abfksw20_cost(bkz).unwrap().log2();
        assert_eq!(cost.round(), 316.0);
    }

    #[test]
    fn test_ablr21_cost() {
        let bkz: BKZ = BKZ::new(500, 1024, None);
        let cost = ablr21_cost(bkz).unwrap().log2();
        assert_eq!(cost.round(), 278.0);
    }

    #[test]
    fn test_adps16_cost() {
        let bkz: BKZ = BKZ::new(500, 1024, None);
        let cost = adbs16_cost(bkz).unwrap().log2();
        assert_eq!(cost.round(), 146.0);
    }

    #[test]
    fn test_kyber_cost() {
        let bkz: BKZ = BKZ::new(500, 1024, None);
        let cost = kyber_cost(bkz, true).unwrap().log2();
        assert_eq!(cost.round(), 177.0);
    }

    
    #[test]
    fn test_kyber_short_vectors(){
        let bkz: BKZ = BKZ::new(100, 500, None);
        let sv = kyber_short_vectors(bkz, None, true);
        assert_almost_eq!(sv.1, 2.736747612813679e19, 1e-2);
        assert_almost_eq!(sv.0, 1.1547, 1e-2);
        assert_eq!(sv.2, 176584);
        assert_eq!(sv.3, 84);
    }

    #[test]
    fn test_kyber_short_vectors_2(){
        let bkz: BKZ = BKZ::new(100, 500, None);
        let sv: (f64, f64, usize, usize) = kyber_short_vectors(bkz, Some(1000), true);
        assert_almost_eq!(sv.1, 2.736747612813679e19, 1e-2);
        assert_almost_eq!(sv.0, 1.1547, 1e-2);
        assert_eq!(sv.2, 176584);
        assert_eq!(sv.3, 84);
    }

    #[test]
    fn test_matzov_short_vectors() {
        let bkz: BKZ = BKZ::new(100, 500, None);
        let sv = matzov_short_vectors(bkz, None, true);
        assert_almost_eq!(sv.1, 9.33915764560094e17, 1e3);
        assert_almost_eq!(sv.0, 1.04228014727497, 1e-2);
        assert_eq!(sv.2, 36150192);
        assert_eq!(sv.3, 121);
    }

    #[test]
    fn test_by_default_l2_param() {
        let falcon512_unf: SIS = SIS::new(512, 12289u64.into(), 5833.9072, 1024, Norm::L2);
        let lambda = falcon512_unf.security_level();
        assert!(lambda >= 128.);
        println!("External : {falcon512_unf} -> lambda: {lambda}");

        let cost_internal: f64 = falcon512_unf.security_level_internal(Estimates::Matzov(true));
        println!("Internal : {falcon512_unf} -> lambda: {cost_internal}");
    }

    #[test]
    fn simulator_cn11(){
        let n = 128;
        let d = 213;
        let q = 2048;
        let beta = 40;

        let vec: Vec<f64> = gsa_simulator(d, n, q, beta, None);
        let mut sum = 0.0;
        for i in 0..vec.len() {
            sum += vec[i].ln();
        }
        assert_almost_eq!(sum, 1296.18522764710, 1e-2);
    }
}
