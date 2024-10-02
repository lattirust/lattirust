use crate::errors::LatticeEstimatorError as E;

//parameter holder for LLL
#[derive(Copy, Clone)]
pub struct LLL{
    d: usize,
    entry_bitsize: usize
}

//cost estimate for LLL reduction, only one possibility
impl LLL {
    fn cost(&self) -> Result<f64, E> {
        if self.entry_bitsize == 0 {
            Ok(f64::powi(self.d as f64, 3))
        } else {
            Ok(f64::powi(self.d as f64, 3) * f64::powi(self.entry_bitsize as f64, 2))
        }
    }
}

//parameter holder for BKZ_like reduction
#[derive(Copy, Clone)]
pub struct BKZ{
    block_size: usize,
    d: usize,
    entry_bitsize: usize
}

//Enum of available cost estimates for lattice reduction
pub enum Estimates{
    LLL,
    CheNgue12
}

pub enum EstimatesParams{
    LLL(LLL),
    CheNgue12(BKZ)
}

pub fn cost(reduction_type: EstimatesParams) -> Result<f64, E> {
    match reduction_type {
        EstimatesParams::LLL(lll) => lll.cost(),
        EstimatesParams::CheNgue12(bkz) => chengue12_cost(bkz)
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

fn chengue12_cost(bkz: BKZ)-> Result<f64, E> {
    let svp_nb: usize;
    if bkz.block_size < bkz.d {
        svp_nb = 8 * bkz.d
    } else {
        svp_nb = 1
    }
    let cost: f64 = 0.270188776350190 as f64 * bkz.block_size as f64 * (bkz.block_size as f64).ln()
                    - 1.0192050451318417 * bkz.block_size as f64 + 16.10253135200765   + (100.0 as f64).log2();


    
    let lll_cost = LLL::new(bkz.d, bkz.entry_bitsize).cost().unwrap();
    let new_cost = lll_cost + (svp_nb as f64) * (2.0 as f64).powf(cost);
    Ok(new_cost)
}

// TO DO A Complete Analysis of the BKZ Lattice
