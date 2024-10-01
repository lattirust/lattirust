pub(crate) struct LLL{
    d: usize,
    entry_bitsize: usize
}

impl LLL {
    pub fn new(d: usize, entry_bitsize: usize) -> Self {
        LLL { d, entry_bitsize }
    }
}

pub trait ReductionCost {
    fn cost(&self) -> f64;
}

impl ReductionCost for LLL {
    fn cost(&self) -> f64 {
        if self.entry_bitsize == 0 {
            f64::powi(self.d as f64, 3)
        } else {
            f64::powi(self.d as f64, 3) * f64::powi(self.entry_bitsize as f64, 2)
        }
    }
}

