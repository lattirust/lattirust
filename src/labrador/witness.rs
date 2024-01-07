#![allow(non_snake_case)]

use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::poly_ring::PolyRing;

pub struct Witness<R: PolyRing> {
    pub s: Vec<Vector<R>>,
}

impl<R: PolyRing> Witness<R> {
    pub fn new_dummy(rank: usize, multiplicity: usize, norm_bound: u64) -> Self {
        Self {
            // Absolutely not secure, but good enough for testing
            s: vec![Vector::<R>::from_fn(multiplicity, |_, _| R::from(1u128)); rank],
        }
    }
}