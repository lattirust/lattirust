#![allow(non_snake_case)]

use ark_ff::Zero;
use super::util::PolyR;

pub struct Ciphertext<const Q: u64, const N: usize> {
    pub c1: PolyR<Q, N>,
    pub c2: PolyR<Q, N>,
    pub modulo: u64,
}

impl<const Q: u64, const N: usize>Default for Ciphertext<Q, N> {
    fn default() -> Self {
         Self {
            c1: PolyR::<Q, N>::zero(),
            c2: PolyR::<Q, N>::zero(),
            modulo: Q,
         }
    }
}