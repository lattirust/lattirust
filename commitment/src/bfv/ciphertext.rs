#![allow(non_snake_case)]

use ark_ff::Zero;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use std::ops::Add;

use super::util::PolyR;
#[derive(
    Clone, Copy, Debug, Default, CanonicalSerialize, CanonicalDeserialize
)]
pub struct Ciphertext<const Q: u64, const N: usize> {
    pub c1: PolyR<Q, N>,
    pub c2: PolyR<Q, N>,
    pub modulo: u64,
}

impl<const Q: u64, const N: usize>Zero for Ciphertext<Q, N> {
    fn zero() -> Self {
         Self {
            c1: PolyR::<Q, N>::zero(),
            c2: PolyR::<Q, N>::zero(),
            modulo: Q,
         }
    }

    fn is_zero(&self) -> bool {
        self.c1.is_zero() && self.c2.is_zero() 
    }
}

impl<const Q: u64, const N: usize> Add for Ciphertext<Q, N> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            c1: self.c1 + other.c1,
            c2: self.c2 + other.c2,
            modulo: Q,
        }
    }
}