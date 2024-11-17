#![allow(non_snake_case)]

use ark_ff::{UniformRand, Zero};
use ark_std::rand;

use super::util::PolyR;
#[derive(
    Clone, Copy, Debug, Default
)]

pub struct Plaintext<const P: u64, const N: usize> {
    pub poly: PolyR<P, N>,
    pub modulo: u64,
}

impl<const P: u64, const N: usize> Plaintext<P, N> {
    pub fn rand_message() -> Self {
        Self {
            poly: PolyR::<P, N>::rand(&mut rand::thread_rng()),
            modulo: P,
        }
    }
    
    pub fn zero() -> Self {
        Self {
            poly: PolyR::<P, N>::zero(),
            modulo: P,
        }
    }
}
