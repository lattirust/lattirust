#![allow(non_snake_case)]

use ark_ff::Zero;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use std::ops::Add;
use super::util::PolyR;
use lattirust_arithmetic::nimue::serialization::{FromBytes, ToBytes};
#[derive(
    Clone, Copy, Debug, Default, CanonicalSerialize, CanonicalDeserialize,
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

// impl<const Q: u64, const N: usize> FromBytes for Ciphertext<Q, N> {
//     type FromBytesError = SerializationError;
//     fn from_bytes(buffer: &[u8]) -> Result<Self, Self::FromBytesError> {
//         // 1. compute the sizes of the members of Ciphertext
//         // 2. slice the bytes according to the sizes
//         // 3. do from_bytes one by one
//         let size_poly = PolyR::<Q, N>::default().to_bytes().unwrap().len();

//         let c1 = PolyR::<Q, N>::from_bytes(&buffer[..size_poly])?;
//         let c2 = PolyR::<Q, N>::from_bytes(&buffer[size_poly..2*size_poly])?;
//         let modulo = u64::from_bytes(&buffer[2*size_poly..])?;

//         Ok(Self { c1, c2, modulo })
//     }
// }