

use std::ops::Mul;
use crate::lattice_arithmetic::matrix::Vector;

use crate::lattice_arithmetic::ring::Ring;
use crate::lattice_arithmetic::traits::{IntegerDiv, Modulus, Normed, WithConjugationAutomorphism, WithLog2};

pub trait PolyRing: Ring
+ Mul<Self::BaseRing, Output=Self>
+ Normed<u64>
+ From<Vec<Self::BaseRing>>
+ WithConjugationAutomorphism
+ From<Self::BaseRing>
{
    type BaseRing: Ring + IntegerDiv + WithLog2 + Modulus;
    fn coeffs(&self) -> Vec<Self::BaseRing>;
    fn flattened(vec: &Vector<Self>) -> Vector<Self::BaseRing> {
        Self::flattened_coeffs(vec).into()
    }
    fn flattened_coeffs(vec: &Vector<Self>) -> Vec<Self::BaseRing> {
        vec.as_slice().into_iter().flat_map(|x| x.coeffs()).collect::<Vec<Self::BaseRing>>()
    }
}
