#![allow(non_snake_case)]
#![allow(unused_imports)]

use std::default;
use std::ops::{Add, Mul};

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::linear_algebra::{Matrix, Vector};
use lattirust_arithmetic::ntt::ntt_modulus;

use lattirust_arithmetic::ring::{ConvertibleRing, PolyRing, Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT, SignedRepresentative, UnsignedRepresentative, Zq};
use rand::{CryptoRng, RngCore};
use rand_distr::Normal;
use rand::distributions::Distribution;
use ark_ff::{One, UniformRand, Zero};
use ark_std::rand;
use lattirust_arithmetic::traits::FromRandomBytes;

pub type PolyR<const M: u64, const N: usize> = Pow2CyclotomicPolyRing<Zq<M>, N>;

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct TuplePolyR<const Q: u64, const N: usize>(pub PolyR<Q, N>, pub PolyR<Q, N>, pub PolyR<Q, N>);

impl<const Q: u64, const N: usize> Default for TuplePolyR<Q, N> {
    fn default() -> Self {
        TuplePolyR(PolyR::zero(), PolyR::zero(), PolyR::zero())
    }
}

impl<const Q: u64, const N: usize> Zero for TuplePolyR<Q, N> {
    fn zero() -> Self {
        TuplePolyR::<Q, N>::default()
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero() && self.1.is_zero() && self.2.is_zero()
    }
}

impl<const Q: u64, const N: usize> Add for TuplePolyR<Q, N> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        TuplePolyR(self.0 + rhs.0, self.1 + rhs.1, self.2 + rhs.2)
    }
}

impl<const Q: u64, const N: usize> Mul<PolyR<Q, N>> for TuplePolyR<Q, N>  {
    type Output = Self;
    fn mul(self, rhs: PolyR<Q, N>) -> Self::Output {
        TuplePolyR(self.0 * rhs, self.1 * rhs, self.2 * rhs)
    }
}
// TODO: use the same rng everywhere

// TODO: toy sampling, need to use OpenFHE code 
pub fn get_gaussian_vec<Rng: RngCore + CryptoRng, const Q: u64>(  
    std_dev: f64, 
    dimension: usize, 
    rng: &mut Rng
) -> Vec<Zq<Q>> {
    let gaussian = Normal::new(0.0, std_dev).unwrap();
    let val: Vec<Zq<Q>> = (0..dimension)
        .map(|_| Zq::<Q>::from(gaussian.sample(rng) as i64))
        .collect();

    val
}

pub fn get_gaussian<Rng: RngCore + CryptoRng, const Q: u64, const N: usize>(std_dev: f64, dimension: usize, rng: &mut Rng) -> PolyR<Q, N> {
    let rand_vec: Vec<Zq<Q>> = get_gaussian_vec(std_dev, dimension, rng);
    let rand_vec: [Zq<Q>; N] = rand_vec.try_into().expect("Bad format");
    
    PolyR::<Q, N>::from(rand_vec)
}

pub fn convert_ring<const SOURCE: u64, const TARGET: u64, const N: usize>(poly: PolyR<SOURCE, N>) -> PolyR<TARGET, N> {
    let coeffs: Vec<Zq<SOURCE>> = poly.coeffs();
    let coeffs: Vec<Zq<TARGET>> = coeffs.into_iter()
        .map(|x| <Zq<TARGET>>::from(SignedRepresentative::from(x).0))
        .collect();
    let coeffs: [Zq<TARGET>; N] = coeffs.try_into().expect("Bad format");

    PolyR::<TARGET, N>::from(coeffs)
}

pub fn rand_ternary_poly<Rng: RngCore + CryptoRng, const P: u64, const N: usize>(size: usize, rng: &mut Rng) -> PolyR<P, N> {
    let bytes: Vector<u8> = Vector::<u8>::rand(size, rng);
    WeightedTernaryChallengeSet::<PolyR::<P, N>>::try_from_random_bytes(bytes.as_slice()).unwrap()
}

pub fn rand_tuple<const Q: u64, const N: usize>(factor: Option<PolyR<Q, N>>) ->  
TuplePolyR<Q, N> {
    let mut rng = rand::thread_rng();
    let size = WeightedTernaryChallengeSet::<PolyR<Q, N>>::byte_size();
    let f = factor.unwrap_or_else(|| PolyR::<Q, N>::one());

    TuplePolyR(f * rand_ternary_poly(size, &mut rng), 
    f * get_gaussian(0.0, N, &mut rng), 
    f * get_gaussian(0.0, N, &mut rng))
}
