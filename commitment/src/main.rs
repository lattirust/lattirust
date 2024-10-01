#[allow(dead_code)]
use std::ops::Neg;

use ark_ff::{One, UniformRand, Zero};
use ark_std::rand::prelude::SliceRandom;
use ark_std::rand;
// use commitments::ppk::get_gaussian_vec;
use lattirust_arithmetic::{challenge_set::{ternary, weighted_ternary::WeightedTernaryChallengeSet}, linear_algebra::{Matrix, Scalar, Vector}, ntt::ntt_modulus, ring::{ConvertibleRing, Pow2CyclotomicPolyRing, Zq}, traits::FromRandomBytes};
use rand::{CryptoRng, RngCore, SeedableRng};
use commitments::ppk;

const N: usize = 128;
const Q: u64 = ntt_modulus::<N>(16);
type R = Zq<Q>; 
type PolyR = Pow2CyclotomicPolyRing<R, N>;

fn main() {
    let rng = &mut rand::thread_rng();
    // // let t = 12;         // Plaintext modulus
    // // let q = 65536;      // Ciphertext modulus
    // let std_dev = 1.0;  // Standard deviation for generating the error
    // // let degree = 4;     // Degree of polynomials used for encoding and encrypting messages

    // let v = get_gaussian_vec(std_dev, N, &mut rng);
    // println!("{:?}", v);

    // // construct the polynomial
    // let tmp: Vec<R> = v.into_iter().map(|x| R::from(x)).collect();
    // let poly = PolyR::from(tmp);
    // println!("{:?}", poly);
    // let mut rng = rand::thread_rng();
    // let coeff: Matrix<R> = Matrix::<R>::rand_ternary(1, N, &mut rng);
    // let coeff: Vector<R> = Vector::from_fn(N, |_, _| {
    //         [-R::one(), R::zero(), R::one()]
    //             .choose(&mut rng)
    //             .unwrap()
    //             .clone()
    //     });
    
    // let coeff = coeff.iter().map(|x| R::from(x)).collect();
    
    let vec = Vector::<u8>::rand(N, rng);
    let vec_slice = vec.as_slice();
    let vec = WeightedTernaryChallengeSet::<PolyR>::try_from_random_bytes(vec_slice).unwrap();
    println!("{:?}", vec);
}