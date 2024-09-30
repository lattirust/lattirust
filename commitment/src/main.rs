#[allow(dead_code)]
use std::ops::Neg;

use ark_ff::{One, UniformRand, Zero};
use ark_std::rand::prelude::SliceRandom;
use ark_std::rand;
use commitments::ppk::get_gaussian_vec;
use lattirust_arithmetic::{challenge_set::ternary, linear_algebra::{Matrix, Scalar}, ntt::ntt_modulus, ring::{ConvertibleRing, Pow2CyclotomicPolyRingNTT, Zq}};
use rand::{CryptoRng, RngCore, SeedableRng};
use commitments::ppk::secret_ternary;

const N: usize = 128;
const Q: u64 = ntt_modulus::<N>(16);
type R = Zq<Q>; 
type PolyR = Pow2CyclotomicPolyRingNTT<Q, N>;

fn main() {
    let mut rng = rand::thread_rng();
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
    let mat: Vec<R> = secret_ternary(&mut rng);
    println!("{:?}", mat);
}