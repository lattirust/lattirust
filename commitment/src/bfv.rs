#![allow(non_snake_case)]
#![allow(unused_imports)]

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

type PolyR<const M: u64, const N: usize> = Pow2CyclotomicPolyRing<Zq<M>, N>;
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

pub fn rand_tuple<const Q: u64, const N: usize>(factor: Option<PolyR<Q, N>>) ->  (PolyR<Q, N>, PolyR<Q, N>, PolyR<Q, N>) {
    let mut rng = rand::thread_rng();
    let size = WeightedTernaryChallengeSet::<PolyR<Q, N>>::byte_size();
    let f = factor.unwrap_or_else(|| PolyR::<Q, N>::one());

    (f * rand_ternary_poly(size, &mut rng), 
    f * get_gaussian(0.0, N, &mut rng), 
    f * get_gaussian(0.0, N, &mut rng))
}

pub struct PublicKey<const Q: u64, const P: u64, const N: usize> {
    pub poly1: PolyR<Q, N>, 
    pub poly2: PolyR<Q, N>,
    pub modulo: u64,
}
pub struct SecretKey<const Q: u64, const P: u64, const N: usize> {
    pub poly: PolyR<P, N>,  
    pub modulo: u64,
}

pub struct Plaintext<const P: u64, const N: usize> {
    pub poly: PolyR<P, N>,
    pub modulo: u64,
}

pub struct Ciphertext<const Q: u64, const N: usize> {
    pub c1: PolyR<Q, N>,
    pub c2: PolyR<Q, N>,
    pub modulo: u64,
}

impl<const Q: u64, const P: u64, const N: usize> SecretKey<Q, P, N> {
    pub fn new() -> Self {
        // generate random bytes to sample a ternary secret
        let size = WeightedTernaryChallengeSet::<PolyR<P, N>>::byte_size();

        Self {
            poly: rand_ternary_poly(size, &mut rand::thread_rng()),
            modulo: P,
        }
    }

    pub fn pk_gen(&self) -> PublicKey<Q, P, N> {
        let mut rng = rand::thread_rng();

        // e should have small std_dev (how small?) for the correctness, TODO: check the parameters
        // TODO: use OpenFHE DGS
        let e: PolyR<Q, N> = get_gaussian(0.0, N, &mut rng);

        // convert sk in PolyR to PolyR in order to perform the operation in PolyR
        let sk_zq: PolyR<Q, N> = convert_ring::<P, Q, N>(self.poly.clone());

        // compute the actual pk pair
        let poly2: PolyR<Q, N> = PolyR::rand(&mut rng);
        let poly1: PolyR<Q, N> = -(poly2.clone() * sk_zq + e);

        PublicKey {
            poly1,
            poly2,
            modulo: Q,
        }
    }
    
    pub fn decrypt(&self, c: Ciphertext<Q, N>) -> Plaintext<P, N> {
        let c1: PolyR<Q, N> = c.c1.clone();
        let c2: PolyR<Q, N> = c.c2.clone();
        let sk_zq: PolyR<Q, N> = convert_ring::<P, Q, N>(self.poly.clone());
        let raw: PolyR<Q, N> = c1 + c2 * sk_zq;

        let p = self.modulo as f64;
        let q = c.modulo as f64;
        let delta = p/q; 

        let coeffs: Vec<Zq<P>> = raw
            .coeffs()
            .into_iter()
            .map(|x| <Zq<P>>::from((delta * (SignedRepresentative::from(x).0 as f64)) as i128))
            .collect();
        let coeffs: [Zq<P>; N] = coeffs.try_into().expect("Bad format");

        Plaintext {
            poly: PolyR::from(coeffs),
            modulo: P,
        }
    }
}

impl<const Q: u64, const P: u64, const N: usize> PublicKey<Q, P, N> {
    pub fn encrypt(
        &self, 
        m: &Plaintext<P, N>, 
        r:  (PolyR<Q, N>, PolyR<Q, N>, PolyR<Q, N>)
        ) -> Ciphertext<Q, N> {
        let (pk1, pk2) = (self.poly1.clone(), self.poly2.clone());
        let (r0, r1, r2) = r;
        
        let p = m.modulo;
        let q = self.modulo;
        let delta = (q as f64 / p as f64).floor() as i128; // round off

        // mutliply each coeff of m with delta as bigint, and convert the polynomial into Zq, or convert m into PolyR and multiply, but attention overflow
        // TODO: define it as scalar - poly multiplication 
        let coeffs_zq: Vec<Zq<Q>> = m
            .poly
            .coeffs()
            .into_iter()
            .map(|x| <Zq<Q>>::from(delta * SignedRepresentative::from(x).0))
            .collect();
        let coeffs_zq: [Zq<Q>; N] = coeffs_zq.try_into().expect("Bad format");
        let m_delta: PolyR<Q, N> = PolyR::<Q, N>::from(coeffs_zq);
        // println!("m_delta: \n{m_delta:?}"); // OK
        
        // compute a, b
        let c1: PolyR<Q, N> = pk1 * r0.clone() + r1 + m_delta;
        let c2: PolyR<Q, N> = pk2 * r0.clone() + r2;

        // return the ciphertext
        Ciphertext {
            c1,
            c2,
            modulo: Q,
        }
    }
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

impl<const Q: u64, const N: usize>Default for Ciphertext<Q, N> {
    fn default() -> Self {
         Self {
            c1: PolyR::<Q, N>::zero(),
            c2: PolyR::<Q, N>::zero(),
            modulo: Q,
         }
    }
}