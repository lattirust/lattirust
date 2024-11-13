#![allow(non_snake_case)]

use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::ring::{PolyRing, Pow2CyclotomicPolyRing, SignedRepresentative, Zq};
use ark_ff::UniformRand;
use ark_std::rand;
use lattirust_arithmetic::traits::FromRandomBytes;

use super::ciphertext::Ciphertext;
use super::plaintext::Plaintext;
use super::public_key::PublicKey;
use super::util::{convert_ring, get_gaussian, rand_ternary_poly};

type PolyR<const M: u64, const N: usize> = Pow2CyclotomicPolyRing<Zq<M>, N>;
// TODO: use the same rng everywhere

pub struct SecretKey<const Q: u64, const P: u64, const N: usize> {
    pub poly: PolyR<P, N>,  
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

        // TODO: try to devide by q then multiply by p
        // or the opposite
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
