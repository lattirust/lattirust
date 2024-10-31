#![allow(non_snake_case)]

use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::linear_algebra::{Matrix, Vector};
use lattirust_arithmetic::ntt::ntt_modulus;

// use fhe::bfv::Ciphertext;
use lattirust_arithmetic::ring::{ConvertibleRing, PolyRing, Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT, SignedRepresentative, UnsignedRepresentative, Zq};
use rand::{CryptoRng, RngCore};
use rand_distr::Normal;
use rand::distributions::{Distribution};
use ark_ff::{One, UniformRand, Zero};
use ark_std::rand::prelude::SliceRandom;
use ark_std::rand;
use lattirust_arithmetic::challenge_set::ternary;
use lattirust_arithmetic::traits::FromRandomBytes;
use crate::bfv::{Ciphertext, Plaintext, PublicKey, SecretKey, rand_tuple};

pub struct Prover<const Q: u64, const P: u64, const N: usize> {
    // params: ParamsBFV,
    secret_key: SecretKey<Q, P, N>, // polynomial w/ coefficients in {-1, 0, +1}
    public_key: PublicKey<Q, P, N>, // ring mod q
    verifier_pk: PublicKey<Q, P, N>,
    message: Plaintext<P, N>, // ring mod p
}

pub struct Verifier<const Q: u64, const P: u64, const N: usize> {
    // params: ParamsBFV,
    secret_key: SecretKey<Q, P, N>,
    public_key: PublicKey<Q, P, N>,
    commitment: Option<Vec<Pow2CyclotomicPolyRing<Zq<Q>, N>>>
}

impl<const Q: u64, const P: u64, const N: usize> Prover<Q, P, N> {
    // instantiation
    fn new(verifier_pk: PublicKey<Q, P, N>, message: Plaintext<P, N>) -> Self {
        let secret_key =  SecretKey::new(); 
        let public_key = secret_key.pk_gen();
        
        Self {
            secret_key, 
            public_key, 
            verifier_pk,
            message,
        }
    }

    // generate-phase
    fn generate(&self) -> Ciphertext<Q, N> {
        let two = Pow2CyclotomicPolyRing::<Zq<Q>, N>::from(Zq::<Q>::from(2));
        
        // todo: use DGS
        let r = rand_tuple::<Q, P, N>(Some(two));
        let c = self.verifier_pk.encrypt(&self.message, r);
        
        c
    }

    // prove-phase, which is a sigma protocol with soundness amplification
    fn commit(&self, l: usize) -> Vec<Ciphertext<Q, N>> {
        let two = Pow2CyclotomicPolyRing::<Zq<Q>, N>::from(Zq::<Q>::from(2));

        let mut u: Vec<Pow2CyclotomicPolyRing<Zq<P>, N>> = Vec::new(); // vec poly of size l
        let mut y: Vec<(Pow2CyclotomicPolyRing<Zq<Q>, N>, Pow2CyclotomicPolyRing<Zq<Q>, N>, Pow2CyclotomicPolyRing<Zq<Q>, N>)> = Vec::new(); // vec of 3-tuple poly of size l
        let mut w: Vec<Ciphertext<Q, N>> = Vec::new(); // vec of ciphertexts

        for i in 0..l {
            u.push(Pow2CyclotomicPolyRing::<Zq<P>, N>::rand(&mut rand::thread_rng()));
            y.push(rand_tuple::<Q, P, N>(Some(two)));
            let ui = Plaintext::<P, N> {
                poly: u[i].clone(),
                modulo: P,
            };
            w.push(self.verifier_pk.encrypt(&ui, y[i].clone()))
        }

        w
    }

    // or a tuple of two arrays
    fn reveal(response: Vec<Pow2CyclotomicPolyRing<Zq<P>, N>>) -> Vec<(Pow2CyclotomicPolyRing<Zq<P>, N>, Pow2CyclotomicPolyRing<Zq<P>, N>)> {
        todo!();
    }
}

impl<const Q: u64, const P: u64, const N: usize> Verifier<Q, P, N> {
    // instantiation
    // fn new() -> Self {
    //     Self {
    //         secret_key: Pow2CyclotomicPolyRing<Rp, N>::new(10),
    //         public_key: secret_key.generate_pk(),
    //         commitment: None
    //     }
    // }

    // fn new(secret_key: SecretKey) -> Self {
    //     Self {
    //         secret_key,
    //         public_key: secret_key.generate_pk(),
    //         commitment: None
    //     }
    // }

    fn challenge(&self, commitment: Vec<Pow2CyclotomicPolyRing<Zq<P>, N>>) -> Vec<Pow2CyclotomicPolyRing<Zq<P>, N>> {
        todo!();
        // self.commitment = Some(commitment);
    }

    fn verify(response: Vec<(Pow2CyclotomicPolyRing<Zq<P>, N>, Pow2CyclotomicPolyRing<Zq<P>, N>)>) -> bool {
        todo!();
    }
}
