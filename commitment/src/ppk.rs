#![allow(non_snake_case)]
#![allow(unused_imports)]

use std::default;

use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::linear_algebra::{Matrix, Vector};
use lattirust_arithmetic::ntt::ntt_modulus;

// use fhe::bfv::Ciphertext;
use lattirust_arithmetic::ring::{ConvertibleRing, PolyRing, Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT, SignedRepresentative, UnsignedRepresentative, Zq};
use rand::{CryptoRng, RngCore, Rng};
use rand_distr::num_traits::zero;
use rand_distr::Normal;
use rand::distributions::{Distribution};
use ark_ff::{One, UniformRand, Zero};
use ark_std::rand::prelude::SliceRandom;
use ark_std::rand;
use lattirust_arithmetic::challenge_set::ternary;
use lattirust_arithmetic::traits::FromRandomBytes;
use crate::bfv::{Ciphertext, Plaintext, PublicKey, SecretKey, rand_tuple, convert_ring};

type PolyR<const M: u64, const N: usize> = Pow2CyclotomicPolyRing<Zq<M>, N>;
type TuplePolyR<const Q: u64, const N: usize> = (PolyR<Q, N>, PolyR<Q, N>, PolyR<Q, N>);

pub struct Prover<const Q: u64, const P: u64, const N: usize> {
    // params: ParamsBFV,
    secret_key: SecretKey<Q, P, N>, // polynomial w/ coefficients in {-1, 0, +1}
    public_key: PublicKey<Q, P, N>, // ring mod q
    message: Plaintext<P, N>, // ring mod p
    r: TuplePolyR<Q, N>,
    u: Vec<PolyR<P, N>>,
    y: Vec<TuplePolyR<Q, N>>,
    l: usize,
}

pub struct Verifier<const Q: u64, const P: u64, const N: usize> {
    // params: ParamsBFV,
    secret_key: SecretKey<Q, P, N>,
    public_key: PublicKey<Q, P, N>,
    prover_pk: PublicKey<Q, P, N>,
    l: usize,
    c: Ciphertext<Q, N>, 
    w: Vec<Ciphertext<Q, N>>,
    gamma: Vec<PolyR<P, N>>,
}

impl<const Q: u64, const P: u64, const N: usize> Prover<Q, P, N> {
    // instantiation
    // TODO: cannot use verifier pk to encrypt, would leak m
    pub fn new(message: Plaintext<P, N>, l: usize) -> Self {
        let secret_key =  SecretKey::new(); 
        let public_key = secret_key.pk_gen();
        
        Self {
            secret_key, 
            public_key, 
            message,
            // TODO: create defaults
            r: (PolyR::<Q, N>::zero(), PolyR::<Q, N>::zero(), PolyR::<Q, N>::zero()),
            u: Vec::default(),
            y: Vec::default(),
            l
        }
    }

    // TODO: move PublicKey out and derive clone
    pub fn return_pk(&self) -> PublicKey<Q, P, N> {
        PublicKey {
            poly1: self.public_key.poly1.clone(),
            poly2: self.public_key.poly2.clone(),
            modulo: Q,
        }
    }

    // generate-phase
    // TODO: use RefCell and make self immutable (?)
    pub fn generate(&mut self) -> Ciphertext<Q, N> {
        let two = PolyR::<Q, N>::from(Zq::<Q>::from(2));

        // todo: use DGS
        let r = rand_tuple::<Q, N>(None);
        self.r = r.clone();
        let r = (r.0.clone() * two, r.1.clone() * two, r.2.clone() * two);

        self.public_key.encrypt(&self.message, r)
    }

    // prove-phase, which is a sigma protocol with soundness amplification
    pub fn commit(&mut self, l: usize) -> Vec<Ciphertext<Q, N>> {
        let two = PolyR::<Q, N>::from(Zq::<Q>::from(2));
        let mut rng = rand::thread_rng();
        let mut u: Vec<PolyR<P, N>> = Vec::new(); // vec poly of size l
        let mut y_times_2: Vec<TuplePolyR<Q, N>> = Vec::new(); 
        let mut w: Vec<Ciphertext<Q, N>> = Vec::new(); // vec of ciphertexts

        for i in 0..l {
            u.push(PolyR::<P, N>::rand(&mut rng));

            let yi: TuplePolyR<Q, N> = rand_tuple::<Q, N>(None);
            self.y.push(yi);
            y_times_2.push((yi.0 * two, yi.1 * two, yi.2 * two));
            
            // transform ui into a Plaintext so we can encrypt it
            let ui = Plaintext::<P, N> {
                poly: u[i].clone(),
                modulo: P,
            };
            // TODO: create a helper method to use self.encrypt instead of self.pk.encrypt
            w.push(self.public_key.encrypt(&ui, y_times_2[i].clone())) 
        }

        self.u = u;
        
        w
    }

    pub fn reveal(&self, challenge: Vec<PolyR<P, N>>) -> (Vec<PolyR<P, N>>, Vec<TuplePolyR<Q, N>>) {
        let mut v: Vec<PolyR<P, N>> = Vec::new();
        let mut z: Vec<TuplePolyR<Q, N>> = Vec::new();
        let l = self.l;
        let m = self.message.poly.clone();
        let r = self.r.clone();

        for i in 0..l {
            let u_i = self.u[i].clone();
            let y_i = self.y[i].clone();
            let gamma_i = challenge[i].clone();

            let v_i: PolyR<P, N> = u_i + gamma_i * m;

            let gamma_i: PolyR<Q, N> = convert_ring::<P, Q, N>(gamma_i);
            let z_i = (
                y_i.0.clone() + gamma_i.clone() * r.0.clone(), 
                y_i.1.clone() + gamma_i.clone() * r.1.clone(),
                y_i.2.clone() + gamma_i.clone() * r.2.clone()
            );

            v.push(v_i);
            z.push(z_i);
        }


        (v, z)
    }
}

impl<const Q: u64, const P: u64, const N: usize> Verifier<Q, P, N> {
    // instantiation
    pub fn new(prover_pk: PublicKey<Q, P, N>, l: usize) -> Self {
        let secret_key =  SecretKey::new(); 
        let public_key = secret_key.pk_gen();

        Self {
            secret_key, 
            public_key, 
            prover_pk,
            l, 
            c: Ciphertext::default(),
            w: Vec::default(),
            gamma: Vec::default(),
        }
    }
    
    // instantiation with no prover_pk and l, 
    // pub fn new() -> Self {
    //     let secret_key =  SecretKey::new(); 
    //     let public_key = secret_key.pk_gen();

    //     Self {
    //         secret_key, 
    //         public_key, 
    //         prover_pk,
    //         l, 
    //         c: Ciphertext::default(),
    //         w: Vec::default(),
    //         gamma: Vec::default(),
    //     }
    // }

    pub fn end_genenerate(&mut self, c: Ciphertext<Q, N>) {
        assert_ne!(self.l, 0); // not sure if it's necessary
        self.c = c;
    }

    pub fn challenge(&mut self, w: Vec<Ciphertext<Q, N>>, l: usize) -> Vec<PolyR<P, N>> {
        self.w = w;

        // sample l monomials
        let mut rng = rand::thread_rng();
        let mut gamma = Vec::new();

        for _ in 0..l {
            let mut coeffs: [Zq<P>; N] = [Zq::<P>::zero(); N];
            let rand_index: usize = rng.gen_range(0..N);
            coeffs[rand_index] = Zq::<P>::one();
            gamma.push(PolyR::from(coeffs));
        }
        self.gamma = gamma.clone();

        gamma
    }

    pub fn verify(&self, opening: (Vec<PolyR<P, N>>, Vec<TuplePolyR<Q, N>>)) -> bool {
        let two = PolyR::<Q, N>::from(Zq::<Q>::from(2));
        let (v, z) = (opening.0, opening.1);
        let (w, c, gamma) = (&self.w, &self.c, &self.gamma.clone());
        let pk= &self.prover_pk;

        for i in 0..self.l {
            let vi = Plaintext{ poly: v[i].clone(), modulo: P };
            let ctxt = pk.encrypt(&vi, (two * z[i].0.clone(), two * z[i].1.clone(), two * z[i].2.clone()));
            let left = (ctxt.c1, ctxt.c2);
            let gamma_zq = convert_ring::<P, Q, N>(gamma[i].clone());
            let right = (w[i].c1 + gamma_zq * c.c1.clone(), w[i].c2 + gamma_zq * c.c2.clone());

            if left != right {
                return false;
            }

            // TODO: check the second condition
        }

        return true;
    }
}
