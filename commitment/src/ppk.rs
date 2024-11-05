#![allow(non_snake_case)]
#![allow(unused_imports)]

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
    verifier_pk: PublicKey<Q, P, N>,
    message: Plaintext<P, N>, // ring mod p
    r: TuplePolyR<Q, N>,
    u: Vec<PolyR<P, N>>,
    y: Vec<TuplePolyR<Q, N>>,
}

pub struct Verifier<const Q: u64, const P: u64, const N: usize> {
    // params: ParamsBFV,
    secret_key: SecretKey<Q, P, N>,
    public_key: PublicKey<Q, P, N>,
    prover_pk: PublicKey<Q, P, N>,
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
            r: (PolyR::<Q, N>::zero(), PolyR::<Q, N>::zero(), PolyR::<Q, N>::zero()),
            u: Vec::new(),
            y: Vec::new(),
        }
    }

    // generate-phase
    fn generate(&mut self) -> Ciphertext<Q, N> {
        let two = PolyR::<Q, N>::from(Zq::<Q>::from(2));
        
        // todo: use DGS
        let r = rand_tuple::<Q, N>(Some(two));
        let c = self.verifier_pk.encrypt(&self.message, r);
        
        self.r = r;

        c
    }

    // prove-phase, which is a sigma protocol with soundness amplification
    fn commit(&mut self, l: usize) -> Vec<Ciphertext<Q, N>> {
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
            
            let ui = Plaintext::<P, N> {
                poly: u[i].clone(),
                modulo: P,
            };
            w.push(self.verifier_pk.encrypt(&ui, y_times_2[i].clone())) 
        }

        self.u = u;
        
        w
    }

    fn reveal(&self, challenge: Vec<PolyR<P, N>>) -> (Vec<PolyR<P, N>>, Vec<TuplePolyR<Q, N>>) {
        let l = challenge.len();
        assert_eq!(l, self.u.len());
        let mut v: Vec<PolyR<P, N>> = Vec::new();
        let mut z: Vec<TuplePolyR<Q, N>> = Vec::new();

        for i in 0..l {
            let u_i = self.u[i].clone();
            let y_i = self.y[i].clone();
            let r = self.r.clone();
            let m = self.message.poly.clone();
            let gamma_i = challenge[i].clone();

            let v_i = u_i + gamma_i * m;

            let gamma_i = convert_ring::<P, Q, N>(gamma_i);
            let z_i = (
                y_i.0 + gamma_i * r.0, 
                y_i.1 + gamma_i * r.1,
                y_i.2 + gamma_i * r.2
            );

            v.push(v_i);
            z.push(z_i);
        }


        (v, z)
    }
}

impl<const Q: u64, const P: u64, const N: usize> Verifier<Q, P, N> {
    // instantiation
    fn new(prover_pk: PublicKey<Q, P, N>) -> Self {
        let secret_key =  SecretKey::new(); 
        let public_key = secret_key.pk_gen();

        Self {
            secret_key, 
            public_key, 
            prover_pk,
        }
    }

    fn challenge(&self, commitment: Vec<PolyR<P, N>>, l: usize) -> Vec<PolyR<P, N>> {
        // sample l monomials
        let mut rng = rand::thread_rng();
        let mut vec = Vec::new();

        for _ in 0..l {
            let mut coeffs: [Zq<P>; N] = [Zq::<P>::zero(); N];
            let index: usize = rng.gen_range(0..N);
            coeffs[index] = Zq::<P>::one();
            let coeffs = coeffs; // immutate it
            vec.push(PolyR::from(coeffs));
        }
        vec
    }

    fn verify(response: Vec<(PolyR<P, N>, PolyR<P, N>)>) -> bool {
        todo!();
    }
}
