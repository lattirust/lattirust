#![allow(non_snake_case)]

use lattirust_arithmetic::{challenge_set::labrador_challenge_set::LabradorChallengeSet, nimue::serialization::ToBytes, ring::Zq};
use nimue::ByteChallenges;
use rand::Rng;
use ark_ff::{One, Zero};
use ark_std::rand;
use relations::{ajtai_cm::Witness, principal_relation::PrincipalRelation};
use crate::bfv::{ciphertext::Ciphertext, plaintext::Plaintext, secret_key::SecretKey, public_key::PublicKey, util::{convert_ring, PolyR, TuplePolyR}};

use super::nizk::new_ppk_io;

pub struct Verifier<const Q: u64, const P: u64, const N: usize> {
    // params: ParamsBFV,
    _secret_key: SecretKey<Q, P, N>,
    _public_key: PublicKey<Q, P, N>,
    prover_pk: PublicKey<Q, P, N>,
    l: usize,
    c: Ciphertext<Q, N>, 
    w: Vec<Ciphertext<Q, N>>,
    gamma: Vec<PolyR<P, N>>,
}

impl<const Q: u64, const P: u64, const N: usize> Verifier<Q, P, N> {
    // instantiation
    pub fn new(prover_pk: PublicKey<Q, P, N>, l: usize) -> Self {
        let _secret_key =  SecretKey::new(); 
        let _public_key = _secret_key.pk_gen();

        Self {
            _secret_key, 
            _public_key, 
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
            let ctxt = pk.encrypt(&vi, z[i].clone() * two);
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

    pub fn nizk(&self, commitment: Vec<u8>, challenge_prover: Vec<u8>) -> bool {
        let io = new_ppk_io(32, 32, 32, 32);
        // instantiate a prover
        let w = commitment.to_bytes().unwrap().clone();
        let mut merlin = io.to_merlin(&w);
        let challenge_verifier = merlin.challenge_bytes::<32>().unwrap();
        
        challenge_prover == challenge_verifier
    }
}
