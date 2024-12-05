#![allow(non_snake_case)]
#![allow(unused_must_use)]

use std::ptr::read;

use ark_serialize::CanonicalSerialize;
use lattirust_arithmetic::{challenge_set::{self, ppk_challenge_set::PPKChallengeSet}, nimue::{arthur::SerArthur, serialization::ToBytes, traits::ChallengeFromRandomBytes}, ring::Zq};
use ark_ff::{UniformRand, Zero};
use ark_std::rand;
use nimue::{Arthur, ByteChallenges, ByteWriter, IOPattern};
use crate::{bfv::{ciphertext::Ciphertext, plaintext::Plaintext, public_key::PublicKey, secret_key::SecretKey, util::{convert_ring, rand_tuple, PolyR, TuplePolyR}}, ppk::util::PublicParameters};
use super::nizk::new_ppk_io;

pub struct Prover<const Q: u64, const P: u64, const N: usize> {
    // params: ParamsBFV,
    _secret_key: SecretKey<Q, P, N>, // polynomial w/ coefficients in {-1, 0, +1}
    public_key: PublicKey<Q, P, N>, // ring mod q
    message: Plaintext<P, N>, // ring mod p
    r: TuplePolyR<Q, N>,
    u: Vec<PolyR<P, N>>,
    y: Vec<TuplePolyR<Q, N>>,
    l: usize,
}

impl<const Q: u64, const P: u64, const N: usize> Prover<Q, P, N> {
    // instantiation
    pub fn new(message: Plaintext<P, N>, l: usize) -> Self {
        let _secret_key =  SecretKey::new(); 
        let public_key = _secret_key.pk_gen();
        
        Self {
            _secret_key, 
            public_key, 
            message,
            r: TuplePolyR::<Q, N>::default(),
            u: Vec::default(),
            y: Vec::default(),
            l
        }
    }

    pub fn return_pk(&self) -> PublicKey<Q, P, N> {
        PublicKey {
            poly1: self.public_key.poly1.clone(),
            poly2: self.public_key.poly2.clone(),
            modulo: Q,
        }
    }

    pub fn public_parameters(&self) -> PublicParameters<Q, P, N> {
        PublicParameters {
            pk: self.public_key.clone(),
            q: Q,
            p: P,
            n: N,
        }
    }

    // generate-phase
    // TODO: use RefCell and make self immutable (?)
    pub fn generate(&mut self) -> Ciphertext<Q, N> {
        let two = PolyR::<Q, N>::from(Zq::<Q>::from(2));

        // todo: use DGS
        let r = rand_tuple::<Q, N>(None);
        self.r = r.clone();
        let r = r * two;

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
            self.y.push(yi.clone());
            y_times_2.push(yi.clone() * two);
            
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
            let z_i = y_i + r.clone() * gamma_i;

            v.push(v_i);
            z.push(z_i);
        }

        (v, z)
    }

    pub fn nizk_prove(&mut self, arthur: &mut Arthur) -> Vec<u8> {
        let l = 6;
        let pp = PublicParameters { 
            pk: self.public_key,
            q: Q,
            p: P,
            n: N,
        };
        arthur
            .absorb_serializable(&pp)
            .expect("error in absorbing the public parameters");
        arthur.ratchet();

        let ctxt = self.generate();
        arthur
            .absorb_serializable(&ctxt.clone())
            .expect("error in absorbing the ciphertext");
        debug_assert_eq!(Ciphertext::<Q, N>::zero().to_bytes().unwrap().len(), ctxt.clone().to_bytes().unwrap().len());
        arthur.ratchet().unwrap();

        let commitment = self.commit(l);
        debug_assert_eq!(vec![Ciphertext::<Q, N>::zero(); l].to_bytes().unwrap().len(), commitment.to_bytes().unwrap().len());

        arthur
            .absorb_serializable(&commitment.clone())
            .expect("error in absorbing the commitment");
        
        // get the challenge
        let challenge: Vec<PolyR<P, N>> = arthur.challenge_vec::<PolyR<P, N>, PPKChallengeSet<PolyR<P, N>>>(l).unwrap();
        debug_assert_eq!(vec![PolyR::<P, N>::zero(); l].to_bytes().unwrap().len(), challenge.to_bytes().unwrap().len());

        let opening = self.reveal(challenge);
        debug_assert_eq!((vec![PolyR::<P, N>::zero(); l], vec![TuplePolyR::<Q, N>::zero(); l]).to_bytes().unwrap().len(), opening.to_bytes().unwrap().len());
        arthur
            .absorb_serializable(&opening.clone())
            .expect("error in absorbing the opening");

        let transcript = arthur.transcript().to_vec();
        transcript
    }

    pub fn test_nizk(&mut self, arthur: &mut Arthur) -> Vec<u8> {
        let ctxt = PolyR::<Q, N>::default();
        arthur.absorb_canonical_serializable::<PolyR::<Q, N>>(&ctxt);
        // arthur.absorb_serializable(&ctxt);
        arthur.transcript().to_vec()
    }
}
