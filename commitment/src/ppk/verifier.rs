#![allow(non_snake_case)]

use ark_serialize::CanonicalSerialize;
use lattirust_arithmetic::{challenge_set::{labrador_challenge_set::LabradorChallengeSet, ppk_challenge_set::PPKChallengeSet}, linear_algebra::Vector, nimue::{arthur::SerArthur, merlin::SerMerlin, serialization::{FromBytes, ToBytes}, traits::ChallengeFromRandomBytes}, ring::Zq, traits::FromRandomBytes};
use nimue::{ByteChallenges, BytePublic, ByteReader, IOPattern, IOPatternError, Merlin};
use rand::Rng;
use ark_ff::{One, Zero};
use ark_std::rand;
use relations::{ajtai_cm::Witness, principal_relation::PrincipalRelation};
use crate::bfv::{ciphertext::Ciphertext, plaintext::Plaintext, secret_key::SecretKey, public_key::PublicKey, util::{convert_ring, PolyR, TuplePolyR}};

use super::{nizk::new_ppk_io, util::PublicParameters};

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
    pub fn new_null() -> Self {
        let _secret_key =  SecretKey::new(); 
        let _public_key = _secret_key.pk_gen();

        Self {
            _secret_key, 
            _public_key, 
            prover_pk: PublicKey::default(),
            l: usize::default(), 
            c: Ciphertext::default(),
            w: Vec::default(),
            gamma: Vec::default(),
        }
    }

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

    pub fn nizk_verify(&mut self, merlin: &mut Merlin, l: usize) -> Result<bool, IOPatternError> { 
        // TODO: need to use create function like `public_parameters`?
        let pp = merlin.next_like(&PublicParameters::<Q, P, N>::default()).unwrap();
        println!("pp: {pp:?}");
        merlin.ratchet()?;

        let ctxt = merlin.next_like(&Ciphertext::<Q, N>::zero()).unwrap();
        merlin.ratchet()?;
        
        let w = merlin.next_vec::<Ciphertext::<Q, N>>(l).unwrap();
        let gamma = merlin.challenge_vec::<PolyR<P, N>, PPKChallengeSet<PolyR<P, N>>>(l).unwrap();
        
        let opening = merlin.next_like(&(vec![PolyR::<P, N>::zero(); l], vec![TuplePolyR::<Q, N>::zero(); l])).unwrap();

        // prepare all the ingredients for the verify        
        self.prover_pk = pp.pk;
        self.l = l;
        self.c = ctxt;
        self.w = w;
        self.gamma = gamma;

        Ok(self.verify(opening))
    }

    pub fn test_nizk(&mut self, merlin: &mut Merlin) -> Result<bool, IOPatternError> { 
        let poly = merlin.next_like_canonical_serializable::<PolyR::<Q, N>>(&PolyR::<Q, N>::zero()).unwrap();
        Ok(true)
    }
}
