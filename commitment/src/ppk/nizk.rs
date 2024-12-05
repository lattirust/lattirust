#![allow(non_snake_case)]

use ark_ff::{UniformRand, Zero};
use lattirust_arithmetic::challenge_set::ppk_challenge_set::PPKChallengeSet;
use lattirust_arithmetic::nimue::iopattern::{SerIOPattern, SqueezeFromRandomBytes, RatchetIOPattern};
use lattirust_arithmetic::nimue::arthur::SerArthur;
use lattirust_arithmetic::ring::PolyRing;
use lattirust_arithmetic::traits::FromRandomBytes;
use nimue::{DefaultHash, IOPattern};

use crate::bfv::util::{PolyR, TuplePolyR};
use crate::bfv::{ciphertext, Ciphertext};

use super::util::PublicParameters;

pub struct PPKRelation<const Q: u64, const P: u64, const N: usize> {
    m: PolyR<P, N>,
    r: TuplePolyR<Q, N>,
    c: Ciphertext<Q, N>,
}

pub struct BaseTranscript<'a, const Q: u64, const P: u64, const N: usize> {
    pub instance: &'a PPKRelation<Q, P, N>,
    pub commitment: &'a Vec<Ciphertext<Q, N>>,
    pub challenge: &'a Vec<PolyR<P, N>>,
    pub opening: &'a Vec<(PolyR<P, N>, PolyR<P, N>)>,
}

impl<'a, const Q: u64, const P: u64, const N: usize> BaseTranscript<'a, Q, P, N> {
    pub fn new(
        instance: &'a PPKRelation<Q, P, N>,
        commitment: &'a Vec<Ciphertext<Q, N>>,
        challenge: &'a Vec<PolyR<P, N>>,
        opening: &'a Vec<(PolyR<P, N>, PolyR<P, N>)>,
    ) -> Self {
        BaseTranscript {
            instance,
            commitment,
            challenge,
            opening,
        }
    }
}


// TODO: check again
// generate the IOPattern
pub fn new_ppk_io<const Q: u64, const P: u64, const N: usize>(l: usize) -> IOPattern {
    IOPattern::<DefaultHash>::new("PPK_protocol")
    // public parameters
    .absorb_serializable_like(&PublicParameters::<Q, P, N>::default(), "pp")
    .ratchet()
    // ciphertext
    .absorb_canonical_serializable_like(&Ciphertext::<Q, N>::zero(), "ciphertext") // c: Ciphertext<Q, N>
    .ratchet()
    // commitment
    .absorb_serializable_like(&vec![Ciphertext::<Q, N>::zero(); l], "commitment") // w_i: Vec<Ciphertext<Q, N>>
    // challenge
    .squeeze(PPKChallengeSet::<PolyR<Q, N>>::byte_size() * l , "challenge") // gamma_i
    // opening
    .absorb_serializable_like(&(vec![PolyR::<P, N>::zero(); l], vec![TuplePolyR::<Q, N>::zero(); l]), "opening") // (v_i, z_i): Vec<PolyR<P, N>>, Vec<TuplePolyR<Q, N>>
}

pub fn test_io<const Q: u64, const P: u64, const N: usize>() -> IOPattern {
    IOPattern::<DefaultHash>::new("test")
    .absorb_canonical_serializable_like::<PolyR::<Q, N>>(&PolyR::<Q, N>::default(), "test")
    // .absorb_canonical_serializable_like(&Ciphertext::<Q, N>::default(), "test2")
}
