#![allow(non_snake_case)]

use lattirust_arithmetic::nimue::iopattern::{SerIOPattern, SqueezeFromRandomBytes, RatchetIOPattern};
use nimue::{DefaultHash, IOPattern};

use crate::bfv::util::{PolyR, TuplePolyR};
use crate::bfv::{ciphertext, Ciphertext};

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

// pub trait PPKIOPattern
// where 
//     Self: SerIOPattern + SqueezeFromRandomBytes + RatchetIOPattern,
// {}


// TODO: check again
// generate the IOPattern
pub fn new_ppk_io(ctxt_size: usize, comm_size: usize, chall_size: usize, resp_size: usize) -> IOPattern {
    IOPattern::<DefaultHash>::new("PPK_protocol")
    .absorb(ctxt_size, "ciphertext") // c
    .absorb(comm_size, "commitment") // w_i
    .squeeze(chall_size, "challenge") // gamma_i
    .absorb(resp_size, "response") // (v_i, z_i)
}
