#![allow(non_snake_case)]

use lattirust_arithmetic::ring::{PolyRing, Pow2CyclotomicPolyRing, SignedRepresentative, Zq};
use super::ciphertext::Ciphertext;
use super::plaintext::Plaintext;
use super::util::TuplePolyR;

type PolyR<const M: u64, const N: usize> = Pow2CyclotomicPolyRing<Zq<M>, N>;
// TODO: use the same rng everywhere

pub struct PublicKey<const Q: u64, const P: u64, const N: usize> {
    pub poly1: PolyR<Q, N>, 
    pub poly2: PolyR<Q, N>,
    pub modulo: u64,
}

impl<const Q: u64, const P: u64, const N: usize> PublicKey<Q, P, N> {
    pub fn encrypt(
        &self, 
        m: &Plaintext<P, N>, 
        r:  TuplePolyR<Q, N>
        ) -> Ciphertext<Q, N> {
        let (pk1, pk2) = (self.poly1.clone(), self.poly2.clone());
        let TuplePolyR(r0, r1, r2) = r;
        
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
