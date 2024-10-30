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

// TODO: toy sampling, need to use OpenFHE code 
pub fn get_gaussian_vec<
    T: RngCore + CryptoRng, 
    const Q: u64
>(  std_dev: f64, 
    dimension: usize, 
    rng: &mut T
) -> Vec<Zq<Q>> {
    // TODO: modulo the coefficients
    let gaussian = Normal::new(0.0, std_dev).unwrap();
    let val: Vec<Zq<Q>> = (0..dimension)
        .map(|_| Zq::<Q>::from(gaussian.sample(rng) as i64))
        .collect();

    val
}


// TODO: do I need both modules everywhere?
// move it to a dedicated library later on
pub struct PublicKey<const Q: u64, const P: u64, const N: usize> {
    // params: ParamsBFV,
    pub poly1: Pow2CyclotomicPolyRing<Zq<Q>, N>, 
    pub poly2: Pow2CyclotomicPolyRing<Zq<Q>, N>,
    pub modulo: u64,
}
pub struct SecretKey<const Q: u64, const P: u64, const N: usize> {
    // params: ParamsBFV,
    pub poly: Pow2CyclotomicPolyRing<Zq<P>, N>,  
    pub modulo: u64,
}

pub struct Plaintext<const P: u64, const N: usize> {
    // params: ParamsBFV,
    pub poly: Pow2CyclotomicPolyRing<Zq<P>, N>,
    pub modulo: u64,
}

pub struct Ciphertext<const Q: u64, const N: usize> {
    // params: ParamsBFV,
    pub c1: Pow2CyclotomicPolyRing<Zq<Q>, N>,
    pub c2: Pow2CyclotomicPolyRing<Zq<Q>, N>,
    pub modulo: u64,
}

pub fn q_to_p_ring<const Q: u64, const P: u64, const N: usize>(poly: Pow2CyclotomicPolyRing<Zq<Q>, N>) -> Pow2CyclotomicPolyRing<Zq<P>, N> {
    let coeffs: Vec<Zq<Q>> = poly.coeffs();
    let coeffs: Vec<Zq<P>> = coeffs.into_iter()
        .map(|x| <Zq<P>>::from(UnsignedRepresentative::from(x).0))
        .collect();
    let coeffs: [Zq<P>; N] = coeffs.try_into().expect("Bad format");
    Pow2CyclotomicPolyRing::<Zq<P>, N>::from(coeffs)
}
pub fn p_to_q_ring<const Q: u64, const P: u64, const N: usize>(poly: Pow2CyclotomicPolyRing<Zq<P>, N>) -> Pow2CyclotomicPolyRing<Zq<Q>, N> {
    let coeffs: Vec<Zq<P>> = poly.coeffs();
    let coeffs: Vec<Zq<Q>> = coeffs.into_iter()
        .map(|x| <Zq<Q>>::from(UnsignedRepresentative::from(x).0))
        .collect();
    let coeffs: [Zq<Q>; N] = coeffs.try_into().expect("Bad format");
    Pow2CyclotomicPolyRing::<Zq<Q>, N>::from(coeffs)
}

pub fn rand_ternary_poly<const P: u64, const N: usize>(size: usize) -> Pow2CyclotomicPolyRing<Zq<P>, N> {
    let bytes = Vector::<u8>::rand(size, &mut rand::thread_rng());
    
    WeightedTernaryChallengeSet::<Pow2CyclotomicPolyRing::<Zq<P>, N>>::try_from_random_bytes(bytes.as_slice()).unwrap()
}

impl<const Q: u64, const P: u64, const N: usize> SecretKey<Q, P, N> {
    pub fn new() -> Self {
        // generate random bytes to sample a ternary secret
        let size = WeightedTernaryChallengeSet::<Pow2CyclotomicPolyRing<Zq<P>, N>>::byte_size();

        Self {
            // params,
            poly: rand_ternary_poly(size),
            modulo: P,
        }
    }

    pub fn pk_gen(&self) -> PublicKey<Q, P, N> {
        let mut rng = rand::thread_rng();
        type Rq<const Q: u64, const N: usize> = Pow2CyclotomicPolyRing::<Zq<Q>, N>;
        // let size = Pow2CyclotomicPolyRing::<Zq<Q>, N>::byte_size();
        
        // e should have small std_dev (how small?), TODO: check the parameters
        // TODO: use OpenFHE DGS
        // let e = get_poly_from_gaussian(1.0, size, &mut rng);
        let e: Pow2CyclotomicPolyRing::<Zq<Q>, N> = Rq::rand(&mut rng);

        // convert sk in Rp to Rq in order to perform the operation in Rq
        let sk_zq = p_to_q_ring(self.poly.clone());
        // compute the actual pk pair
        let pk2: Pow2CyclotomicPolyRing<Zq<Q>, N> = Rq::rand(&mut rng);
        let pk1 = -(pk2 * sk_zq) + e;

        PublicKey {
            // params,
            poly1: pk1,
            poly2: pk2,
            modulo: Q,
        }
    }
    
    pub fn decrypt(&self, c: Ciphertext<Q, N>) -> Plaintext<P, N> {
        type Rp<const P: u64, const N: usize> = Pow2CyclotomicPolyRing::<Zq<P>, N>;
        let a = c.c1.clone();
        let b = c.c2.clone();

        // convert the elements to the other field
        // TODO: DRY
        let a_zp = q_to_p_ring(a);
        let b_zp = q_to_p_ring(b);

        let p = self.modulo as f64;
        let q = c.modulo as f64;
        let delta = p/q; // round off

        let raw: Pow2CyclotomicPolyRing<Zq<P>, N> = a_zp + b_zp * self.poly;
        let coeffs: Vec<Zq<P>> = raw.coeffs();
        let coeffs: Vec<Zq<P>> = coeffs.into_iter()
            .map(|x| <Zq<P>>::from((delta * (UnsignedRepresentative::from(x).0 as f64)) as u128))
            .collect();
        let coeffs: [Zq<P>; N] = coeffs.try_into().expect("Bad format");
        let poly = Rp::from(coeffs);

        // println!("poly: {poly:?} \n");
        Plaintext {
            poly,
            modulo: P,
        }

    }
}

impl<const Q: u64, const P: u64, const N: usize> PublicKey<Q, P, N> {
    pub fn encrypt(&self, m: &Plaintext<P, N>, r: (Pow2CyclotomicPolyRing<Zq<Q>, N>, Pow2CyclotomicPolyRing<Zq<Q>, N>, Pow2CyclotomicPolyRing<Zq<Q>, N>)) -> Ciphertext<Q, N> {
        // let mut rng = rand::thread_rng();
        let pk1 = self.poly1.clone();
        let pk2 = self.poly2.clone();
        let (u, e1, e2) = r;
        
        let p = m.modulo;
        let q = self.modulo;

        let delta = (q as f64 / p as f64).floor() as u128; // round off

        // mutliply each coeff of m with delta as bigint, and convert the polynomial into Zq, or convert m into Rq and multiply, but attention overflow
        let coeffs: Vec<Zq<P>> = m.poly.coeffs();
        let coeffs_zq: Vec<Zq<Q>> = coeffs.into_iter()
            .map(|x| <Zq<Q>>::from(delta * UnsignedRepresentative::from(x).0))
            .collect();
        let coeffs_zq: [Zq<Q>; N] = coeffs_zq.try_into().expect("Bad format");
        let m_delta = Pow2CyclotomicPolyRing::<Zq<Q>, N>::from(coeffs_zq);

        // compute a, b
        let c1: Pow2CyclotomicPolyRing<Zq<Q>, N> = pk1 * u + e1 + m_delta;
        let c2: Pow2CyclotomicPolyRing<Zq<Q>, N> = pk2 * u + e2;

        // return the ciphertext
        Ciphertext {
            c1,
            c2,
            modulo: Q,
        }
    }

    pub fn rand_tuple(factor: Option<Pow2CyclotomicPolyRing<Zq<Q>, N>>) 
    -> (Pow2CyclotomicPolyRing<Zq<Q>, N>, 
        Pow2CyclotomicPolyRing<Zq<Q>, N>, 
        Pow2CyclotomicPolyRing<Zq<Q>, N>) {
        type Rq<const Q: u64, const N: usize> = Pow2CyclotomicPolyRing::<Zq<Q>, N>;
        let mut rng = rand::thread_rng();
        let size = WeightedTernaryChallengeSet::<Pow2CyclotomicPolyRing<Zq<Q>, N>>::byte_size();

        match factor {
            Some(f) => (
                f * rand_ternary_poly(size), 
                f * Rq::rand(&mut rng), 
                f * Rq::rand(&mut rng)),

            None => (
                rand_ternary_poly(size), 
                Rq::rand(&mut rng), 
                Rq::rand(&mut rng)),
        }
    }

}

impl<const P: u64, const N: usize> Plaintext<P, N> {
    pub fn rand_message() -> Self {
        Self {
            poly: Pow2CyclotomicPolyRing::<Zq<P>, N>::rand(&mut rand::thread_rng()),
            modulo: P,
        }
    }
}

#[test]
fn test() {
    const N: usize = 128;
    const Q: u64 = ntt_modulus::<N>(17);
    const P: u64 = ntt_modulus::<N>(15);
    // type Rq = Zq<Q>;
    // type Rp = Zq<P>;
    // type PolyRq = Pow2CyclotomicPolyRing<Rq, N>;
    // type PolyRp = Pow2CyclotomicPolyRing<Rp, N>;

    let sk: SecretKey<Q, P, N> = SecretKey::new();
    let pk = sk.pk_gen();
    let ptxt = Plaintext::rand_message();
    let ctxt = pk.encrypt(&ptxt, PublicKey::<Q, P, N>::rand_tuple(None));
    
    let act = sk.decrypt(ctxt).poly;
    // println!("act {act:?}");
    let exp = ptxt.poly;
    // println!("exp {exp:?}");
    assert_eq!(act, exp);
}