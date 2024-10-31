#![allow(non_snake_case)]
#![allow(unused_imports)]

use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::linear_algebra::{Matrix, Vector};
use lattirust_arithmetic::ntt::ntt_modulus;

use lattirust_arithmetic::ring::{ConvertibleRing, PolyRing, Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT, SignedRepresentative, UnsignedRepresentative, Zq};
use rand::{CryptoRng, RngCore};
use rand_distr::Normal;
use rand::distributions::Distribution;
use ark_ff::{One, UniformRand, Zero};
use ark_std::rand;
use lattirust_arithmetic::traits::FromRandomBytes;

type PolyRq<const Q: u64, const N: usize> = Pow2CyclotomicPolyRing<Zq<Q>, N>;
type PolyRp<const P: u64, const N: usize> = Pow2CyclotomicPolyRing<Zq<P>, N>;
// TODO: use the same rng everywhere

// TODO: toy sampling, need to use OpenFHE code 
pub fn get_gaussian_vec<Rng: RngCore + CryptoRng, const Q: u64>(  
    std_dev: f64, 
    dimension: usize, 
    rng: &mut Rng
) -> Vec<Zq<Q>> {
    let gaussian = Normal::new(0.0, std_dev).unwrap();
    let val: Vec<Zq<Q>> = (0..dimension)
        .map(|_| Zq::<Q>::from(gaussian.sample(rng) as i64))
        .collect();

    val
}

pub fn get_gaussian<Rng: RngCore + CryptoRng, const Q: u64, const N: usize>(std_dev: f64, dimension: usize, rng: &mut Rng) -> PolyRq<Q, N> {
    let rand_vec: Vec<Zq<Q>> = get_gaussian_vec(std_dev, dimension, rng);
    let rand_vec: [Zq<Q>; N] = rand_vec.try_into().expect("Bad format");
    
    Pow2CyclotomicPolyRing::<Zq<Q>, N>::from(rand_vec)
}

pub fn convert_ring<const SOURCE: u64, const TARGET: u64, const N: usize>(poly: PolyRq<SOURCE, N>) -> PolyRp<TARGET, N> {
    let coeffs: Vec<Zq<SOURCE>> = poly.coeffs();
    let coeffs: Vec<Zq<TARGET>> = coeffs.into_iter()
        .map(|x| <Zq<TARGET>>::from(UnsignedRepresentative::from(x).0))
        .collect();
    let coeffs: [Zq<TARGET>; N] = coeffs.try_into().expect("Bad format");

    Pow2CyclotomicPolyRing::<Zq<TARGET>, N>::from(coeffs)
}

pub fn rand_ternary_poly<Rng: RngCore + CryptoRng, const P: u64, const N: usize>(size: usize, rng: &mut Rng) -> PolyRp<P, N> {
    let bytes: Vector<u8> = Vector::<u8>::rand(size, rng);
    WeightedTernaryChallengeSet::<Pow2CyclotomicPolyRing::<Zq<P>, N>>::try_from_random_bytes(bytes.as_slice()).unwrap()
}

pub fn rand_tuple<const Q: u64, const P: u64, const N: usize>(factor: Option<PolyRq<Q, N>>) ->  (PolyRq<Q, N>, PolyRq<Q, N>, PolyRq<Q, N>) {
    let mut rng = rand::thread_rng();
    let size = WeightedTernaryChallengeSet::<PolyRq<Q, N>>::byte_size();
    let f = factor.unwrap_or_else(|| Pow2CyclotomicPolyRing::<Zq<Q>, N>::one());

    (f * rand_ternary_poly(size, &mut rng), 
    f * get_gaussian(3.2, N, &mut rng), 
    f * get_gaussian(3.2, N, &mut rng))
}

pub struct PublicKey<const Q: u64, const P: u64, const N: usize> {
    pub poly1: PolyRq<Q, N>, 
    pub poly2: PolyRq<Q, N>,
    pub modulo: u64,
}
pub struct SecretKey<const Q: u64, const P: u64, const N: usize> {
    pub poly: PolyRp<P, N>,  
    pub modulo: u64,
}

pub struct Plaintext<const P: u64, const N: usize> {
    pub poly: PolyRp<P, N>,
    pub modulo: u64,
}

pub struct Ciphertext<const Q: u64, const N: usize> {
    pub c1: PolyRq<Q, N>,
    pub c2: PolyRq<Q, N>,
    pub modulo: u64,
}

impl<const Q: u64, const P: u64, const N: usize> SecretKey<Q, P, N> {
    pub fn new() -> Self {
        // generate random bytes to sample a ternary secret
        let size = WeightedTernaryChallengeSet::<PolyRp<P, N>>::byte_size();

        Self {
            poly: rand_ternary_poly(size, &mut rand::thread_rng()),
            modulo: P,
        }
    }

    pub fn pk_gen(&self) -> PublicKey<Q, P, N> {
        let mut rng = rand::thread_rng();

        // e should have small std_dev (how small?) for the correctness, TODO: check the parameters
        // TODO: use OpenFHE DGS
        let e: PolyRq<Q, N> = get_gaussian(3.2, N, &mut rng);

        // convert sk in PolyRp to PolyRq in order to perform the operation in PolyRq
        let sk_zq: PolyRq<Q, N> = convert_ring::<P, Q, N>(self.poly.clone());

        // compute the actual pk pair
        let poly2: PolyRq<Q, N> = PolyRq::rand(&mut rng);
        let poly1: PolyRq<Q, N> = -(poly2.clone() * sk_zq + e);

        PublicKey {
            poly1,
            poly2,
            modulo: Q,
        }
    }
    
    pub fn decrypt(&self, c: Ciphertext<Q, N>) -> Plaintext<P, N> {
        let c1: PolyRq<Q, N> = c.c1.clone();
        let c2: PolyRq<Q, N> = c.c2.clone();
        let sk_zq: PolyRq<Q, N> = convert_ring::<P, Q, N>(self.poly.clone());
        let raw: PolyRq<Q, N> = c1 + c2 * sk_zq;

        let p = self.modulo as f64;
        let q = c.modulo as f64;
        let delta = p/q; 

        let coeffs: Vec<Zq<P>> = raw
            .coeffs()
            .into_iter()
            .map(|x| <Zq<P>>::from((delta * (UnsignedRepresentative::from(x).0 as f64)) as u128))
            .collect();
        let coeffs: [Zq<P>; N] = coeffs.try_into().expect("Bad format");

        Plaintext {
            poly: PolyRp::from(coeffs),
            modulo: P,
        }

    }
}

impl<const Q: u64, const P: u64, const N: usize> PublicKey<Q, P, N> {
    pub fn encrypt(
        &self, 
        m: &Plaintext<P, N>, 
        r:  (PolyRq<Q, N>, PolyRq<Q, N>, PolyRq<Q, N>)
        ) -> Ciphertext<Q, N> {
        let (pk1, pk2) = (self.poly1.clone(), self.poly2.clone());
        let (u, e1, e2) = r;
        
        let p = m.modulo;
        let q = self.modulo;
        let delta = (q as f64 / p as f64).floor() as u128; // round off

        // mutliply each coeff of m with delta as bigint, and convert the polynomial into Zq, or convert m into PolyRq and multiply, but attention overflow
        // TODO: define it as scalar - poly multiplication 
        let coeffs_zq: Vec<Zq<Q>> = m
            .poly
            .coeffs()
            .into_iter()
            .map(|x| <Zq<Q>>::from(delta * UnsignedRepresentative::from(x).0))
            .collect();
        let coeffs_zq: [Zq<Q>; N] = coeffs_zq.try_into().expect("Bad format");
        let m_delta= Pow2CyclotomicPolyRing::<Zq<Q>, N>::from(coeffs_zq);

        // compute a, b
        let c1: PolyRq<Q, N> = pk1 * u.clone() + e1 + m_delta;
        let c2: PolyRq<Q, N> = pk2 * u.clone() + e2;

        // return the ciphertext
        Ciphertext {
            c1,
            c2,
            modulo: Q,
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
    // type PolyRq = Zq<Q>;
    // type PolyRp = Zq<P>;
    // type Pol PolyRq = Pow2CyclotomicPolyRing PolyRq, N>;
    // type PolyPolyRp = Pow2CyclotomicPolyRing<PolyRp, N>;

    let sk: SecretKey<Q, P, N> = SecretKey::new();
    let pk = sk.pk_gen();
    let ptxt = Plaintext::rand_message();
    let ctxt = pk.encrypt(&ptxt, rand_tuple::<Q, P, N>(None));
    
    let act = sk.decrypt(ctxt).poly;
    // println!("act {act:?}");
    let exp = ptxt.poly;
    // println!("exp {exp:?}");
    assert_eq!(act, exp);
}