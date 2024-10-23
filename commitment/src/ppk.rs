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

pub struct ParamsBFV {
    P: u64, // ptxt modulo 
    Q: u64, // ctxt modulo
    N: usize, // degree of the polynomial
}

// TODO: toy sampling, need to use OpenFHE code 
pub fn get_gaussian_vec<T: RngCore + CryptoRng, const Q: u64>(std_dev: f64, dimension: usize, rng: &mut T) -> Vec<Zq<Q>> {
    // TODO: modulo the coefficients
    let gaussian = Normal::new(0.0, std_dev).unwrap();
    let val: Vec<Zq<Q>> = (0..dimension)
        .map(|_| Zq::<Q>::from(gaussian.sample(rng) as i64))
        .collect();

    val
}

// pub fn get_poly_from_gaussian<Rng: RngCore + CryptoRng, const Q: u64, const N: usize>(std_dev: f64, dimension: usize, rng: &mut Rng) -> Pow2CyclotomicPolyRing<Zq<Q>, N> {
//     // let rand_vec = get_gaussian_vec(std_dev, dimension, rng);

//     // Pow2CyclotomicPolyRing::<Zq<Q>, N>::from(rand_vec)
//     // Pow2CyclotomicPolyRing::<Zq<Q>, N>::rand(&mut rng)
// }

// define Q, P, N before using 
pub struct Prover<const Q: u64, const P: u64, const N: usize> {
    // params: ParamsBFV,
    secret_key: SecretKey<Q, P, N>, // polynomial w/ coefficients in {-1, 0, +1}
    public_key: PublicKey<Q, P, N>, // ring mod q
    verifier_pk: PublicKey<Q, P, N>,
    message: Pow2CyclotomicPolyRing<Zq<Q>, N>, // ring mod p
}

pub struct Verifier<const Q: u64, const P: u64, const N: usize> {
    // params: ParamsBFV,
    secret_key: SecretKey<Q, P, N>,
    public_key: PublicKey<Q, P, N>,
    commitment: Option<Vec<Pow2CyclotomicPolyRing<Zq<Q>, N>>>
}

// TODO: do I need both modules everywhere?
// move it to a dedicated library later on
pub struct PublicKey<const Q: u64, const P: u64, const N: usize> {
    // params: ParamsBFV,
    poly1: Pow2CyclotomicPolyRing<Zq<Q>, N>, 
    poly2: Pow2CyclotomicPolyRing<Zq<Q>, N>,
    pub modulo: u64,
}
pub struct SecretKey<const Q: u64, const P: u64, const N: usize> {
    // params: ParamsBFV,
    poly: Pow2CyclotomicPolyRing<Zq<P>, N>,  
    pub modulo: u64,
}

pub struct Plaintext<const P: u64, const N: usize> {
    // params: ParamsBFV,
    poly: Pow2CyclotomicPolyRing<Zq<P>, N>,
    pub modulo: u64,
}

pub struct Ciphertext<const Q: u64, const N: usize> {
    // params: ParamsBFV,
    c1: Pow2CyclotomicPolyRing<Zq<Q>, N>,
    c2: Pow2CyclotomicPolyRing<Zq<Q>, N>,
    pub modulo: u64,
}

impl<const Q: u64, const P: u64, const N: usize> SecretKey<Q, P, N> {
    pub fn new() -> Self {
        let mut rng = rand::thread_rng();
        
        // generate random bytes to sample a ternary secret
        // TODO: make a function
        let size = WeightedTernaryChallengeSet::<Pow2CyclotomicPolyRing<Zq<P>, N>>::byte_size();
        let bytes = Vector::<u8>::rand(size, &mut rng);

        let sk: Pow2CyclotomicPolyRing::<Zq<P>, N> = WeightedTernaryChallengeSet::<Pow2CyclotomicPolyRing<Zq<P>, N>>::try_from_random_bytes(bytes.as_slice()).unwrap();

        Self {
            // params,
            poly: sk,
            modulo: P,
        }
    }

    pub fn pk_gen(&self) -> PublicKey<Q, P, N> {
        let mut rng = rand::thread_rng();

        let size = Pow2CyclotomicPolyRing::<Zq<Q>, N>::byte_size();
        
        // e should have small std_dev (how small?), TODO: check the parameters
        // TODO: use OpenFHE DGS
        // let e = get_poly_from_gaussian(1.0, size, &mut rng);
        let e = Pow2CyclotomicPolyRing::<Zq<Q>, N>::rand(&mut rng);
        println!("{e:?}");

        // let bytes: Vector<u8> = Vector::<u8>::rand(size, &mut rng);
        // println!("{bytes:?}");

        // convert sk in Rp to Rq in order to perform the operation in Rq
        let coeffs: Vec<Zq<P>> = self.poly.coeffs();
        let coeffs: Vec<Zq<Q>> = coeffs.into_iter()
            .map(|x| <Zq<Q>>::from(UnsignedRepresentative::from(x).0))
            .collect();
        
        let sk_zq = Pow2CyclotomicPolyRing::<Zq<Q>, N>::from(coeffs);
        
        // compute the actual pk pair
        // let pk2: Pow2CyclotomicPolyRing::<Zq<Q>, N> = Pow2CyclotomicPolyRing::<Zq<Q>, N>::try_from_random_bytes(bytes.as_slice()).unwrap();
        let pk2 = Pow2CyclotomicPolyRing::<Zq<Q>, N>::rand(&mut rng);
        let pk1: Pow2CyclotomicPolyRing::<Zq<Q>, N> = -(pk2 * sk_zq) + e;

        PublicKey {
            // params,
            poly1: pk1,
            poly2: pk2,
            modulo: Q,
        }
    }
    
    fn decrypt(&self, c: Ciphertext<Q, N>) -> Plaintext<P, N> {
        let a = c.c1;
        let b = c.c2;

        // convert the elements to the other field
        let coeffs: Vec<Zq<Q>> = a.coeffs();
        let coeffs: Vec<Zq<P>> = coeffs.into_iter()
        .map(|x| <Zq<P>>::from(SignedRepresentative::from(x).0))
        .collect();
        let a_zp = Pow2CyclotomicPolyRing::<Zq<P>, N>::from(coeffs);
        
        let coeffs: Vec<Zq<Q>> = b.coeffs();
        let coeffs: Vec<Zq<P>> = coeffs.into_iter()
        .map(|x| <Zq<P>>::from(SignedRepresentative::from(x).0))
        .collect();
        let b_zp = Pow2CyclotomicPolyRing::<Zq<P>, N>::from(coeffs);

        let p = self.modulo;
        let q = c.modulo;

        let delta = ((p as f64) / (q as f64)).floor() as u128; // round off
        let delta = Pow2CyclotomicPolyRing::<Zq<P>, N>::from(delta);

        let poly = delta * (a_zp + b_zp * self.poly);

        Plaintext {
            poly,
            modulo: p,
        }

    }
}

impl<const Q: u64, const P: u64, const N: usize> PublicKey<Q, P, N> {
    pub fn encrypt(&self, m: &Plaintext<P, N>) -> Ciphertext<Q, N> {
        let mut rng = rand::thread_rng();
        let pk1 = self.poly1;
        let pk2 = self.poly2;

        // samples u, e1, e2 in Rq
        let u = Pow2CyclotomicPolyRing::<Zq<Q>, N>::rand(&mut rng);
        let e1 = Pow2CyclotomicPolyRing::<Zq<Q>, N>::rand(&mut rng);
        let e2 = Pow2CyclotomicPolyRing::<Zq<Q>, N>::rand(&mut rng);
        
        let p = m.modulo;
        let q = self.modulo;

        let delta = ((q as f64) / (p as f64)).floor() as u128; // round off

        // mutliply each coeff of m with delta as bigint, and convert the polynomial into Zq, or convert m into Rq and multiply, but attention overflow
        let coeffs: Vec<Zq<P>> = m.poly.coeffs();
        let coeffs_zq: Vec<Zq<Q>> = coeffs.into_iter()
            .map(|x| <Zq<Q>>::from(delta * UnsignedRepresentative::from(x).0))
            .collect();
        let m_delta = Pow2CyclotomicPolyRing::<Zq<Q>, N>::from(coeffs_zq);

        // compute a, b
        let c1 = pk1 * u + e1 + m_delta;
        let c2 = pk2 * u + e2;

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
        let mut rng = rand::thread_rng();
        let poly = Pow2CyclotomicPolyRing::<Zq<P>, N>::rand(&mut rng);

        Self {
            poly, 
            modulo: P,
        }
    }
}

impl<const Q: u64, const P: u64, const N: usize> Prover<Q, P, N> {
    // instantiation
    fn new(verifier_pk: PublicKey<Q, P, N>, message: Pow2CyclotomicPolyRing<Zq<Q>, N>) -> Self {
        let secret_key =  SecretKey::new(); 
        let public_key = secret_key.pk_gen();
        
        Self {
            secret_key, 
            public_key, 
            verifier_pk,
            message,
        }
    }

    // generate-phase
    fn generate(&self) -> Pow2CyclotomicPolyRing<Zq<P>, N> {
        todo!();
        // todo: use DGS
        // let mut Rp: (Pow2CyclotomicPolyRing<Rp, N>, Pow2CyclotomicPolyRing<Rp, N>, Pow2CyclotomicPolyRing<Rp, N>) = (0.0, 0.0, 0.0);
        // todo: define scalar-vector poly multiplication(?)
        // let c = self.verifier_pk.encrypt(self.message, 2 * Rp);
        
        // todo: send message
    }

    // prove-phase, which is a sigma protocol
    fn commit(&self, l: usize) -> Vec<Pow2CyclotomicPolyRing<Zq<P>, N>> {
        todo!();

        // let mut u: Vec<Pow2CyclotomicPolyRing<Rp, N>> = Vec::new(); // vec poly of size l
        // let mut y: Vec<Pow2CyclotomicPolyRing<Rp, N>> = Vec::new(); // vec of 3-tuple poly of size l
        // let mut w: Vec<Ciphertext> = Vec::new(); // vec of ciphertexts

        // for i in 0..l {
            // u.push(uar element)
            // y.push(DGS element)
            // let ui: Plaintext = u[i].clone()
            // w.push(self.verifier_pk.encrypt(ui, 2*y[i].clone()))
        // }

        // send the commitment

    }

    // or a tuple of two arrays
    fn reveal(response: Vec<Pow2CyclotomicPolyRing<Zq<P>, N>>) -> Vec<(Pow2CyclotomicPolyRing<Zq<P>, N>, Pow2CyclotomicPolyRing<Zq<P>, N>)> {
        todo!();
    }
}

impl<const Q: u64, const P: u64, const N: usize> Verifier<Q, P, N> {
    // instantiation
    // fn new() -> Self {
    //     Self {
    //         secret_key: Pow2CyclotomicPolyRing<Rp, N>::new(10),
    //         public_key: secret_key.generate_pk(),
    //         commitment: None
    //     }
    // }

    // fn new(secret_key: SecretKey) -> Self {
    //     Self {
    //         secret_key,
    //         public_key: secret_key.generate_pk(),
    //         commitment: None
    //     }
    // }

    fn challenge(&self, commitment: Vec<Pow2CyclotomicPolyRing<Zq<P>, N>>) -> Vec<Pow2CyclotomicPolyRing<Zq<P>, N>> {
        todo!();
        // self.commitment = Some(commitment);
    }

    fn verify(response: Vec<(Pow2CyclotomicPolyRing<Zq<P>, N>, Pow2CyclotomicPolyRing<Zq<P>, N>)>) -> bool {
        todo!();
    }
}

#[test]
fn test() {
    const N: usize = 128;
    const Q: u64 = ntt_modulus::<N>(16);
    const P: u64 = ntt_modulus::<N>(15);
    // type Rq = Zq<Q>;
    // type Rp = Zq<P>;
    // type PolyRq = Pow2CyclotomicPolyRing<Rq, N>;
    // type PolyRp = Pow2CyclotomicPolyRing<Rp, N>;

    let sk: SecretKey<Q, P, N> = SecretKey::new();
    let a = sk.poly;
    println!("{a:?}");
    let pk = sk.pk_gen();
    // let a = pk.poly1;
    // let b = pk.poly2;
    // println!("{a:?}");
    // let ptxt = Plaintext::rand_message();
    // let a = ptxt.poly;
    // println!("{a:?}");
    // let ctxt = pk.encrypt(&ptxt);

    // let act = sk.decrypt(ctxt);
    // let act = act.poly;
    // let exp = ptxt.poly;
    // assert!(act.eq(&exp));
}