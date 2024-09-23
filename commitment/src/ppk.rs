#![allow(non_snake_case)]
use lattirust_arithmetic::ntt::ntt_modulus;
use lattirust_arithmetic::ring::{PolyRing, Pow2CyclotomicPolyRing, Zq};
use rand::{CryptoRng, RngCore};
use rand::{random, rngs, Rng, SeedableRng, thread_rng, rngs::OsRng};
use rand_distr::Normal;
use rand::distributions::{Distribution, Alphanumeric, Uniform, Standard};
// use fhe::bfv::{BfvParametersBuilder, Ciphertext, Encoding, Plaintext, PublicKey, SecretKey};
// use fhe_traits;
use std::error::Error;

// todo: 
const N: usize = 128;
const Q: u64 = ntt_modulus::<N>(16);
const P: u64 = ntt_modulus::<N>(16);
type R = Zq<Q>;

// from cathieyun/bfv12
pub fn get_gaussian_vec<T: RngCore + CryptoRng>(std_dev: f64, dimension: usize, rng: &mut T) -> Vec<i64> {
    let gaussian = Normal::new(0.0, std_dev).unwrap();
    let val: Vec<i64> = (0..dimension)
        .map(|_| gaussian.sample(rng) as i64)
        .collect();
    return val;
}

// pub fn get_poly_from_gaussian<T: RngCore + CryptoRng>(std_dev: f64, dimension: usize, rng: &mut T) -> Pow2CyclotomicPolyRing<R, N> {
//     let v = get_gaussian_vec(std_dev, dimension, &mut rng);

// }

struct Prover {
    secret_key: SecretKey, // polynomial w/ coefficients in {-1, 0, +1}
    public_key: PublicKey, // ring mod q
    message: Pow2CyclotomicPolyRing<R, N>, // ring mod p
    verifier_pk: PublicKey
}

struct Verifier {
    secret_key: SecretKey,
    public_key: PublicKey,
    commitment: Option<Vec<Pow2CyclotomicPolyRing<R, N>>>
}

// struct PublicKey {
//     poly1: Pow2CyclotomicPolyRing<R, N>,
//     poly2: Pow2CyclotomicPolyRing<R, N>,
// }
// struct SecretKey {
//     poly: Pow2CyclotomicPolyRing<R, N>,
// }

// move it to a dedicated library later on
struct PublicKey(Pow2CyclotomicPolyRing<R, N>, Pow2CyclotomicPolyRing<R, N>);
struct SecretKey(Pow2CyclotomicPolyRing<R, N>);

type Plaintext = Pow2CyclotomicPolyRing<R, N>;
type Ciphertext = (Pow2CyclotomicPolyRing<R, N>, Pow2CyclotomicPolyRing<R, N>);

impl SecretKey {
    // fn new(degree: usize) -> Self {
    //     let std_dev = 1; // coefficients in {-1, 0, +1}
    //     let mut rng = rand::thread_rng();

        
    // }

    fn from_poly(poly: Pow2CyclotomicPolyRing<R, N>) -> Self {
        todo!();
    }

    fn generate_pk() -> PublicKey {
        todo!();
    }
    
    fn decrypt(c: (Pow2CyclotomicPolyRing<R, N>, Pow2CyclotomicPolyRing<R, N>)) -> Pow2CyclotomicPolyRing<R, N> {
        todo!();
    }
}

impl PublicKey {
    fn encrypt(m: Pow2CyclotomicPolyRing<R, N>, r: (Pow2CyclotomicPolyRing<R, N>, Pow2CyclotomicPolyRing<R, N>, Pow2CyclotomicPolyRing<R, N>)) -> Pow2CyclotomicPolyRing<R, N> {
        todo!();
    }
}

impl Prover {
    // instantiation
    // fn new(verifier_pk: PublicKey, message: Pow2CyclotomicPolyRing<R, N>) -> Self {

        // Self {
        //     secret_key: Pow2CyclotomicPolyRing<R, N>,
        //     public_key: verifier_pk.generate_pk(),
        //     verifier_pk,
        //     message
        // }
    // }

    // fn new(secret_key: SecretKey, verifier_pk: PublicKey, message: Pow2CyclotomicPolyRing<R, N>) -> Self {
    //     Self {
    //         secret_key,
    //         public_key: secret_key.generate_pk(),
    //         verifier_pk,
    //         message
    //     }
    // }

    // generate-phase
    fn generate(&self) -> Pow2CyclotomicPolyRing<R, N> {
        todo!();
        // todo: use DGS
        // let mut r: (Pow2CyclotomicPolyRing<R, N>, Pow2CyclotomicPolyRing<R, N>, Pow2CyclotomicPolyRing<R, N>) = (0.0, 0.0, 0.0);
        // todo: define scalar-vector poly multiplication(?)
        // let c = self.verifier_pk.encrypt(self.message, 2 * r);
        
        // todo: send message
    }

    // prove-phase, which is a sigma protocol
    fn commit(&self, l: usize) -> Vec<Pow2CyclotomicPolyRing<R, N>> {
        todo!();

        // let mut u: Vec<Pow2CyclotomicPolyRing<R, N>> = Vec::new(); // vec poly of size l
        // let mut y: Vec<Pow2CyclotomicPolyRing<R, N>> = Vec::new(); // vec of 3-tuple poly of size l
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
    fn reveal(response: Vec<Pow2CyclotomicPolyRing<R, N>>) -> Vec<(Pow2CyclotomicPolyRing<R, N>, Pow2CyclotomicPolyRing<R, N>)> {
        todo!();
    }
}

impl Verifier {
    // instantiation
    // fn new() -> Self {
    //     Self {
    //         secret_key: Pow2CyclotomicPolyRing<R, N>::new(10),
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

    fn challenge(&self, commitment: Vec<Pow2CyclotomicPolyRing<R, N>>) -> Vec<Pow2CyclotomicPolyRing<R, N>> {
        todo!();
        // self.commitment = Some(commitment);
    }

    fn verify(response: Vec<(Pow2CyclotomicPolyRing<R, N>, Pow2CyclotomicPolyRing<R, N>)>) -> bool {
        todo!();
    }
}
