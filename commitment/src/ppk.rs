#![allow(non_snake_case)]
// use fhe::bfv::Ciphertext;
use lattirust_arithmetic::ntt::ntt_modulus;
use lattirust_arithmetic::ring::{ConvertibleRing, Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT, Zq};
use rand::{CryptoRng, RngCore};
use rand_distr::Normal;
use rand::distributions::{Distribution};
use ark_ff::{One, UniformRand, Zero};
use ark_std::rand::prelude::SliceRandom;
use ark_std::rand;
use lattirust_arithmetic::challenge_set::ternary;

// todo: 
const N: usize = 128;
const Q: u64 = ntt_modulus::<N>(16);
const P: u64 = ntt_modulus::<N>(14);
type Rp = Zq<Q>;
type Rq = Zq<P>;
type PolyRp = Pow2CyclotomicPolyRing<Rp, N>;
type PolyRq = Pow2CyclotomicPolyRing<Rq, N>;

pub fn secret_ternary<T: ConvertibleRing, Rng: rand::Rng + ?Sized>(rng: &mut Rng) -> Vec<T> {
    let mut secret = Vec::new();
    for _ in 0..3 {
        secret.push([-T::one(), T::zero(), T::one()]
            .choose(rng)
            .unwrap()
            .clone());
    }
    return secret;
}

// TODO: move out
pub fn get_gaussian_vec<T: RngCore + CryptoRng>(std_dev: f64, dimension: usize, rng: &mut T) -> Vec<i64> {
    // TODO: modulo the coefficients
    let gaussian = Normal::new(0.0, std_dev).unwrap();
    let val: Vec<i64> = (0..dimension)
        .map(|_| gaussian.sample(rng) as i64)
        .collect();
    return val;
}

pub fn get_poly_from_gaussian<T: RngCore + CryptoRng, const q: u64>(std_dev: f64, dimension: usize, rng: &mut T) -> Pow2CyclotomicPolyRing<Zq<q>, N> {
    let rand_vec = get_gaussian_vec(std_dev, dimension, rng);
    let rand_vec: Vec<Zq<q>> = rand_vec.into_iter().map(|x| Zq::<q>::from(x)).collect();

    Pow2CyclotomicPolyRing::from(rand_vec)
}

struct Prover {
    secret_key: SecretKey, // polynomial w/ coefficients in {-1, 0, +1}
    public_key: PublicKey, // ring mod q
    message: Pow2CyclotomicPolyRing<Rp, N>, // ring mod p
    verifier_pk: PublicKey
}

struct Verifier {
    secret_key: SecretKey,
    public_key: PublicKey,
    commitment: Option<Vec<Pow2CyclotomicPolyRing<Rp, N>>>
}

// struct PublicKey {
//     poly1: Pow2CyclotomicPolyRing<Rp, N>,
//     poly2: Pow2CyclotomicPolyRing<Rp, N>,
// }
// struct SecretKey {
//     poly: Pow2CyclotomicPolyRing<Rp, N>,
// }

// move it to a dedicated library later on
struct PublicKey {
    pk1: Pow2CyclotomicPolyRing<Rq, N>, 
    pk2: Pow2CyclotomicPolyRing<Rq, N>,
    modulo: u64,
}
struct SecretKey {
    sk: Pow2CyclotomicPolyRing<Rp, N>,
    modulo: u64,
}

type Plaintext = Pow2CyclotomicPolyRing<Rp, N>;
type Ciphertext = (Pow2CyclotomicPolyRing<Rq, N>, Pow2CyclotomicPolyRing<Rq, N>);

impl SecretKey {
    //     let small_std_dev = 1.0; // coefficients in {-1, 0, +1}
    //     let std_dev = 5.0
    //     let mut rng = rand::rngs::StdRng::seed_from_u64(17);
    // fn new(degree: usize) -> Self {
    //         get_poly_from_gaussian(small_std_dev, degree, rng)
    // }

    // fn from_poly(poly: Pow2CyclotomicPolyRing<Rp, N>) -> Self {
    //     todo!();
    // }

    fn generate_pk(&self) -> PublicKey {
        todo!();
        // e should have small std_dev 
        // let e = get_poly_from_gaussian(small_std_dev, degree, rng);
        // let pk2 = get_poly_from_gaussian(std_dev, degree, rng);
        // let pk1 = -(pk2 * self + e);
        // (pk1, pk2)
    }
    
    fn decrypt(c: (Pow2CyclotomicPolyRing<Rp, N>, Pow2CyclotomicPolyRing<Rp, N>)) -> Pow2CyclotomicPolyRing<Rp, N> {
        todo!();
    }
}

impl PublicKey {
    fn encrypt(&self, m: Pow2CyclotomicPolyRing<Rp, N>) -> Ciphertext {
        todo!();
        // let mut rng = rand::rngs::StdRng::seed_from_u64(16);
        // let pk1 = self.pk1;
        // let pk2 = self.pk2;
        
        // // samples a random polynomial u in Rq (?)
        // // TODO: make the sampling a function
        // // let uniform = Uniform::new(0, Q as u8); // check
        // // let val: Vec<u8> = (0..N).map(|_| uniform.sample(&mut rng)).collect();
        // // let u = Pow2CyclotomicPolyRingNTT::from_fn(&val); // check
        
        // let u = PolyRq::rand(&mut rng);
        // // samples e1, e2 in Rq
        // let e1 = PolyRq::rand(&mut rng);
        // let e2 = PolyRq::rand(&mut rng);
        
        // let delta = Q/P as u64; // round off

        // let delta = Rp::from(delta);
        // // compute a, b
        // // delta * m.coeffs();
        // // TODO: convert it into NTT
        // let a = pk1 * u + e1 + m * delta;

        // // either cast m into Rq
        // // mutliply each coeff with delta as bigint, and convert the polynomial into Zq
        // // let b = pk2 * u + e2;

        // // return the ciphertext
        // // (a, b)
    }
}

impl Prover {
    // instantiation
    // fn new(verifier_pk: PublicKey, message: Pow2CyclotomicPolyRing<Rp, N>) -> Self {

        // Self {
        //     secret_key: Pow2CyclotomicPolyRing<Rp, N>,
        //     public_key: verifier_pk.generate_pk(),
        //     verifier_pk,
        //     message
        // }
    // }

    // fn new(secret_key: SecretKey, verifier_pk: PublicKey, message: Pow2CyclotomicPolyRing<Rp, N>) -> Self {
    //     Self {
    //         secret_key,
    //         public_key: secret_key.generate_pk(),
    //         verifier_pk,
    //         message
    //     }
    // }

    // generate-phase
    fn generate(&self) -> Pow2CyclotomicPolyRing<Rp, N> {
        todo!();
        // todo: use DGS
        // let mut Rp: (Pow2CyclotomicPolyRing<Rp, N>, Pow2CyclotomicPolyRing<Rp, N>, Pow2CyclotomicPolyRing<Rp, N>) = (0.0, 0.0, 0.0);
        // todo: define scalar-vector poly multiplication(?)
        // let c = self.verifier_pk.encrypt(self.message, 2 * Rp);
        
        // todo: send message
    }

    // prove-phase, which is a sigma protocol
    fn commit(&self, l: usize) -> Vec<Pow2CyclotomicPolyRing<Rp, N>> {
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
    fn reveal(response: Vec<Pow2CyclotomicPolyRing<Rp, N>>) -> Vec<(Pow2CyclotomicPolyRing<Rp, N>, Pow2CyclotomicPolyRing<Rp, N>)> {
        todo!();
    }
}

impl Verifier {
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

    fn challenge(&self, commitment: Vec<Pow2CyclotomicPolyRing<Rp, N>>) -> Vec<Pow2CyclotomicPolyRing<Rp, N>> {
        todo!();
        // self.commitment = Some(commitment);
    }

    fn verify(response: Vec<(Pow2CyclotomicPolyRing<Rp, N>, Pow2CyclotomicPolyRing<Rp, N>)>) -> bool {
        todo!();
    }
}
