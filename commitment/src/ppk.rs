#![allow(non_snake_case)]
use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::linear_algebra::{Matrix, Vector};
// use fhe::bfv::Ciphertext;
use lattirust_arithmetic::ntt::ntt_modulus;
use lattirust_arithmetic::ring::{ConvertibleRing, PolyRing, Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT, UnsignedRepresentative, Zq};
use rand::{CryptoRng, RngCore};
use rand_distr::Normal;
use rand::distributions::{Distribution};
use ark_ff::{One, UniformRand, Zero};
use ark_std::rand::prelude::SliceRandom;
use ark_std::rand;
use lattirust_arithmetic::challenge_set::ternary;
use lattirust_arithmetic::traits::FromRandomBytes;

// TODO: toy sampling, need to use OpenFHE code 
pub fn get_gaussian_vec<T: RngCore + CryptoRng>(std_dev: f64, dimension: usize, rng: &mut T) -> Vec<i64> {
    // TODO: modulo the coefficients
    let gaussian = Normal::new(0.0, std_dev).unwrap();
    let val: Vec<i64> = (0..dimension)
        .map(|_| gaussian.sample(rng) as i64)
        .collect();
    val
}

pub fn get_poly_from_gaussian<Rng: RngCore + CryptoRng, const Q: u64, const N: usize>(std_dev: f64, dimension: usize, rng: &mut Rng) -> Pow2CyclotomicPolyRing<Zq<Q>, N> {
    let rand_vec = get_gaussian_vec(std_dev, dimension, rng);
    let rand_vec: Vec<Zq<Q>> = rand_vec.into_iter().map(|x| Zq::<Q>::from(x)).collect();

    Pow2CyclotomicPolyRing::from(rand_vec)
}


struct Prover<const Q: u64, const P: u64, const N: usize> {
    secret_key: SecretKey<Q, P, N>, // polynomial w/ coefficients in {-1, 0, +1}
    public_key: PublicKey<Q, P, N>, // ring mod q
    message: Pow2CyclotomicPolyRing<Zq<Q>, N>, // ring mod p
    verifier_pk: PublicKey<Q, P, N>
}

struct Verifier<const Q: u64, const P: u64, const N: usize> {
    secret_key: SecretKey<Q, P, N>,
    public_key: PublicKey<Q, P, N>,
    commitment: Option<Vec<Pow2CyclotomicPolyRing<Zq<Q>, N>>>
}

// TODO: do I need both modules everywhere?
// move it to a dedicated library later on
struct PublicKey<const Q: u64, const P: u64, const N: usize> {
    pk1: Pow2CyclotomicPolyRing<Zq<Q>, N>, 
    pk2: Pow2CyclotomicPolyRing<Zq<Q>, N>,
    pub modulo: u64,
}
struct SecretKey<const Q: u64, const P: u64, const N: usize> {
    sk: Pow2CyclotomicPolyRing<Zq<Q>, N>,  // p or q?
    pub modulo: u64,
}

struct Plaintext<const P: u64, const N: usize> {
    ptxt: Pow2CyclotomicPolyRing<Zq<P>, N>,
    pub modulo: u64,
}

struct Ciphertext<const Q: u64, const N: usize> {
    c1: Pow2CyclotomicPolyRing<Zq<Q>, N>,
    c2: Pow2CyclotomicPolyRing<Zq<Q>, N>,
    pub modulo: u64,
}

impl<const Q: u64, const P: u64, const N: usize> SecretKey<Q, P, N> {
    fn new() -> Self {
        let mut rng = rand::thread_rng();
        // let coeff: Matrix<Rp> = Matrix::<Rp>::rand_ternary(1, degree, &mut rng);
        // let coeff = coeff.row(0);
        // println!("{:?}", coeff);

        let bytes = Vector::<u8>::rand(N, &mut rng);
        let vec_slice = bytes.as_slice();
        // TODO: fix the type cause now the ring is not passed
        let sk = WeightedTernaryChallengeSet::<Pow2CyclotomicPolyRing<Zq<Q>, N>>::try_from_random_bytes(vec_slice).unwrap();
        Self {
            sk,
            modulo: Q,
        }
    }

    fn pk_gen(&self) -> PublicKey<Q, P, N> {
        let mut rng = rand::thread_rng();

        // e should have small std_dev, TODO: check the parameters
        let e = get_poly_from_gaussian(2.3, N, &mut rng);
        let a = Vector::<u8>::rand(N, &mut rng);
        let pk2 = Pow2CyclotomicPolyRing::<Zq<Q>, N>::try_from_random_bytes(a.as_slice()).unwrap();
        
        let pk1 = -(pk2 * self.sk) + e;
        PublicKey {
            pk1,
            pk2,
            modulo: Q,
        }
    }
    
    // fn decrypt(c: (Pow2CyclotomicPolyRing<Rp, N>, Pow2CyclotomicPolyRing<Rp, N>)) -> Pow2CyclotomicPolyRing<Rp, N> {
    //     todo!();
    // }
}

impl<const Q: u64, const P: u64, const N: usize> PublicKey<Q, P, N> {
    fn encrypt(&self, m: Plaintext<P, N>) -> Ciphertext<Q, N> {
        let mut rng = rand::thread_rng();
        let pk1 = self.pk1;
        let pk2 = self.pk2;
        
        // samples a random polynomial u in Rq (?)
        // let uniform = Uniform::new(0, Q as u8); // check
        // let val: Vec<u8> = (0..N).map(|_| uniform.sample(&mut rng)).collect();
        // let u = Pow2CyclotomicPolyRingNTT::from_fn(&val); // check
        // let bytes = Vector::<u8>::rand(N, &mut rng);
        // let u = Pow2CyclotomicPolyRing::<Zq<Q>, N>::try_from_random_bytes(bytes.as_slice()).unwrap();

        let u = Pow2CyclotomicPolyRing::<Zq<Q>, N>::rand(&mut rng);
        // samples e1, e2 in Rq
        let e1 = Pow2CyclotomicPolyRing::<Zq<Q>, N>::rand(&mut rng);
        let e2 = Pow2CyclotomicPolyRing::<Zq<Q>, N>::rand(&mut rng);
        
        let p = m.modulo;
        let q = self.modulo;

        let delta = ((q as f64) / (p as f64)).floor() as u128; // round off

        // mutliply each coeff of m with delta as bigint, and convert the polynomial into Zq
        let coeffs: Vec<Zq<P>> = m.ptxt.coeffs();
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

impl<const Q: u64, const P: u64, const N: usize> Prover<Q, P, N> {
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
