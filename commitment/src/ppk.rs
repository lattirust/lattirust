#![allow(non_snake_case)]
use lattirust_arithmetic::ntt::ntt_modulus;
use lattirust_arithmetic::ring::{Pow2CyclotomicPolyRing, Zq};
use rand::{random, rngs, Rng, SeedableRng};
use rand::distributions::{Distribution, Alphanumeric, Uniform, Standard};

// todo: 
const N: usize = 128;
const Q: u64 = ntt_modulus::<N>(16);
type R = Zq<Q>;

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

struct PublicKey(Pow2CyclotomicPolyRing<R, N>, Pow2CyclotomicPolyRing<R, N>);
struct SecretKey(Pow2CyclotomicPolyRing<R, N>);

type Plaintext = Pow2CyclotomicPolyRing<R, N>;
type Ciphertext = (Pow2CyclotomicPolyRing<R, N>, Pow2CyclotomicPolyRing<R, N>);

impl SecretKey {
    fn new(degree: usize) -> Self {
        todo!();
    }

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