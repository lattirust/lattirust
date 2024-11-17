#![allow(non_snake_case)]

use lattirust_arithmetic::ring::Zq;
use ark_ff::{UniformRand, Zero};
use ark_std::rand;
use nimue::Arthur;
use crate::bfv::{ciphertext::{self, Ciphertext}, plaintext::Plaintext, public_key::PublicKey, secret_key::SecretKey, util::{convert_ring, rand_tuple, PolyR, TuplePolyR}};

use super::{shared::PPKIOPattern, util::PublicParameters};

pub struct Prover<const Q: u64, const P: u64, const N: usize> {
    // params: ParamsBFV,
    _secret_key: SecretKey<Q, P, N>, // polynomial w/ coefficients in {-1, 0, +1}
    public_key: PublicKey<Q, P, N>, // ring mod q
    message: Plaintext<P, N>, // ring mod p
    r: TuplePolyR<Q, N>,
    u: Vec<PolyR<P, N>>,
    y: Vec<TuplePolyR<Q, N>>,
    l: usize,
}

impl<const Q: u64, const P: u64, const N: usize> Prover<Q, P, N> {
    // instantiation
    // TODO: cannot use verifier pk to encrypt, would leak m
    pub fn new(message: Plaintext<P, N>, l: usize) -> Self {
        let _secret_key =  SecretKey::new(); 
        let public_key = _secret_key.pk_gen();
        
        Self {
            _secret_key, 
            public_key, 
            message,
            r: TuplePolyR::<Q, N>::default(),
            u: Vec::default(),
            y: Vec::default(),
            l
        }
    }

    // TODO: move PublicKey out and derive clone
    pub fn return_pk(&self) -> PublicKey<Q, P, N> {
        PublicKey {
            poly1: self.public_key.poly1.clone(),
            poly2: self.public_key.poly2.clone(),
            modulo: Q,
        }
    }

    // generate-phase
    // TODO: use RefCell and make self immutable (?)
    pub fn generate(&mut self) -> Ciphertext<Q, N> {
        let two = PolyR::<Q, N>::from(Zq::<Q>::from(2));

        // todo: use DGS
        let r = rand_tuple::<Q, N>(None);
        self.r = r.clone();
        let r = r * two;

        self.public_key.encrypt(&self.message, r)
    }

    // prove-phase, which is a sigma protocol with soundness amplification
    pub fn commit(&mut self, l: usize) -> Vec<Ciphertext<Q, N>> {
        let two = PolyR::<Q, N>::from(Zq::<Q>::from(2));
        let mut rng = rand::thread_rng();
        let mut u: Vec<PolyR<P, N>> = Vec::new(); // vec poly of size l
        let mut y_times_2: Vec<TuplePolyR<Q, N>> = Vec::new(); 
        let mut w: Vec<Ciphertext<Q, N>> = Vec::new(); // vec of ciphertexts

        for i in 0..l {
            u.push(PolyR::<P, N>::rand(&mut rng));

            let yi: TuplePolyR<Q, N> = rand_tuple::<Q, N>(None);
            self.y.push(yi.clone());
            y_times_2.push(yi.clone() * two);
            
            // transform ui into a Plaintext so we can encrypt it
            let ui = Plaintext::<P, N> {
                poly: u[i].clone(),
                modulo: P,
            };
            // TODO: create a helper method to use self.encrypt instead of self.pk.encrypt
            w.push(self.public_key.encrypt(&ui, y_times_2[i].clone())) 
        }

        self.u = u;
        
        w
    }

    pub fn reveal(&self, challenge: Vec<PolyR<P, N>>) -> (Vec<PolyR<P, N>>, Vec<TuplePolyR<Q, N>>) {
        let mut v: Vec<PolyR<P, N>> = Vec::new();
        let mut z: Vec<TuplePolyR<Q, N>> = Vec::new();
        let l = self.l;
        let m = self.message.poly.clone();
        let r = self.r.clone();

        for i in 0..l {
            let u_i = self.u[i].clone();
            let y_i = self.y[i].clone();
            let gamma_i = challenge[i].clone();

            let v_i: PolyR<P, N> = u_i + gamma_i * m;

            let gamma_i: PolyR<Q, N> = convert_ring::<P, Q, N>(gamma_i);
            let z_i = TuplePolyR(
                y_i.0.clone() + gamma_i.clone() * r.0.clone(), 
                y_i.1.clone() + gamma_i.clone() * r.1.clone(),
                y_i.2.clone() + gamma_i.clone() * r.2.clone()
            );

            v.push(v_i);
            z.push(z_i);
        }


        (v, z)
    }

    // pub fn nizk(&self) {
    //     let io = PPKIOPattern::<Q, P, N>::generate(1, 2, 3, 4);
    
    //     // instantiate a prover
    //     let mut arthur = io.to_arthur();
    
    //     // prepare all the ingredients
    //     let l = 6;
        

    // }
}

// pub fn prove_relation<'a, R: PolyRing> (
//     arthur: &'a mut Arthur,
//     pp: &PublicParameters<R>,
//     witness: &Witness
//     witness: &Witness<R>,
// ) -> ProofResult<&'a [u8]>
// where 
//     PPKChallengeSet<R>: FromRandomBytes<R>,
// {
//     // let num_constraints = instance.quad_dot_prod_funcs.len();
//     // let num_ct_constraints = instance.ct_quad_dot_prod_funcs.len();
//     // let mut prover = Prover::new(message, l);
    
//     arthur.absorb_vector()
// }