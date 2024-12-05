use lattirust_arithmetic::{nimue::{iopattern, serialization::ToBytes}, ntt::ntt_modulus};
use crate::{bfv::plaintext::Plaintext, ppk::{prover::Prover, util::PublicParameters, verifier::Verifier}};

use super::{nizk::{new_ppk_io, test_io}, verifier};


const N: usize = 128;
const P: u64 = ntt_modulus::<N>(15);
const Q: u64 = P * 5;

#[test]
fn zero_message() {
    // Generate instances
    let ptxt = Plaintext::<P, N>::zero();
    let l = 6;
    let mut prover: Prover<Q, P, N> = Prover::new(ptxt, l);
    let prover_pk = prover.return_pk();
    // let pp1 = PublicParameters::<Q, P, N>::default();

    // assert_eq!(pp1.to_bytes().unwrap().len(), prover.public_parameters().to_bytes().unwrap().len());

    let mut verifier = Verifier::new(prover_pk, l);
    
    // Generate phase
    let ctxt = prover.generate();
    verifier.end_genenerate(ctxt);

    // Prove phase
    let w = prover.commit(l);
    let gamma = verifier.challenge(w, l);
    let response = prover.reveal(gamma.clone());
    let result = verifier.verify(response);
    assert_eq!(result, true);
}

#[test]
fn random_message() {
    // Generate instances
    let ptxt = Plaintext::<P, N>::rand_message();
    let l = 6;
    let mut prover: Prover<Q, P, N> = Prover::new(ptxt, l);
    let prover_pk = prover.return_pk();

    let mut verifier = Verifier::new(prover_pk, l);
    
    // Generate phase
    let ctxt = prover.generate();
    verifier.end_genenerate(ctxt);

    // Prove phase
    let w = prover.commit(l);
    let gamma = verifier.challenge(w, l);
    // let len = gamma.to_bytes().unwrap().len();
    // println!("{len}");
    let response = prover.reveal(gamma.clone());
    let result = verifier.verify(response);
    assert_eq!(result, true);
}

#[test]
fn nizk_rand_message() {
    let l = 6;
    let iopattern = new_ppk_io::<Q, P, N>(l);
    let mut arthur = iopattern.to_arthur();

    let pp = PublicParameters::<Q, P, N>::default();
    let size = pp.to_bytes().unwrap().len();
    let ptxt = Plaintext::<P, N>::rand_message();
    let mut prover: Prover<Q, P, N> = Prover::new(ptxt, l);
    let transcript  = prover.nizk_prove(&mut arthur);
    
    let mut merlin = iopattern.to_merlin(&transcript);
    let mut verifier = Verifier::<Q, P, N>::new_null();
    assert!(matches!(verifier.nizk_verify(&mut merlin, l), Ok(true)));
}

#[test]
fn random_test() {
    let iopattern = test_io::<Q, P, N>();
    let l = 6;
    let mut arthur = iopattern.to_arthur();

    let pp = PublicParameters::<Q, P, N>::default();
    let size = pp.to_bytes().unwrap().len();
    let ptxt = Plaintext::<P, N>::rand_message();
    let mut prover: Prover<Q, P, N> = Prover::new(ptxt, l);
    let transcript  = prover.test_nizk(&mut arthur);
    
    let mut merlin = iopattern.to_merlin(&transcript);
    let mut verifier = Verifier::<Q, P, N>::new_null();
    assert!(matches!(verifier.test_nizk(&mut merlin), Ok(true)));
}