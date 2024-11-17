use lattirust_arithmetic::ntt::ntt_modulus;
use crate::{bfv::plaintext::Plaintext, ppk::{prover::Prover, verifier::Verifier}};


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
    let response = prover.reveal(gamma.clone());
    let result = verifier.verify(response);
    assert_eq!(result, true);
}
