use log::{debug, error, info, warn};
use nimue::hash::Keccak;
use num_traits::One;
use pretty_env_logger::env_logger;

use crate::labrador::binary_r1cs::prover::prove_binary_r1cs;
use crate::labrador::binary_r1cs::util::{BinaryR1CSCRS, BinaryR1CSInstance, BinaryR1CSWitness, Z2};
use crate::labrador::binary_r1cs::verifier::verify_binary_r1cs;
use crate::labrador::iopattern::LabradorIOPattern;
use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::ntt::ntt_modulus;
use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
use crate::lattice_arithmetic::traits::Modulus;
use crate::nimue::iopattern::LatticeIOPattern;

const Q: u64 = ntt_modulus::<64>(32);
const D: usize = 64;

type R = Pow2CyclotomicPolyRingNTT<Q, 64>;

fn init() {
    let _ = env_logger::builder().is_test(true).try_init();
}

#[test]
#[allow(non_snake_case)]
fn test_prove_binary_r1cs() {
    init();

    let m: usize = (R::modulus().next_power_of_two().ilog2() * 2) as usize;
    let k = m * D;
    let A = Matrix::<Z2>::identity(k, k);
    let B = Matrix::<Z2>::identity(k, k);
    let C = Matrix::<Z2>::identity(k, k);

    let crs = BinaryR1CSCRS::<R>::new(k, k);
    let instance = BinaryR1CSInstance { A, B, C };
    let witness = BinaryR1CSWitness { w: Vector::<Z2>::from_element(k, Z2::one()) };

    let io = LatticeIOPattern::<R, Keccak>::new("labrador_binaryr1cs")
        .labrador_binaryr1cs_io(&instance, &crs)
        .ratchet()
        // .labrador_crs(&crs.pr_crs())
        // .ratchet()
        // .labrador_instance(&PrincipalRelation::<R>::new_empty(&crs.pr_crs()))
        // .ratchet()
        .labrador_io(&crs.pr_crs());
    let mut arthur = io.to_arthur();
    let proof = prove_binary_r1cs(&crs, &mut arthur, &instance, &witness).unwrap();
    println!("Finished proving, proof size = {} bytes", proof.len());

    let mut merlin = io.to_merlin(proof);
    verify_binary_r1cs(&mut merlin, &instance, &crs).unwrap();
}