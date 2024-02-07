use nimue::hash::Keccak;
use num_traits::One;

use crate::labrador::binary_r1cs::prover::{BinaryR1CSCRS, BinaryR1CSInstance, BinaryR1CSWitness, prove_binary_r1cs, Z2};
use crate::labrador::binary_r1cs::verifier::verify_binary_r1cs;
use crate::labrador::iopattern::LabradorIOPattern;
use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::ntt::ntt_modulus;
use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
use crate::nimue::iopattern::LatticeIOPattern;

const Q: u64 = ntt_modulus::<64>(32);
const D: usize = 64;

type R = Pow2CyclotomicPolyRingNTT<Q, 64>;

#[test]
#[allow(non_snake_case)]
fn test_prove_binary_r1cs() {
    let A = Matrix::<Z2>::identity(D, D);
    let B = Matrix::<Z2>::identity(D, D);
    let C = Matrix::<Z2>::identity(D, D);

    let crs = BinaryR1CSCRS::<R>::new(D, D);
    let instance = BinaryR1CSInstance { A, B, C };
    let witness = BinaryR1CSWitness { w: Vector::<Z2>::from_element(D, Z2::one()) };

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