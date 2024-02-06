use nimue::hash::Keccak;

use crate::labrador::iopattern::LabradorIOPattern;
use crate::labrador::prover::{prove_principal_relation, Witness};
use crate::labrador::setup::setup;
use crate::labrador::verifier::verify_principal_relation;
use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::ntt::ntt_modulus;
use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
use crate::lattice_arithmetic::ring::{Ring, Zq};
use crate::lattice_arithmetic::traits::Modulus;
use crate::nimue::iopattern::LatticeIOPattern;
use crate::relations::labrador::principal_relation::PrincipalRelation;

const D: usize = 64;
// R1CS over Z_{2^D+1}
const Q: u64 = ntt_modulus::<{ 64 }>(32);

// Example parameters from Table 3 in the Labrador paper
const N: usize = 325;
const R: usize = 5;

type R = Zq<Q>;

type PolyR = Pow2CyclotomicPolyRingNTT<Q, { 64 }>;

fn get_beta<R: Ring + Modulus>(num_r1cs_constraints: usize, num_r1cs_variables: usize, d: usize) -> f64 {
    assert_eq!(d, 64); // Otherwise we need to recompute l correctly
    let l: usize = 18;// l = round(128 / lp), where lp = log2(p) is the bitsize of the smallest prime factor of 2^D+1, which is 18 bits for D=64

    let temp = (num_r1cs_variables + (3 + l) * num_r1cs_constraints) * D;
    assert!(temp < (0.3 * (Q as f64)) as usize);

    let beta = f64::sqrt(128f64 / 30f64) * f64::sqrt(temp as f64 / 2f64);
    let val: f64 = (temp as f64) * beta + (beta * beta) / 2f64;
    assert!(val < R::modulus() as f64);
    beta
}

#[test]
fn test_setup() {
    let b = R::from(7u64);
    let num_r1cs_constraints: usize = 10;
    let num_r1cs_variables: usize = 8;

    let beta = get_beta::<R>(num_r1cs_constraints, num_r1cs_variables, D);

    let _crs = setup::<PolyR>(10, N, D, beta, 10, 10, 10, 2, 3, b);
}

#[test]
fn test_prove() {
    let b = R::from(7u64);
    let num_r1cs_constraints: usize = 10;
    let num_r1cs_variables: usize = 8;
    let num_quad_constraints: usize = 5; //TODO
    let num_ct_quad_constraints: usize = 7; //TODO

    let beta = get_beta::<R>(num_r1cs_constraints, num_r1cs_variables, D);

    let crs = setup::<PolyR>(R, N, D, beta, 10, 10, 10, num_quad_constraints, num_ct_quad_constraints, b);
    let instance = PrincipalRelation::<PolyR>::new_dummy(R, N, beta, num_quad_constraints, num_ct_quad_constraints);
    let witness = Witness::<PolyR> { s: vec![Vector::<PolyR>::zeros(N); R] };

    let io = LatticeIOPattern::<PolyR, Keccak>::new("labrador")
        .labrador_instance(&instance)
        .ratchet()
        .labrador_io(&crs);
    let mut arthur = io.to_arthur();
    let proof = prove_principal_relation(&mut arthur, &instance, &witness, &crs);
    assert!(proof.is_ok());
}

#[test]
fn test_verify() {
    let b = R::from(7u64);
    let num_r1cs_constraints: usize = 10;
    let num_r1cs_variables: usize = 8;
    let num_quad_constraints: usize = 5; //TODO
    let num_ct_quad_constraints: usize = 7; //TODO

    let beta = get_beta::<R>(num_r1cs_constraints, num_r1cs_variables, D);

    let crs = setup::<PolyR>(R, N, D, beta, 10, 10, 10, num_quad_constraints, num_ct_quad_constraints, b);
    let instance = PrincipalRelation::<PolyR>::new_dummy(R, N, beta, num_quad_constraints, num_ct_quad_constraints);
    let witness = Witness::<PolyR> { s: vec![Vector::<PolyR>::zeros(N); R] };

    let io = LatticeIOPattern::<PolyR, Keccak>::new("labrador")
        .labrador_crs(&crs)
        .ratchet()
        .labrador_instance(&instance)
        .ratchet()
        .labrador_io(&crs);
    let mut arthur = io.to_arthur();
    let proof = prove_principal_relation(&mut arthur, &instance, &witness, &crs);

    let mut merlin = io.to_merlin(proof.unwrap());
    let result = verify_principal_relation(&mut merlin, &instance, &crs);
    assert!(result.is_ok(), "{:?}", result);
}