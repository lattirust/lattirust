use nimue::hash::Keccak;

use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
use crate::lattice_arithmetic::ring::{Ring, Zq};
use crate::lattice_arithmetic::traits::Modulus;
use crate::nimue::iopattern::LatticeIOPattern;
use crate::relations::labrador::principal_relation::PrincipalRelation;
use crate::labrador::iopattern::LabradorIOPattern;
use crate::labrador::prover::{prove_principal_relation, Witness};
use crate::labrador::setup::setup;
use crate::labrador::verifier::verify_principal_relation;

const Q: u64 = 4294967291;
// 2^32-5, prime
const D: usize = 64;
// R1CS over Z_{2^D+1}

const N: usize = 1;

type R = Zq<Q>;

type PolyR = Pow2CyclotomicPolyRing<R, D>;

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

    let _crs = setup::<PolyR>(10, N, D, beta, 10, 10, 10, b);
}

#[test]
fn test_prove() {
    let r = 3;
    let b = R::from(7u64);
    let num_r1cs_constraints: usize = 10;
    let num_r1cs_variables: usize = 8;

    let beta = get_beta::<R>(num_r1cs_constraints, num_r1cs_variables, D);

    let crs = setup::<PolyR>(r, N, D, beta, 10, 10, 10, b);
    let instance = PrincipalRelation::<PolyR>::new_dummy(r, N, beta, 2, 3);
    let witness = Witness::<PolyR> { s: vec![Vector::<PolyR>::zeros(N); r] };

    let io = LatticeIOPattern::<PolyR, Keccak>::new("labrador")
        .labrador_statement()
        .ratchet()
        .labrador_io(&instance, &crs);
    let mut arthur = io.to_arthur();
    let proof = prove_principal_relation(&mut arthur, instance, witness, crs);
    assert!(proof.is_ok());
}

#[test]
fn test_verify() {
    let r = 3;
    let b = R::from(7u64);
    let num_r1cs_constraints: usize = 10;
    let num_r1cs_variables: usize = 8;

    let beta = get_beta::<R>(num_r1cs_constraints, num_r1cs_variables, D);

    let crs = setup::<PolyR>(r, N, D, beta, 10, 10, 10, b);
    let instance = PrincipalRelation::<PolyR>::new_dummy(r, N, beta, 2, 3);
    let witness = Witness::<PolyR> { s: vec![Vector::<PolyR>::zeros(N); r] };

    let io = LatticeIOPattern::<PolyR, Keccak>::new("labrador")
        .labrador_statement()
        .ratchet()
        .labrador_io(&instance, &crs);
    let mut arthur = io.to_arthur();
    let proof = prove_principal_relation(&mut arthur, instance.clone(), witness, crs.clone());

    let mut merlin = io.to_merlin(proof.unwrap());
    let result = verify_principal_relation(&mut merlin, instance.clone(), crs.clone());
    assert!(result.is_ok());
}