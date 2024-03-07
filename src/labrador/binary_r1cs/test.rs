use ark_relations::lc;
use log::{debug, error, info, warn};
use nimue::hash::Keccak;
use num_traits::One;
use pretty_env_logger::env_logger;

use crate::labrador::binary_r1cs::prover::prove_binary_r1cs;
use crate::labrador::binary_r1cs::util::{BinaryR1CSCRS, Z2};
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

    let m: usize = (1 << 15); //(R::modulus().next_power_of_two().ilog2() * 2) as usize;
    let k = m * D;

    debug!("Constructing constraint system...");
    let mut cs = ark_relations::r1cs::ConstraintSystem::<Z2>::new();
    for _ in 0..k {
        let v = cs.new_witness_variable(|| Ok(Z2::one())).unwrap();
        cs.enforce_constraint(lc!() + v, lc!() - v, lc!() + v).unwrap();
    }
    let mats = cs.to_matrices().unwrap();
    debug!("\t{} ≈ 2^{} constraints, {}+{} ≈ 2^{} variables, ({}, {}, {}) non-zero entries", cs.num_constraints,cs.num_constraints.next_power_of_two().ilog2(), cs.num_instance_variables, cs.num_witness_variables, (cs.num_instance_variables+cs.num_witness_variables).next_power_of_two().ilog2(), mats.a_num_non_zero, mats.b_num_non_zero, mats.c_num_non_zero);

    debug!("Constructing common reference string...");
    let crs = BinaryR1CSCRS::<R>::new(k, k);

    debug!("Setting IOPattern...");
    let io = LatticeIOPattern::<R, Keccak>::new("labrador_binaryr1cs")
        .labrador_binaryr1cs_io(&cs, &crs)
        .ratchet()
        // .labrador_crs(&crs.pr_crs())
        // .ratchet()
        // .labrador_instance(&PrincipalRelation::<R>::new_empty(&crs.pr_crs()))
        // .ratchet()
        .labrador_io(&crs.core_crs);
    let mut arthur = io.to_arthur();

    debug!("Proving...");
    let proof = prove_binary_r1cs(&crs, &mut arthur, &cs).unwrap();
    debug!("Finished proving, proof size = {} bytes", proof.len());

    let mut merlin = io.to_merlin(proof);
    verify_binary_r1cs(&mut merlin, &cs, &crs).unwrap();
}