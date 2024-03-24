#![allow(non_snake_case)]

use anyhow::Error;
use ark_relations::r1cs::ConstraintSystem;
use log::debug;
use num_traits::One;

use crate::labrador::binary_r1cs::util::{BinaryR1CSCRS, BinaryR1CSTranscript, reduce, SECPARAM, Z2};
use crate::labrador::prover::prove_principal_relation;
use crate::labrador::r1cs::util::{ark_sparse_matrices, embed, is_satisfied, is_wellformed, lift, mul_dense_sparse};
use crate::labrador::util::concat;
use crate::labrador::witness::Witness;
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::traits::FromRandomBytes;
use crate::nimue::lattice_arthur::LatticeArthur;

pub fn prove_binary_r1cs<'a, R: PolyRing>(crs: &BinaryR1CSCRS<R>, arthur: &'a mut LatticeArthur<R>, mut cs: &ConstraintSystem<Z2>) -> Result<&'a [u8], Error>
    where LabradorChallengeSet<R>: FromRandomBytes<R>, WeightedTernaryChallengeSet<R>: FromRandomBytes<R>
{
    debug!("labrador::binary_r1cs starting BinR1CS -> PrincipalRelation reduction");
    // debug_assert!(is_wellformed(crs, cs));
    // debug_assert!(is_satisfied(crs, cs));

    //cs.set_mode(ark_relations::r1cs::SynthesisMode::Prove { construct_matrices: true });
    let (A, B, C) = ark_sparse_matrices(cs);
    let w = concat(vec![cs.instance_assignment.as_slice(), cs.witness_assignment.as_slice()].as_slice());
    let (k, n) = (cs.num_constraints, cs.num_instance_variables + cs.num_witness_variables);
    assert_eq!(k, n, "the current implementation only support k = n"); // TODO: remove this restriction by splitting a,b,c or w into multiple vectors

    // TODO: add statement

    let a = &A * &w;
    let b = &B * &w;
    let c = &C * &w;

    debug!("labrador::binary_r1cs computed A*w, B*w, C*w");

    let a_R = lift::<R>(a.as_slice());
    let b_R = lift::<R>(b.as_slice());
    let c_R = lift::<R>(c.as_slice());
    let w_R = lift::<R>(w.as_slice());

    let v = concat(&[a_R.as_slice(), b_R.as_slice(), c_R.as_slice(), w_R.as_slice()]);
    let t = &crs.A * &v;

    arthur.absorb_vector(&t).unwrap();

    debug!("labrador::binary_r1cs squeeze alpha in {{0,1}}^{SECPARAM}x{k}");
    let alpha = arthur.squeeze_binary_matrix(SECPARAM, k).unwrap();
    debug!("labrador::binary_r1cs squeeze beta in {{0,1}}^{SECPARAM}x{k}");
    let beta = arthur.squeeze_binary_matrix(SECPARAM, k).unwrap();
    debug!("labrador::binary_r1cs squeeze gamma in {{0,1}}^{SECPARAM}x{k}");
    let gamma = arthur.squeeze_binary_matrix(SECPARAM, k).unwrap();

    // delta_i is computed mod 2, i.e., over Z2
    let delta = mul_dense_sparse(&alpha, &A) + mul_dense_sparse(&beta, &B) + mul_dense_sparse(&gamma, &C);

    // g_i is computed over Zq
    let (a_BR, b_BR, c_BR) = (a.map(embed::<R::BaseRing>), b.map(embed::<R::BaseRing>), c.map(embed::<R::BaseRing>));
    let (alpha_BR, beta_BR, gamma_BR) = (alpha.map(embed::<R::BaseRing>), beta.map(embed::<R::BaseRing>), gamma.map(embed::<R::BaseRing>));
    let (delta_BR, w_BR) = (delta.map(embed::<R::BaseRing>), w.map(embed::<R::BaseRing>));

    let g = &alpha_BR * &a_BR + &beta_BR * &b_BR + &gamma_BR * &c_BR - &delta_BR * &w_BR;

    debug_assert_eq!(g.len(), SECPARAM, "g has length {} but should have length {}", g.len(), SECPARAM);
    debug!("labrador::binary_r1cs absorb g in R_q^{SECPARAM}");
    arthur.absorb_vector_baseringlem(&g).unwrap();

    let a_tilde = R::sigma_vec(&a_R);
    let b_tilde = R::sigma_vec(&b_R);
    let c_tilde = R::sigma_vec(&c_R);
    let w_tilde = R::sigma_vec(&w_R);

    let transcript = BinaryR1CSTranscript { t, alpha, beta, gamma, g, delta };

    let instance_pr = reduce(&crs, &cs, &transcript);

    let witness_pr = Witness::<R> {
        s: vec![a_R, b_R, c_R, w_R, a_tilde, b_tilde, c_tilde, w_tilde] // see definition of indices above
    };

    arthur.ratchet()?;

    debug!("labrador::binary_r1cs finished BinR1CS -> PrincipalRelation reduction");
    prove_principal_relation(arthur, &instance_pr, &witness_pr, &crs.core_crs)
}

