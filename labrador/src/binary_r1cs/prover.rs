#![allow(non_snake_case)]

use ark_relations::r1cs::ConstraintSystem;
use log::debug;
use nimue::{Arthur, ProofResult};

use lattirust_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::nimue::arthur::SerArthur;
use lattirust_arithmetic::nimue::traits::ChallengeFromRandomBytes;
use lattirust_arithmetic::ring::PolyRing;
use lattirust_arithmetic::traits::FromRandomBytes;
use relations::principal_relation::Witness;

use crate::binary_r1cs::util::{reduce, BinaryR1CSCRS, BinaryR1CSTranscript, Z2};
use crate::prover::prove_principal_relation;
use crate::r1cs::util::{embed, lift};
use crate::util::{ark_sparse_matrices, concat};

pub fn prove_binary_r1cs<'a, R: PolyRing>(
    crs: &BinaryR1CSCRS<R>,
    arthur: &'a mut Arthur,
    mut cs: &ConstraintSystem<Z2>,
) -> ProofResult<&'a [u8]>
where
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
{
    debug!("labrador::binary_r1cs starting BinR1CS -> PrincipalRelation reduction");

    //cs.set_mode(ark_relations::r1cs::SynthesisMode::Prove { construct_matrices: true });
    let (A, B, C) = ark_sparse_matrices(cs);
    let w = concat(
        vec![
            cs.instance_assignment.as_slice(),
            cs.witness_assignment.as_slice(),
        ]
        .as_slice(),
    );
    let (k, n) = (
        cs.num_constraints,
        cs.num_instance_variables + cs.num_witness_variables,
    );
    assert_eq!(k, n, "the current implementation only support k = n"); // TODO: remove this restriction by splitting a,b,c or w into multiple vectors

    // TODO: add statement

    let a = &A * &w;
    let b = &B * &w;
    let c = &C * &w;

    debug!("labrador::binary_r1cs computed A*w, B*w, C*w");

    let a_R = lift::<R>(&a);
    let b_R = lift::<R>(&b);
    let c_R = lift::<R>(&c);
    let w_R = lift::<R>(&w);

    let v = concat(&[
        a_R.as_slice(),
        b_R.as_slice(),
        c_R.as_slice(),
        w_R.as_slice(),
    ]);
    let t = &crs.A * &v;

    arthur.absorb_vector(&t).unwrap();

    debug!(
        "labrador::binary_r1cs squeeze alpha in {{0,1}}^{}x{k}",
        crs.security_parameter
    );
    let alpha = arthur
        .challenge_binary_matrix(crs.security_parameter, k)
        .unwrap();
    debug!(
        "labrador::binary_r1cs squeeze beta in {{0,1}}^{}x{k}",
        crs.security_parameter
    );
    let beta = arthur
        .challenge_binary_matrix(crs.security_parameter, k)
        .unwrap();
    debug!(
        "labrador::binary_r1cs squeeze gamma in {{0,1}}^{}x{k}",
        crs.security_parameter
    );
    let gamma = arthur
        .challenge_binary_matrix(crs.security_parameter, k)
        .unwrap();

    // delta_i is computed mod 2, i.e., over Z2
    let delta = &alpha * &A + &beta * &B + &gamma * &C;

    // g_i is computed over Zq
    let (a_BR, b_BR, c_BR) = (
        a.map(embed::<R::BaseRing>),
        b.map(embed::<R::BaseRing>),
        c.map(embed::<R::BaseRing>),
    );
    let (alpha_BR, beta_BR, gamma_BR) = (
        alpha.map(embed::<R::BaseRing>),
        beta.map(embed::<R::BaseRing>),
        gamma.map(embed::<R::BaseRing>),
    );
    let (delta_BR, w_BR) = (delta.map(embed::<R::BaseRing>), w.map(embed::<R::BaseRing>));

    let g = &alpha_BR * &a_BR + &beta_BR * &b_BR + &gamma_BR * &c_BR - &delta_BR * &w_BR;

    debug_assert_eq!(
        g.len(),
        crs.security_parameter,
        "g has length {} but should have length {}",
        g.len(),
        crs.security_parameter
    );
    debug!(
        "labrador::binary_r1cs absorb g in R_q^{}",
        crs.security_parameter
    );
    arthur.absorb_vector_canonical::<R::BaseRing>(&g).unwrap();

    let a_tilde = R::sigma_vec(&a_R);
    let b_tilde = R::sigma_vec(&b_R);
    let c_tilde = R::sigma_vec(&c_R);
    let w_tilde = R::sigma_vec(&w_R);

    let transcript = BinaryR1CSTranscript {
        t,
        alpha,
        beta,
        gamma,
        g,
        delta,
    };

    let instance_pr = reduce(&crs, &cs, &transcript);

    let witness_pr = Witness::<R> {
        s: vec![a_R, b_R, c_R, w_R, a_tilde, b_tilde, c_tilde, w_tilde], // see definition of indices above
    };

    arthur.ratchet()?;

    debug!("labrador::binary_r1cs finished BinR1CS -> PrincipalRelation reduction");
    prove_principal_relation(arthur, &instance_pr, &witness_pr, &crs.core_crs)
}
