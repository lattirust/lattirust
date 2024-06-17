#![allow(non_snake_case)]

use ark_relations::r1cs::ConstraintSystem;
use log::debug;
use nimue::{Merlin, ProofError};

use lattirust_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::nimue::merlin::SerMerlin;
use lattirust_arithmetic::nimue::traits::ChallengeFromRandomBytes;
use lattirust_arithmetic::ring::{PolyRing, UnsignedRepresentative};
use lattirust_arithmetic::traits::FromRandomBytes;
use lattirust_util::{check, check_eq};

use crate::binary_r1cs::util::{reduce, BinaryR1CSCRS, BinaryR1CSTranscript, Z2};
use crate::util::ark_sparse_matrices;
use crate::verifier::verify_principal_relation;

pub fn verify_binary_r1cs<R: PolyRing>(
    merlin: &mut Merlin,
    cs: &ConstraintSystem<Z2>,
    crs: &BinaryR1CSCRS<R>,
) -> Result<(), ProofError>
where
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
{
    //TODO: add crs and statement to transcript

    let (A, B, C) = ark_sparse_matrices(cs);

    let d = R::dimension();
    let (k, n) = (
        cs.num_constraints,
        cs.num_instance_variables + 1 + cs.num_witness_variables,
    );

    let t = merlin.next_vector(crs.m.div_ceil(d))?;

    let alpha = merlin.challenge_binary_matrix(crs.security_parameter, k)?;
    let beta = merlin.challenge_binary_matrix(crs.security_parameter, n)?;
    let gamma = merlin.challenge_binary_matrix(crs.security_parameter, n)?;

    // delta_i is computed mod 2, i.e., over Z2
    let delta = &alpha * &A + &beta * &B + &gamma * &C;

    let g = merlin.next_vector_canonical::<R::BaseRing>(crs.security_parameter)?;

    for i in 0..g.len() {
        // Check that all g_i's are even
        check_eq!(Into::<UnsignedRepresentative>::into(g[i]).0 % 2, 0);
    }

    let transcript = BinaryR1CSTranscript {
        t,
        alpha,
        beta,
        gamma,
        g,
        delta,
    };
    let instance_pr = reduce(&crs, &cs, &transcript);

    merlin.ratchet()?;

    verify_principal_relation(merlin, &instance_pr, &crs.core_crs)
}
