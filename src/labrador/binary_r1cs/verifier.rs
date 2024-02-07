#![allow(non_snake_case)]

use nimue::ProofError;

use crate::{check, check_eq};
use crate::labrador::binary_r1cs::prover::{BinaryR1CSCRS, BinaryR1CSInstance};
use crate::labrador::binary_r1cs::util::{BinaryR1CSTranscript, reduce, SECPARAM};
use crate::labrador::verifier::verify_principal_relation;
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::traits::FromRandomBytes;
use crate::nimue::merlin::LatticeMerlin;

pub fn verify_binary_r1cs<R: PolyRing>(merlin: &mut LatticeMerlin, instance: &BinaryR1CSInstance, crs: &BinaryR1CSCRS<R>) -> Result<(), ProofError>
    where LabradorChallengeSet<R>: FromRandomBytes<R>, WeightedTernaryChallengeSet<R>: FromRandomBytes<R>
{
    //TODO: add crs and statement to transcript

    let (A, B, C) = (&instance.A, &instance.B, &instance.C);

    let d = R::dimension();
    let (k, n) = (crs.num_constraints, crs.num_variables);

    let t = merlin.next_vector::<R>(crs.m.div_ceil(d))?;

    let alpha = merlin.challenge_binary_matrix(SECPARAM, k)?;
    let beta = merlin.challenge_binary_matrix(SECPARAM, n)?;
    let gamma = merlin.challenge_binary_matrix(SECPARAM, n)?;

    // delta_i is computed mod 2, i.e., over Z2
    let delta = &alpha * A + &beta * B + &gamma * C;

    let g = merlin.next_vector::<R::BaseRing>(SECPARAM)?;

    for i in 0..g.len() {
        check_eq!(Into::<i64>::into(g[i]) % 2, 0) // Check that all g_i's are even
    }

    let transcript = BinaryR1CSTranscript { t, alpha, beta, gamma, g, delta };
    let instance_pr = reduce(&crs, &instance, &transcript);

    merlin.ratchet()?;

    verify_principal_relation(merlin, &instance_pr, &crs.pr_crs())
}