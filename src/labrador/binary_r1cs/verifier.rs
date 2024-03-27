#![allow(non_snake_case)]

use ark_relations::r1cs::ConstraintSystem;
use nimue::ProofError;

use crate::{check, check_eq};
use crate::labrador::binary_r1cs::util::{BinaryR1CSCRS, BinaryR1CSTranscript, reduce, SECPARAM, Z2};
use crate::labrador::r1cs::util::{mul_dense_sparse};
use crate::labrador::util::ark_sparse_matrices;
use crate::labrador::verifier::verify_principal_relation;
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use crate::lattice_arithmetic::poly_ring::{PolyRing, UnsignedRepresentative};
use crate::lattice_arithmetic::traits::FromRandomBytes;
use crate::nimue::lattice_merlin::LatticeMerlin;

pub fn verify_binary_r1cs<R: PolyRing>(merlin: &mut LatticeMerlin<R>, cs: &ConstraintSystem<Z2>, crs: &BinaryR1CSCRS<R>) -> Result<(), ProofError>
    where LabradorChallengeSet<R>: FromRandomBytes<R>, WeightedTernaryChallengeSet<R>: FromRandomBytes<R>
{
    //TODO: add crs and statement to transcript

    let (A, B, C) = ark_sparse_matrices(cs);

    let d = R::dimension();
    let (k, n) = (cs.num_constraints, cs.num_instance_variables + 1 + cs.num_witness_variables);

    let t = merlin.next_vector(crs.m.div_ceil(d))?;

    let alpha = merlin.challenge_binary_matrix(SECPARAM, k)?;
    let beta = merlin.challenge_binary_matrix(SECPARAM, n)?;
    let gamma = merlin.challenge_binary_matrix(SECPARAM, n)?;

    // delta_i is computed mod 2, i.e., over Z2
    let delta = mul_dense_sparse(&alpha, &A) + mul_dense_sparse(&beta, &B) + mul_dense_sparse(&gamma, &C);

    let g = merlin.next_vector_baseringelem(SECPARAM)?;

    for i in 0..g.len() {
        // Check that all g_i's are even
        check_eq!(Into::<UnsignedRepresentative>::into(g[i]).0 % 2, 0)
    }

    let transcript = BinaryR1CSTranscript { t, alpha, beta, gamma, g, delta };
    let instance_pr = reduce(&crs, &cs, &transcript);

    merlin.ratchet()?;

    verify_principal_relation(merlin, &instance_pr, &crs.core_crs)
}