#![allow(non_snake_case)]


use nimue::{Arthur, ProofError, ProofResult};
use num_traits::Zero;

use lattirust_arithmetic::balanced_decomposition::DecompositionFriendlySignedRepresentative;
use lattirust_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::nimue::arthur::SerArthur;
use lattirust_arithmetic::nimue::traits::ChallengeFromRandomBytes;
use lattirust_arithmetic::ring::representatives::WithSignedRepresentative;
use lattirust_arithmetic::ring::PolyRing;
use lattirust_arithmetic::traits::FromRandomBytes;
use lattirust_util::check;
use relations::principal_relation::{Index, Instance};
use relations::Relation;

use crate::binary_r1cs::util::{reduce, BinaryR1CSCRS, BinaryR1CSTranscript};
use crate::binary_r1cs::BinaryR1CS;
use crate::verifier::verify_principal_relation;

pub fn verify_reduction_binaryr1cs_labradorpr<R: PolyRing>(
    arthur: &mut Arthur,
    pp: &BinaryR1CSCRS<R>,
    index: &<BinaryR1CS as Relation>::Index,
    instance: &<BinaryR1CS as Relation>::Instance,
) -> ProofResult<(Index<R>, Instance<R>)>
where
    <R as PolyRing>::BaseRing: WithSignedRepresentative,
{
    let (A, B, C) = (&index.a, &index.b, &index.c);

    let (k, n) = (pp.num_constraints, pp.num_variables);

    let t = arthur.next_vector(pp.A.nrows())?;

    let alpha = arthur.challenge_binary_matrix(pp.security_parameter, k)?;
    let beta = arthur.challenge_binary_matrix(pp.security_parameter, n)?;
    let gamma = arthur.challenge_binary_matrix(pp.security_parameter, n)?;

    // delta_i is computed mod 2, i.e., over Z2
    let delta = &alpha * A + &beta * B + &gamma * C;

    let g = arthur.next_vector_canonical::<R::BaseRing>(pp.security_parameter)?;

    for g_i in &g {
        // Check that all g_i's are even
        let two = R::BaseRing::try_from(2u128).unwrap().as_signed_representative();
        check!(
            (g_i.as_signed_representative() % two).is_zero()
        );
    }

    // TODO: ratchet to make sure we consumed everything?
    let transcript = BinaryR1CSTranscript {
        t,
        alpha,
        beta,
        gamma,
        g,
        delta,
    };
    
    let instance_pr = reduce(pp, &transcript);
    Ok((pp.core_crs.clone(), instance_pr))
}

pub fn verify_binary_r1cs<R: PolyRing>(
    arthur: &mut Arthur,
    pp: &BinaryR1CSCRS<R>,
    index: &<BinaryR1CS as Relation>::Index,
    instance: &<BinaryR1CS as Relation>::Instance,
) -> Result<(), ProofError>
where
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
    <R as PolyRing>::BaseRing: WithSignedRepresentative,
    <<R as PolyRing>::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,


{
    //TODO: add crs and statement to transcript
    let (index_pr, instance_pr) =
        verify_reduction_binaryr1cs_labradorpr(arthur, pp, index, instance)?;
    
    debug_assert_eq!(index_pr, pp.core_crs);

    arthur.ratchet()?;

    verify_principal_relation(arthur, &index_pr, &instance_pr)
}
