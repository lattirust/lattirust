#![allow(non_snake_case)]

use log::debug;
use nimue::{Arthur, ProofError, ProofResult};
use num_traits::ToPrimitive;

use lattirust_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::decomposition::DecompositionFriendlySignedRepresentative;
use lattirust_arithmetic::nimue::arthur::SerArthur;
use lattirust_arithmetic::nimue::traits::ChallengeFromRandomBytes;
use lattirust_arithmetic::ring::representatives::WithSignedRepresentative;
use lattirust_arithmetic::ring::PolyRing;
use lattirust_arithmetic::traits::{FromRandomBytes, WithL2Norm};
use lattirust_util::{check, check_eq};
use relations::principal_relation::{Index, Instance, PrincipalRelation, Witness};
use relations::Relation;

use crate::common_reference_string::CommonReferenceString;
use crate::shared::{compute_a__, compute_phi, compute_phi__, fold_instance, TranscriptView};

pub fn verify_principal_relation_oneround<'a, R: PolyRing>(
    arthur: &mut Arthur,
    crs: &'a CommonReferenceString<R>,
    index: &'a Index<R>,
    instance: &'a Instance<R>,
) -> ProofResult<(Index<R>, Instance<R>)>
where
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
    <R as PolyRing>::BaseRing: WithSignedRepresentative,
    <<R as PolyRing>::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
{
    let transcript = verify_core(crs, index, instance, arthur)?;
    let (index_next, instance_next) = fold_instance(&crs, &instance, &transcript);
    Ok((index_next, instance_next))
}

/// Verify consistency for one instance of the core Labrador protocol, used in each step of the recursion
pub fn verify_core<'a, R: PolyRing>(
    crs: &'a CommonReferenceString<R>,
    index: &'a Index<R>,
    instance: &'a Instance<R>,
    arthur: &mut Arthur,
) -> ProofResult<TranscriptView<R>>
where
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
{
    let (n, r) = (index.n, index.r);
    let num_constraints = instance.quad_dot_prod_funcs.len();
    let num_ct_constraints = instance.ct_quad_dot_prod_funcs.len();

    let u_1 = arthur
        .next_vector(crs.k1)
        .expect("error extracting prover message 1 from transcript");

    let num_projections = 256;
    let Pi = arthur
        .challenge_matrices::<R, WeightedTernaryChallengeSet<R>>(num_projections, n, r)
        .expect("error extracting verifier message 1 from transcript");

    let p = arthur
        .next_vector_canonical::<R::BaseRing>(num_projections)
        .expect("error extracting prover message 2 from transcript");
    let norm_p_sq = p.l2_norm_squared();
    let p_norm_bound_sq = 128f64 * index.norm_bound_squared;
    check!(
        norm_p_sq.to_f64().unwrap() <= p_norm_bound_sq,
        format!(
            "||p||_2^2 = {} must be <= 128*beta^2 = {}",
            norm_p_sq, p_norm_bound_sq
        )
    );

    let psi = arthur
        .challenge_vectors::<R::BaseRing, R::BaseRing>(num_ct_constraints, crs.num_aggregs)
        .expect("error extracting verifier message 2 (psi) from transcript");
    let omega = arthur
        .challenge_vectors::<R::BaseRing, R::BaseRing>(num_projections, crs.num_aggregs)
        .expect("error extracting verifier message 2 (omega) from transcript");

    let b__ = arthur
        .next_vec::<R>(crs.num_aggregs)
        .expect("error extracting prover message 3 from transcript");

    for k in 0..crs.num_aggregs {
        let mut rhs_k = omega[k].dot(&p);
        for l in 0..num_ct_constraints {
            rhs_k += psi[k][l] * instance.ct_quad_dot_prod_funcs[l].b;
        }
        check_eq!(b__[k].coefficients()[0], rhs_k);
    }

    let alpha = arthur
        .challenge_vector::<R, R>(num_constraints)
        .expect("error extracting verifier message 3 (alpha) from transcript");
    let beta = arthur
        .challenge_vector::<R, R>(crs.num_aggregs)
        .expect("error extracting verifier message 3 (beta) from transcript");

    let u_2 = arthur
        .next_vector(crs.k2)
        .expect("error extracting prover message 4 from transcript");

    let c = arthur
        .challenge_vec::<R, LabradorChallengeSet<R>>(crs.r)
        .expect("error extracting verifier message 4 from transcript");

    // Compute phi
    let phi__ = compute_phi__(crs, index, instance, &Pi, &psi, &omega);
    let phi = compute_phi(crs, instance, &alpha, &beta, &phi__);
    let a__ = compute_a__(crs, instance, &psi);

    Ok(TranscriptView {
        u_1,
        b__,
        alpha,
        beta,
        u_2,
        c,
        phi,
        a__,
    })
}

pub fn verify_principal_relation<R: PolyRing>(
    arthur: &mut Arthur,
    mut crs: &CommonReferenceString<R>,
    index: &Index<R>,
    instance: &Instance<R>,
) -> Result<(), ProofError>
where
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
    <R as PolyRing>::BaseRing: WithSignedRepresentative,
    <<R as PolyRing>::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
{
    let mut index_curr = index.clone();
    let mut instance_curr = instance.clone();

    while crs.next_crs.is_some() {
        (index_curr, instance_curr) =
            verify_principal_relation_oneround(arthur, crs, &index_curr, &instance_curr)?;
        crs = crs.next_crs.as_ref().unwrap();
    }
    let s = arthur.next_vectors(index_curr.n, index_curr.r)?;
    let witness = Witness::<R>::new(s);
    match PrincipalRelation::<R>::is_satisfied_err(&index_curr, &instance_curr, &witness) {
        Ok(_) => {
            debug!("└ Verifier::verify_principal_relation: OK");
            Ok(())
        }
        Err(e) => {
            debug!("└ Verifier::verify_principal_relation: ERROR {}", e);
            Err(ProofError::InvalidProof)
        }
    }
}
