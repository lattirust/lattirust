#![allow(non_snake_case)]

use nimue::{Merlin, ProofResult};
use tracing::{event, instrument, Level};

use lattirust_arithmetic::balanced_decomposition::DecompositionFriendlySignedRepresentative;
use lattirust_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::nimue::merlin::SerMerlin;
use lattirust_arithmetic::nimue::traits::ChallengeFromRandomBytes;
use lattirust_arithmetic::ring::PolyRing;
use lattirust_arithmetic::ring::representatives::WithSignedRepresentative;
use lattirust_arithmetic::traits::FromRandomBytes;
use relations::{principal_relation, Relation};

use crate::binary_r1cs::BinaryR1CS;
use crate::binary_r1cs::util::{BinaryR1CSCRS, BinaryR1CSTranscript, reduce};
use crate::prover::prove_principal_relation;
use crate::util::{concat, embed, lift};

#[instrument(
    name = "BinR1CS -> PR",
    level = "info",
    skip(pp, merlin, index, instance, witness)
)]
pub fn prove_reduction_binaryr1cs_labradorpr<'a, R: PolyRing>(
    pp: &BinaryR1CSCRS<R>,
    merlin: &'a mut Merlin,
    index: &<BinaryR1CS as Relation>::Index,
    instance: &<BinaryR1CS as Relation>::Instance,
    witness: &<BinaryR1CS as Relation>::Witness,
) -> (
    principal_relation::Index<R>,
    principal_relation::Instance<R>,
    principal_relation::Witness<R>,
)
where
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
{
    let (A, B, C) = (&index.a, &index.b, &index.c);

    let w = concat(vec![instance.0.as_slice(), witness.0.as_slice()].as_slice());
    let (k, n) = (pp.num_constraints, pp.num_variables);
    assert_eq!(k, n, "the current implementation only support k = n"); // TODO: remove this restriction by splitting a,b,c or w into multiple vectors

    // TODO: add statement to merlin

    let a = A * &w;
    let b = B * &w;
    let c = C * &w;

    event!(Level::DEBUG, "computing A*w, B*w, C*w");

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
    let t = &pp.A * &v;

    merlin.absorb_vector(&t).unwrap();

    event!(
        Level::DEBUG,
        "squeezing alpha in {{0,1}}^{}x{k}",
        pp.security_parameter
    );
    let alpha = merlin
        .challenge_binary_matrix(pp.security_parameter, k)
        .unwrap();
    event!(
        Level::DEBUG,
        "squeezing beta in {{0,1}}^{}x{k}",
        pp.security_parameter
    );
    let beta = merlin
        .challenge_binary_matrix(pp.security_parameter, k)
        .unwrap();
    event!(
        Level::DEBUG,
        "squeezing gamma in {{0,1}}^{}x{k}",
        pp.security_parameter
    );
    let gamma = merlin
        .challenge_binary_matrix(pp.security_parameter, k)
        .unwrap();

    // delta_i is computed mod 2, i.e., over Z2
    let delta = &alpha * A + &beta * B + &gamma * C;

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
        pp.security_parameter,
        "g has length {} but should have length {}",
        g.len(),
        pp.security_parameter
    );
    event!(Level::DEBUG, "absorbing g in R_q^{}", pp.security_parameter);
    merlin.absorb_vector_canonical::<R::BaseRing>(&g).unwrap();

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

    let index_pr = pp.core_crs.clone(); // TODO: this should be derived separately

    let instance_pr = reduce(pp, &transcript);
    let i =  reduce(pp, &transcript);
    assert_eq!(i, instance_pr, "reduction should be deterministic");

    let witness_pr = principal_relation::Witness::<R> {
        s: vec![a_R, b_R, c_R, w_R, a_tilde, b_tilde, c_tilde, w_tilde], // see definition of indices above
    };

    (index_pr, instance_pr, witness_pr)
}

pub fn prove_binary_r1cs<'a, R: PolyRing>(
    pp: &BinaryR1CSCRS<R>,
    merlin: &'a mut Merlin,
    index: &<BinaryR1CS as Relation>::Index,
    instance: &<BinaryR1CS as Relation>::Instance,
    witness: &<BinaryR1CS as Relation>::Witness,
) -> ProofResult<&'a [u8]>
where
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
    <R as PolyRing>::BaseRing: WithSignedRepresentative,
    <<R as PolyRing>::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
{
    let (index_pr, instance_pr, witness_pr) =
        prove_reduction_binaryr1cs_labradorpr(pp, merlin, index, instance, witness);

    merlin.ratchet()?;

    prove_principal_relation(merlin, &index_pr, &instance_pr, &witness_pr, &pp.core_crs)
}
