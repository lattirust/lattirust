#![allow(dead_code)]
use nimue::{Arthur, IOPattern, Merlin, ProofResult};
use tracing_subscriber::fmt::format;
use tracing_subscriber::fmt::format::FmtSpan;

use lattirust_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::decomposition::DecompositionFriendlySignedRepresentative;
use lattirust_arithmetic::nimue::iopattern::{SerIOPattern, SqueezeFromRandomBytes};
use lattirust_arithmetic::ring::{PolyRing, Pow2CyclotomicPolyRingNTT};
use lattirust_arithmetic::ring::representatives::WithSignedRepresentative;
use lattirust_arithmetic::ring::Zq2;
use lattirust_arithmetic::traits::FromRandomBytes;
use relations::{test_completeness_with_init, test_soundness_with_init};
use relations::principal_relation::PrincipalRelation;
use relations::principal_relation::Size;
use relations::reduction::Reduction;

use crate::common_reference_string::CommonReferenceString;
use crate::prover::prove_principal_relation_oneround;
use crate::verifier::verify_principal_relation_oneround;

// Q = 2^64+1
const Q1: u64 = 274177;
const Q2: u64 = 67280421310721;
pub type Z64 = Zq2<Q1, Q2>;
const D: usize = 64;

type R = Pow2CyclotomicPolyRingNTT<Z64, D>;

fn init() {
    let _ = tracing_subscriber::fmt::fmt()
        .with_span_events(FmtSpan::ENTER | FmtSpan::CLOSE)
        .event_format(format().compact())
        .with_env_filter("none,labrador=trace")
        .try_init();
}

pub struct Labrador<R: PolyRing> {
    _marker: std::marker::PhantomData<R>,
}

impl<R: PolyRing> Reduction<PrincipalRelation<R>, PrincipalRelation<R>, CommonReferenceString<R>>
    for Labrador<R>
where
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
    <R as PolyRing>::BaseRing: WithSignedRepresentative,
    <<R as PolyRing>::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
{
    fn iopattern(
        crs: &CommonReferenceString<R>,
        _index_in: &Self::IndexIn,
        _instance_in: &Self::InstanceIn,
    ) -> IOPattern {
        let log_q = R::modulus().bits() as f64;
        let num_aggregs = (128. / log_q).ceil() as usize;
        IOPattern::new("reduction_binaryr1cs_principalrelation")
            .absorb_vector::<R>(crs.k1, "prover message 1")
            .squeeze_matrices::<R, WeightedTernaryChallengeSet<R>>(
                256,
                crs.n,
                crs.r,
                "verifier message 1",
            )
            .absorb_vector_canonical::<R::BaseRing>(256, "prover message 2")
            .squeeze_vectors::<R::BaseRing, R::BaseRing>(
                crs.num_constant_constraints,
                num_aggregs,
                "verifier message 2 (psi)",
            )
            .squeeze_vectors::<R::BaseRing, R::BaseRing>(
                256,
                num_aggregs,
                "verifier message 2 (omega)",
            )
            .absorb_vec::<R>(num_aggregs, "prover message 3")
            .squeeze_vector::<R, R>(crs.num_constraints, "verifier message 3 (alpha)")
            .squeeze_vector::<R, R>(num_aggregs, "verifier message 3 (beta)")
            .absorb_vector::<R>(crs.k2, "prover message 4")
            .squeeze_vec::<R, LabradorChallengeSet<R>>(crs.r, "verifier message 4")
            .absorb_vector::<R>(crs.n, "prover message 5 (z)")
            .absorb_vectors::<R>(crs.k, crs.r, "prover message 5 (t)")
            .absorb_symmetric_matrix::<R>(crs.r, "prover message 5 (G)")
            .absorb_symmetric_matrix::<R>(crs.r, "prover message 5 (H)")
    }

    fn prove(
        pp: &CommonReferenceString<R>,
        index: &Self::IndexIn,
        instance: &Self::InstanceIn,
        witness: &Self::WitnessIn,
        merlin: &mut Merlin,
    ) -> ProofResult<(Self::IndexOut, Self::InstanceOut, Self::WitnessOut)> {
        prove_principal_relation_oneround(merlin, pp, index, instance, witness)
    }

    fn verify(
        pp: &CommonReferenceString<R>,
        index_in: &Self::IndexIn,
        instance_in: &Self::InstanceIn,
        arthur: &mut Arthur,
    ) -> ProofResult<(Self::IndexOut, Self::InstanceOut)> {
        verify_principal_relation_oneround(arthur, pp, index_in, instance_in)
    }
}

const TEST_SIZE: Size = Size {
    num_witnesses: 2,
    witness_len: 64,
    norm_bound_sq: 4294967296., // sqrt(Q1*Q2)
    num_constraints: 1,
    num_constant_constraints: 1,
};

test_completeness_with_init!(
    Labrador<R>,
    CommonReferenceString::new_for_size(TEST_SIZE),
    TEST_SIZE,
    init
);

test_soundness_with_init!(
    Labrador<R>,
    CommonReferenceString::new_for_size(TEST_SIZE),
    TEST_SIZE,
    init
);
