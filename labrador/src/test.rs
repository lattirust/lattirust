// use nimue::{Arthur, IOPattern, Merlin, ProofResult};
// use tracing_subscriber::fmt::format;
// use tracing_subscriber::fmt::format::FmtSpan;
// use lattirust_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
// use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
// 
// use lattirust_arithmetic::ring::{PolyRing, Pow2CyclotomicPolyRingNTT};
// use lattirust_arithmetic::ring::representatives::WithSignedRepresentative;
// use lattirust_arithmetic::ring::Zq2;
// use lattirust_arithmetic::traits::FromRandomBytes;
// use relations::r1cs::Size;
// use relations::reduction::Reduction;
// use relations::{test_completeness, test_soundness};
// use relations::principal_relation::PrincipalRelation;
// use crate::binary_r1cs::{BinaryR1CS, ReductionBinaryR1CSPrincipalRelation};
// use crate::binary_r1cs::prover::prove_reduction_binaryr1cs_labradorpr;
// use crate::binary_r1cs::util::BinaryR1CSCRS;
// use crate::binary_r1cs::verifier::verify_reduction_binaryr1cs_labradorpr;
// 
// pub type Z64 = Zq2<274177, 67280421310721>; // Q = 2^64+1
// const D: usize = 64;
// 
// type R = Pow2CyclotomicPolyRingNTT<Z64, D>;
// 
// fn init() {
//     let _ = tracing_subscriber::fmt::fmt()
//         .with_span_events(FmtSpan::ENTER | FmtSpan::CLOSE)
//         .event_format(format().compact())
//         .with_env_filter("none,labrador=trace")
//         .try_init();
// }
// 
// type TestReduction = Labrador<R>;
// 
// pub struct Labrador<R: PolyRing> {
//     _marker: std::marker::PhantomData<R>,
// }
// 
// impl<R: PolyRing> Reduction<BinaryR1CS, PrincipalRelation<R>, BinaryR1CSCRS<R>>
// for Labrador<R>
// where
//     LabradorChallengeSet<R>: FromRandomBytes<R>,
//     WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
//     <R as PolyRing>::BaseRing: WithSignedRepresentative,
// {
//     fn iopattern(
//         pp: &BinaryR1CSCRS<R>,
//         index_in_: &Self::IndexIn,
//         instance_in_: &Self::InstanceIn,
//     ) -> IOPattern {
//         let k = pp.num_constraints;
//         let secparam = pp.security_parameter;
//         IOPattern::new("reduction_binaryr1cs_principalrelation")
//             .absorb_vector::<R>(pp.A.nrows(), "prover message 1 (t)")
//             .squeeze_binary_matrix(secparam, k, "verifier message 1 (alpha)")
//             .squeeze_binary_matrix(secparam, k, "verifier message 1 (beta)")
//             .squeeze_binary_matrix(secparam, k, "verifier message 1 (gamma)")
//             .absorb_vector_canonical::<R::BaseRing>(secparam, "prover message 2 (g)")
//     }
// 
//     fn prove(
//         pp: &BinaryR1CSCRS<R>,
//         index: &Self::IndexIn,
//         instance: &Self::InstanceIn,
//         witness: &Self::WitnessIn,
//         merlin: &mut Merlin,
//     ) -> ProofResult<(Self::IndexOut, Self::InstanceOut, Self::WitnessOut)> {
//         Ok(prove_reduction_binaryr1cs_labradorpr(
//             pp, merlin, index, instance, witness,
//         ))
//     }
// 
//     fn verify(
//         pp: &BinaryR1CSCRS<R>,
//         index_in: &Self::IndexIn,
//         instance_in: &Self::InstanceIn,
//         arthur: &mut Arthur,
//     ) -> ProofResult<(Self::IndexOut, Self::InstanceOut)>
//     {
//         verify_reduction_binaryr1cs_labradorpr(arthur, pp, index_in, instance_in)
//     }
// }
// 
// const TEST_SIZE: Size = Size {
//     num_constraints: D * 4,
//     num_instance_variables: D,
//     num_witness_variables: D * 3,
// };
// 
// test_completeness!(
//     TestReduction,
//     BinaryR1CSCRS::new(
//         TEST_SIZE.num_constraints,
//         TEST_SIZE.num_instance_variables + TEST_SIZE.num_witness_variables
//     ),
//     TEST_SIZE
// );
// 
// test_soundness!(
//     TestReduction,
//     BinaryR1CSCRS::new(
//         TEST_SIZE.num_constraints,
//         TEST_SIZE.num_instance_variables + TEST_SIZE.num_witness_variables
//     ),
//     TEST_SIZE
// );
