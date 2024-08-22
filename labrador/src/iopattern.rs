use ark_relations::r1cs::{ConstraintSystem, ConstraintSystemRef};
use nimue::{DuplexHash, IOPattern};

use lattirust_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::nimue::iopattern::{
    RatchetIOPattern, SerIOPattern, SqueezeFromRandomBytes,
};
use lattirust_arithmetic::ring::PolyRing;
use lattirust_arithmetic::traits::FromRandomBytes;

use crate::binary_r1cs::util::{BinaryR1CSCRS, Z2};
use crate::common_reference_string::CommonReferenceString;

pub trait LabradorIOPattern<R, H>:
    SerIOPattern + SqueezeFromRandomBytes + RatchetIOPattern
where
    R: PolyRing,
    H: DuplexHash<u8>,
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
{
    // fn labrador_crs(self, crs: &CommonReferenceString<R>) -> Self {
    //     self.absorb_serializable_like(crs, "labrador_principalrelation_crs")
    // }

    // fn labrador_instance(self, instance: &PrincipalRelation<R>) -> Self {
    //     self.absorb_serializable_like(instance, "labrador_principalrelation_crs")
    // }

    fn labrador_io(self, crs: &CommonReferenceString<R>) -> Self {
        let log_q = R::modulus().bits() as f64;
        let num_aggregs = (128. / log_q).ceil() as usize;
        self.absorb_vector::<R>(crs.k1, "prover message 1")
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
    fn labrador_binaryr1cs_io(self, r1cs: &ConstraintSystemRef<Z2>, crs: &BinaryR1CSCRS<R>) -> Self {
        let k = r1cs.num_constraints();
        let secparam = crs.security_parameter;
        self.absorb_vector::<R>(crs.A.nrows(), "prover message 1 (t)")
            .squeeze_binary_matrix(secparam, k, "verifier message 1 (alpha)")
            .squeeze_binary_matrix(secparam, k, "verifier message 1 (beta)")
            .squeeze_binary_matrix(secparam, k, "verifier message 1 (gamma)")
            .absorb_vector_canonical::<R::BaseRing>(secparam, "prover message 2 (g)")
    }
}

impl<R, H> LabradorIOPattern<R, H> for IOPattern<H>
where
    R: PolyRing,
    H: DuplexHash<u8>,
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
    Self: SerIOPattern + SqueezeFromRandomBytes + RatchetIOPattern,
{
}
