use ark_relations::r1cs::ConstraintSystem;
use nimue::DuplexHash;
use nimue::hash::Unit;

use crate::labrador::binary_r1cs::util::{BinaryR1CSCRS, Z2};
use crate::labrador::common_reference_string::CommonReferenceString;
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::traits::FromRandomBytes;
use crate::nimue::iopattern::{LatticeIOPattern, SerIOPattern, SqueezeFromRandomBytes};
use crate::relations::labrador::principal_relation::PrincipalRelation;

pub trait LabradorIOPattern<R, H>
    where
        R: PolyRing,
        H: DuplexHash<u8>,
        LabradorChallengeSet<R>: FromRandomBytes<R>,
        WeightedTernaryChallengeSet<R>: FromRandomBytes<R>
{
    fn labrador_crs(self, crs: &CommonReferenceString<R>) -> Self;
    fn labrador_instance(self, instance: &PrincipalRelation<R>) -> Self;
    fn labrador_io(self, crs: &CommonReferenceString<R>) -> Self;
    fn labrador_binaryr1cs_io(self, r1cs: &ConstraintSystem<Z2>, crs: &BinaryR1CSCRS<R>) -> Self;
}

impl<R, H> LabradorIOPattern<R, H> for LatticeIOPattern<R, H> where
    R: PolyRing,
    H: DuplexHash<u8>,
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>
{
    fn labrador_crs(self, crs: &CommonReferenceString<R>) -> Self {
        self.absorb_serializable_like(crs, "labrador_principalrelation_crs")
    }

    fn labrador_instance(self, instance: &PrincipalRelation<R>) -> Self {
        self.absorb_serializable_like(instance, "labrador_principalrelation_crs")
    }

    fn labrador_io(self, crs: &CommonReferenceString<R>) -> Self {
        let log_q = R::modulus().next_power_of_two().ilog2() as f64;
        let num_aggregs = (128. / log_q).ceil() as usize;
        self.absorb_vector(crs.k1, "prover message 1")
            .squeeze_matrices::<R, WeightedTernaryChallengeSet<R>>(256, crs.n, crs.r, "verifier message 1")
            .absorb_vector_baseringelem(256, "prover message 2")
            .squeeze_vectors::<R::BaseRing, R::BaseRing>(crs.num_constant_constraints, num_aggregs, "verifier message 2 (psi)")
            .squeeze_vectors::<R::BaseRing, R::BaseRing>(256, num_aggregs, "verifier message 2 (omega)")
            .absorb_vec(num_aggregs, "prover message 3")
            .squeeze_vector::<R, R>(crs.num_constraints, "verifier message 3 (alpha)")
            .squeeze_vector::<R, R>(num_aggregs, "verifier message 3 (beta)")
            .absorb_vector(crs.k2, "prover message 4")
            .squeeze_vec::<R, LabradorChallengeSet<R>>(crs.r, "verifier message 4")
            .absorb_vector(crs.n, "prover message 5 (z)")
            .absorb_vectors(crs.k, crs.r, "prover message 5 (t)")
            .absorb_lower_triangular_matrix(crs.r, "prover message 5 (G)")
            .absorb_lower_triangular_matrix(crs.r, "prover message 5 (H)")
    }

    fn labrador_binaryr1cs_io(self, r1cs: &ConstraintSystem<Z2>, crs: &BinaryR1CSCRS<R>) -> Self {
        let k = r1cs.num_constraints;
        let secparam = 128;
        self.absorb_vector(k, "prover message 1 (t)")
            .squeeze_binary_matrix(secparam, k, "verifier message 1 (alpha)")
            .squeeze_binary_matrix(secparam, k, "verifier message 1 (beta)")
            .squeeze_binary_matrix(secparam, k, "verifier message 1 (gamma)")
            .absorb_vector_baseringelem(secparam, "prover message 2 (g)")
    }
}