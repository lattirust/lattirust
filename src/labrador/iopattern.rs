use nimue::DuplexHash;
use nimue::hash::Unit;

use crate::labrador::binary_r1cs::prover_binary_r1cs::{BinaryR1CSCRS, BinaryR1CSInstance};
use crate::labrador::common_reference_string::CommonReferenceString;
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::traits::{FromRandomBytes, WithLog2};
use crate::nimue::iopattern::LatticeIOPattern;
use crate::relations::labrador::principal_relation::PrincipalRelation;

pub trait LabradorIOPattern<R, H, U = u8>
    where
        R: PolyRing,
        H: DuplexHash<U>,
        U: Unit,
        LabradorChallengeSet<R>: FromRandomBytes<R>,
        WeightedTernaryChallengeSet<R>: FromRandomBytes<R>
{
    fn labrador_statement(self) -> Self;
    fn labrador_io(self, instance: &PrincipalRelation<R>, crs: &CommonReferenceString<R>) -> Self;
    fn labrador_binaryr1cs_io(self, r1cs: &BinaryR1CSInstance, crs: &BinaryR1CSCRS<R>) -> Self;
}

impl<R, H, U> LabradorIOPattern<R, H, U> for LatticeIOPattern<R, H, U> where
    R: PolyRing,
    H: DuplexHash<U>,
    U: Unit,
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>
{
    fn labrador_statement(self) -> Self {
        self // TODO: what is the statement for a Labrador principal relation?
    }

    fn labrador_io(self, instance: &PrincipalRelation<R>, crs: &CommonReferenceString<R>) -> Self {
        let num_aggregs = (128. / R::BaseRing::log2_q()).ceil() as usize;
        self.absorb_vector(crs.k1, "prover message 1")
            .squeeze_matrices::<R, WeightedTernaryChallengeSet<R>>(256, crs.n, crs.r, "verifier message 1")
            .absorb_vector_baseringelem(256, "prover message 2")
            .squeeze_vectors::<R::BaseRing, R::BaseRing>(instance.ct_quad_dot_prod_funcs.len(), num_aggregs, "verifier message 2 (psi)")
            .squeeze_vectors::<R::BaseRing, R::BaseRing>(256, num_aggregs, "verifier message 2 (omega)")
            .absorb_vec(num_aggregs, "prover message 3")
            .squeeze_vector::<R, R>(instance.quad_dot_prod_funcs.len(), "verifier message 3 (alpha)")
            .squeeze_vector::<R, R>(num_aggregs, "verifier message 3 (beta)")
            .absorb_vector(crs.k2, "prover message 4")
            .squeeze_vec::<R, LabradorChallengeSet<R>>(crs.r, "verifier message 4")
            .absorb_vector(crs.n, "prover message 5 (z)")
            .absorb_vectors(crs.k1, crs.r, "prover message 5 (t)")
            .absorb_lower_triangular_matrix(crs.r, "prover message 5 (G)")
            .absorb_lower_triangular_matrix(crs.r, "prover message 5 (H)")
    }

    fn labrador_binaryr1cs_io(self, r1cs: &BinaryR1CSInstance, crs: &BinaryR1CSCRS<R>) -> Self {
        let k = r1cs.A.nrows();
        let secparam = 128;
        self.absorb_vector(crs.A.nrows(), "prover message 1 (t)")
            .squeeze_binary_matrix(secparam, k, "verifier message 1 (alpha)")
            .squeeze_binary_matrix(secparam, k, "verifier message 1 (beta)")
            .squeeze_binary_matrix(secparam, k, "verifier message 1 (gamma)")
            .absorb_vector_baseringelem(secparam, "prover message 2 (g)")
    }
}