use log::debug;
use nimue::{Merlin, ProofError, ProofResult};

use lattirust_arithmetic::balanced_decomposition::{
    recompose_left_right_symmetric_matrix, recompose_matrix,
};
use lattirust_arithmetic::challenge_set::ternary::{
    mul_f_trit, mul_trit_transpose_sym_trit, TernaryChallengeSet, Trit,
};
use lattirust_arithmetic::linear_algebra::SymmetricMatrix;
use lattirust_arithmetic::nimue::merlin::SerMerlin;
use lattirust_arithmetic::nimue::traits::ChallengeFromRandomBytes;
use lattirust_arithmetic::ring::{ConvertibleRing, SignedRepresentative};
use lattirust_util::{check, check_eq};

use crate::util::{Instance, PublicParameters};

pub struct Verifier<F: ConvertibleRing> {
    _marker: std::marker::PhantomData<F>,
}

impl<F: ConvertibleRing> Verifier<F> {
    #[tracing::instrument]
    pub fn merge(
        merlin: &mut Merlin,
        pp: &PublicParameters<F>,
        instance_1: Instance<F>,
        instance_2: Instance<F>,
    ) -> ProofResult<Instance<F>> {
        debug_assert_eq!(instance_1.commitment.ncols(), pp.security_parameter);
        debug_assert_eq!(instance_2.commitment.ncols(), pp.security_parameter);

        let cross_terms = merlin
            .next_matrix::<SignedRepresentative>(pp.security_parameter, pp.security_parameter)
            .unwrap();
        merlin.ratchet().unwrap();

        let mut commitment = instance_1.commitment;
        commitment.extend(instance_2.commitment.column_iter());

        let inner_products = SymmetricMatrix::from_blocks(
            instance_1.inner_products,
            cross_terms,
            instance_2.inner_products,
        );

        Ok(Instance {
            commitment,
            inner_products,
        })
    }

    #[tracing::instrument]
    pub fn reduce(
        merlin: &mut Merlin,
        pp: &PublicParameters<F>,
        instance: Instance<F>,
    ) -> Result<Instance<F>, ProofError> {
        debug_assert_eq!(instance.commitment.ncols(), 2 * pp.security_parameter);
        let committed_decomp_witness = merlin
            .next_matrix(
                pp.commitment_mat.nrows(),
                2 * pp.security_parameter * pp.decomposition_length,
            )
            .unwrap();
        merlin.ratchet().unwrap();

        let inner_products_decomp = merlin
            .next_symmetric_matrix::<SignedRepresentative>(
                2 * pp.security_parameter * pp.decomposition_length,
            )
            .unwrap();
        merlin.ratchet().unwrap();

        let challenge = merlin
            .challenge_matrix::<Trit, TernaryChallengeSet<Trit>>(
                2 * pp.security_parameter * pp.decomposition_length,
                pp.security_parameter,
            )
            .unwrap();
        merlin.ratchet().unwrap();

        // Check G^T * inner_products_decomp * G == instance.inner_products (over the integers)
        let inner_products_recomp: SymmetricMatrix<SignedRepresentative> =
            recompose_left_right_symmetric_matrix(
                &inner_products_decomp,
                pp.powers_of_basis_int().as_slice(),
            );

        check_eq!(
            inner_products_recomp,
            instance.inner_products,
            "inner products match"
        );

        check_eq!(
            recompose_matrix(&committed_decomp_witness, &pp.powers_of_basis().as_slice()),
            instance.commitment,
            "commitments match"
        );

        // Compute new instance (commitment and inner products)
        let commitment_new = mul_f_trit(&committed_decomp_witness, &challenge);
        let inner_products_new = mul_trit_transpose_sym_trit(&inner_products_decomp, &challenge);
        Ok(Instance {
            commitment: commitment_new,
            inner_products: inner_products_new,
        })
    }

    #[tracing::instrument]
    pub fn fold(
        merlin: &mut Merlin,
        pp: &PublicParameters<F>,
        instance_1: Instance<F>,
        instance_2: Instance<F>,
    ) -> ProofResult<Instance<F>> {
        let merged_instance = Self::merge(merlin, pp, instance_1, instance_2).unwrap();
        let folded_instance = Self::reduce(merlin, pp, merged_instance).unwrap();
        Ok(folded_instance)
    }
}
