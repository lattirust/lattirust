use std::ops::Mul;
use log::debug;
use nimue::{Arthur, ProofError, ProofResult};

use lattirust_arithmetic::balanced_decomposition::{
    recompose_left_right_symmetric_matrix, recompose_matrix,
};
use lattirust_arithmetic::challenge_set::ternary::{
    mul_f_trit, mul_trit_transpose_sym_trit, TernaryChallengeSet, Trit,
};
use lattirust_arithmetic::linear_algebra::SymmetricMatrix;
use lattirust_arithmetic::nimue::arthur::SerArthur;
use lattirust_arithmetic::nimue::traits::ChallengeFromRandomBytes;
use lattirust_arithmetic::ring::representatives::WithSignedRepresentative;
use lattirust_arithmetic::ring::Ring;
use lattirust_util::{check, check_eq};

use crate::util::{Instance, PublicParameters};

pub struct Verifier<F: Ring> {
    _marker: std::marker::PhantomData<F>,
}

impl<F: Ring + WithSignedRepresentative> Verifier<F>
where for<'a> &'a F: Mul<&'a F, Output = F>,
{
    #[tracing::instrument]
    pub fn merge(
        arthur: &mut Arthur,
        pp: &PublicParameters<F>,
        instance_1: Instance<F>,
        instance_2: Instance<F>,
    ) -> ProofResult<Instance<F>> {
        debug!("┌ Verifier::merge");
        debug_assert_eq!(instance_1.commitment.ncols(), pp.inner_security_parameter);
        debug_assert_eq!(instance_2.commitment.ncols(), pp.inner_security_parameter);

        let cross_terms = arthur
            .next_matrix::<F>(
                pp.inner_security_parameter,
                pp.inner_security_parameter,
            )
            .unwrap();
        arthur.ratchet().unwrap();

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
        arthur: &mut Arthur,
        pp: &PublicParameters<F>,
        instance: Instance<F>,
    ) -> Result<Instance<F>, ProofError> {
        debug!("┌ Verifier::reduce");
        debug_assert_eq!(instance.commitment.ncols(), 2 * pp.inner_security_parameter);
        let committed_decomp_witness = arthur
            .next_matrix(
                pp.commitment_mat.nrows(),
                2 * pp.inner_security_parameter * pp.decomposition_length,
            )
            .unwrap();
        arthur.ratchet().unwrap();

        let inner_products_decomp = arthur
            .next_symmetric_matrix::<F>(
                2 * pp.inner_security_parameter * pp.decomposition_length,
            )
            .unwrap();
        arthur.ratchet().unwrap();

        let challenge = arthur
            .challenge_matrix::<Trit, TernaryChallengeSet<Trit>>(
                2 * pp.inner_security_parameter * pp.decomposition_length,
                pp.inner_security_parameter,
            )
            .unwrap();
        arthur.ratchet().unwrap();

        // Check G^T * inner_products_decomp * G == instance.inner_products (over the integers)
        let inner_products_recomp: SymmetricMatrix<F> =
            recompose_left_right_symmetric_matrix(
                &inner_products_decomp,
                pp.powers_of_basis().as_slice(),
            );

        check_eq!(
            inner_products_recomp,
            instance.inner_products,
            "inner products match"
        );

        check_eq!(
            recompose_matrix(&committed_decomp_witness, pp.powers_of_basis().as_slice()),
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
        arthur: &mut Arthur,
        pp: &PublicParameters<F>,
        instance_1: Instance<F>,
        instance_2: Instance<F>,
    ) -> ProofResult<Instance<F>> {
        let merged_instance = Self::merge(arthur, pp, instance_1, instance_2).unwrap();
        let folded_instance = Self::reduce(arthur, pp, merged_instance).unwrap();
        Ok(folded_instance)
    }
}
