use log::log_enabled;
use log::Level::Debug;
use nimue::{Arthur, ProofResult};
use tracing::debug;

use lattirust_arithmetic::balanced_decomposition::{decompose_matrix, recompose_matrix};
use lattirust_arithmetic::challenge_set::ternary::{mul_f_trit, TernaryChallengeSet, Trit};
use lattirust_arithmetic::linear_algebra::inner_products::inner_products_mat;
use lattirust_arithmetic::linear_algebra::SymmetricMatrix;
use lattirust_arithmetic::nimue::arthur::SerArthur;
use lattirust_arithmetic::nimue::traits::ChallengeFromRandomBytes;
use lattirust_arithmetic::ring::{ConvertibleRing, SignedRepresentative};

use crate::util::{norm_l2_columnwise, to_integers, PublicParameters, Witness};

pub struct Prover<F: ConvertibleRing> {
    _marker: std::marker::PhantomData<F>,
}

impl<F: ConvertibleRing> Prover<F> {
    #[tracing::instrument]
    pub fn merge(
        arthur: &mut Arthur,
        pp: &PublicParameters<F>,
        witness_1: Witness<F>,
        witness_2: Witness<F>,
    ) -> ProofResult<Witness<F>> {
        debug!("┌ Prover::merge");
        debug_assert_eq!(witness_1.ncols(), pp.inner_security_parameter);
        debug_assert_eq!(witness_2.ncols(), pp.inner_security_parameter);

        // Compute cross inner products witness_2^T * witness_1 (over the integers)
        let witness_2_int = witness_2.map(|x| Into::<SignedRepresentative>::into(x));
        let witness_1_int = witness_1.map(|x| Into::<SignedRepresentative>::into(x));
        let cross_terms = &witness_2_int.transpose() * &witness_1_int;
        arthur.absorb_matrix::<SignedRepresentative>(&cross_terms)?;
        arthur.ratchet()?;

        let mut witness = witness_1;
        witness.extend(witness_2.column_iter());

        Ok(witness)
    }

    #[tracing::instrument]
    pub fn reduce(
        arthur: &mut Arthur,
        pp: &PublicParameters<F>,
        witness: Witness<F>,
    ) -> ProofResult<Witness<F>> {
        debug!("┌ Prover::reduce");
        debug_assert_eq!(witness.ncols(), 2 * pp.inner_security_parameter);

        // Decompose witness matrix into wider matrix with small (column-wise) norms
        let decomp_witness =
            decompose_matrix(&witness, pp.decomposition_basis, pp.decomposition_length);
        debug_assert_eq!(decomp_witness.nrows(), witness.nrows());
        debug_assert_eq!(
            decomp_witness.ncols(),
            witness.ncols() * pp.decomposition_length
        );

        if log_enabled!(Debug) {
            let norm_decomp = PublicParameters::<F>::decomposed_norm_max(
                pp.decomposition_basis,
                pp.witness_len(),
            );
            debug!(
            "  (Assuming uniform witness norms), expected norm should be should be at most max = {norm_decomp}"
            );
            let decomp_norms = norm_l2_columnwise(&decomp_witness);
            debug!(
                "  Columns of decomposed witness have norms (min, mean, max) = ({}, {}, {})",
                decomp_norms.iter().min_by(|a, b| a.total_cmp(b)).unwrap(),
                decomp_norms.iter().sum::<f64>() / decomp_norms.len() as f64,
                decomp_norms.iter().max_by(|a, b| a.total_cmp(b)).unwrap()
            );
            debug_assert!(decomp_norms
                .into_iter()
                .all(|l2_norm| l2_norm <= norm_decomp));
        }

        // Commit to decomposed witness matrix (mod q)
        let committed_decomp_witness = &pp.commitment_mat * &decomp_witness;
        debug_assert_eq!(
            recompose_matrix(&committed_decomp_witness, pp.powers_of_basis().as_slice())
                .map(|x| Into::<SignedRepresentative>::into(x).0),
            (&pp.commitment_mat * &witness).map(|x| Into::<SignedRepresentative>::into(x).0),
            "commitments match"
        );

        // Add to FS transcript
        arthur.absorb_matrix::<F>(&committed_decomp_witness)?;
        arthur.ratchet()?;

        // Compute inner products over the integers
        let decomp_witness_int = &to_integers(&decomp_witness);
        let inner_products_decomp: SymmetricMatrix<SignedRepresentative> =
            inner_products_mat(decomp_witness_int);
        debug_assert_eq!(
            inner_products_decomp.size(),
            witness.ncols() * pp.decomposition_length
        );

        // Add to FS transcript
        arthur.absorb_symmetric_matrix::<SignedRepresentative>(&inner_products_decomp)?;
        arthur.ratchet()?;

        // Get challenge
        let challenge = arthur.challenge_matrix::<Trit, TernaryChallengeSet<Trit>>(
            2 * pp.inner_security_parameter * pp.decomposition_length,
            pp.inner_security_parameter,
        )?;
        arthur.ratchet()?;

        // Compute next witness
        // TODO: We don't actually have to do this mod q, but with the current implementation we're barely doing modular arithmetic, not sure if it makes sense to work over the integers instead here.
        let new_witness = mul_f_trit(&decomp_witness, &challenge);

        if log_enabled!(Debug) {
            // Check that the new witness does not grow in norm
            let new_norms = norm_l2_columnwise(&new_witness);
            debug!(
                "Columns of folded witness have norms (min, mean, max) = ({}, {}, {})",
                new_norms.iter().min_by(|a, b| a.total_cmp(b)).unwrap(),
                new_norms.iter().sum::<f64>() / new_norms.len() as f64,
                new_norms.iter().max_by(|a, b| a.total_cmp(b)).unwrap()
            );
            debug!(
                "(Assuming uniform witness norms), expected norm should be should give (mean, max) = ({}, {})",
                (2 * pp.decomposition_length * pp.inner_security_parameter * pp.decomposition_basis as usize) as f64 / 3.,
                (2 * pp.decomposition_length * pp.inner_security_parameter * pp.decomposition_basis as usize) as f64
            );
            debug_assert!(new_norms
                .into_iter()
                .all(|l2_norm| l2_norm <= pp.norm_bound));
        }
        Ok(new_witness)
    }

    #[tracing::instrument]
    pub fn fold(
        arthur: &mut Arthur,
        pp: &PublicParameters<F>,
        witness_1: Witness<F>,
        witness_2: Witness<F>,
    ) -> ProofResult<Witness<F>> {
        debug!("┌ Prover::fold");
        let merged_witness = Prover::merge(arthur, pp, witness_1, witness_2)?;
        let folded_witness = Prover::reduce(arthur, pp, merged_witness)?;
        debug!("└ Stop  Prover::fold");
        Ok(folded_witness)
    }
}
