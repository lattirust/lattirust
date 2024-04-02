use nimue::{Arthur, IOPatternError};

use crate::lattice_arithmetic::balanced_decomposition::decompose_matrix;
use crate::lattice_arithmetic::challenge_set::ternary::{mul_f_trit, TernaryChallengeSet, Trit};
use crate::lattice_arithmetic::matrix::Matrix;
use crate::lattice_arithmetic::poly_ring::ConvertibleField;
use crate::lattice_arithmetic::traits::WithL2Norm;
use crate::lova::util::{Instance, norm_l2_squared_columnwise, PublicParameters, SECPARAM, Witness};
use crate::nimue::arthur::SerArthur;
use crate::nimue::traits::ChallengeFromRandomBytes;

pub fn prove_folding<F: ConvertibleField>(arthur: &mut Arthur, pp: &PublicParameters<F>, instance: &Instance<F>, witness: &Witness<F>) -> Result<(Matrix<F>), IOPatternError> {
    debug_assert_eq!(witness.ncols(), 2 * SECPARAM);
    debug_assert!(norm_l2_squared_columnwise(witness).into_iter().all(|l2_norm_sq| l2_norm_sq as f64 <= pp.norm_bound));

    // in {x in F | ||x|| <= decomposition_basis}^(n x 2*SEC*decomposition_length)
    let decomp_witness = decompose_matrix(&witness, pp.decomposition_basis, pp.decomposition_length);

    // in F^(2*SEC*decomposition_length x m)
    let committed_decomp_witness = &pp.commitment_mat * &decomp_witness;

    // Add to FS transcript
    arthur.absorb_matrix(&committed_decomp_witness)?;

    // in {x in F | ||x|| <= norm_bound^2}^(2*SEC*decomposition_length x 2*SEC*decomposition_length)
    let inner_products = &decomp_witness.transpose() * &decomp_witness; // TODO: compute only lower triangular part

    // Add to FS transcript
    arthur.absorb_matrix(&inner_products)?;

    // Get challenge
    let challenge = arthur.challenge_matrix::<Trit, TernaryChallengeSet<Trit>>(2 * SECPARAM * pp.decomposition_length, SECPARAM)?;

    // Compute next witness
    let new_witness = mul_f_trit(&decomp_witness, &challenge);

    // Check that the new witness does not grow in norm
    debug_assert!(norm_l2_squared_columnwise(&new_witness).into_iter().all(|l2_norm_sq| l2_norm_sq as f64 <= pp.norm_bound));

    Ok(new_witness)
}