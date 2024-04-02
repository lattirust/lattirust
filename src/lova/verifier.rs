use nimue::{Merlin, ProofError};

use crate::{check, check_eq};
use crate::labrador::util::mul_basescalar_vector;
use crate::lattice_arithmetic::balanced_decomposition::{decompose_matrix, recompose_matrix};
use crate::lattice_arithmetic::challenge_set::ternary::{mul_f_trit, mul_f_trit_sym, mul_trit_f, mul_trit_transpose_sym_trit, TernaryChallengeSet, Trit};
use crate::lattice_arithmetic::matrix::Matrix;
use crate::lattice_arithmetic::poly_ring::ConvertibleField;
use crate::lattice_arithmetic::traits::WithL2Norm;
use crate::lova::util::{Instance, PublicParameters, SECPARAM, Witness};
use crate::nimue::merlin::SerMerlin;
use crate::nimue::traits::ChallengeFromRandomBytes;

pub fn verify_folding<F: ConvertibleField>(merlin: &mut Merlin, crs: &PublicParameters<F>, instance: &Instance<F>) -> Result<Instance<F>, ProofError> {
    let committed_decomp_witness = merlin.next_matrix(crs.commitment_mat.nrows(), 2 * SECPARAM * crs.decomposition_length)?;

    let inner_products = merlin.next_symmetric_matrix::<F>(2 * SECPARAM)?.into();

    let challenge = merlin.challenge_matrix::<Trit, TernaryChallengeSet<Trit>>(2 * SECPARAM * crs.decomposition_length, SECPARAM)?;

    // Check WG = T (mod q)
    check_eq!(recompose_matrix(&committed_decomp_witness, &crs.powers_of_basis().as_slice()), instance.commitment);

    let commitment_new = mul_f_trit(&committed_decomp_witness, &challenge);
    let inner_products_new = mul_trit_transpose_sym_trit(&inner_products, &challenge);
    Ok(Instance{commitment: commitment_new, inner_products: inner_products_new})
}