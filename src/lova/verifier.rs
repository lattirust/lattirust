use nimue::{Merlin, ProofError};

use crate::{check, check_eq};
use crate::lattice_arithmetic::balanced_decomposition::recompose_matrix;
use crate::lattice_arithmetic::challenge_set::ternary::{mul_f_trit, mul_trit_transpose_sym_trit, TernaryChallengeSet, Trit};
use crate::lattice_arithmetic::poly_ring::ConvertibleField;
use crate::lova::util::{Instance, PublicParameters};
use crate::nimue::merlin::SerMerlin;
use crate::nimue::traits::ChallengeFromRandomBytes;

pub fn verify_folding<F: ConvertibleField>(merlin: &mut Merlin, pp: &PublicParameters<F>, instance_1: &Instance<F>, instance_2: &Instance<F>) -> Result<Instance<F>, ProofError> {
    let committed_decomp_witness = merlin.next_matrix(pp.commitment_mat.nrows(), 2 * pp.security_parameter * pp.decomposition_length)?;
    merlin.ratchet()?;

    let inner_products_decomp = merlin.next_symmetric_matrix::<i128>(2 * pp.security_parameter * pp.decomposition_length)?.into();
    merlin.ratchet()?;

    let challenge = merlin.challenge_matrix::<Trit, TernaryChallengeSet<Trit>>(2 * pp.security_parameter * pp.decomposition_length, pp.security_parameter)?;
    merlin.ratchet()?;

    // Check G^T * inner_products_decomp_1 * G == instance_1.inner_products
    // TODO

    // Check G^T * inner_products_decomp_2 * G == instance_2.inner_products
    // TODO

    // Check committed_decomp_witness * G == (instance_1.commitment || instance_2.commitment) (mod q)
    let mut commitment = instance_1.commitment.clone();
    commitment.extend(instance_2.commitment.column_iter());
    check_eq!(recompose_matrix(&committed_decomp_witness, &pp.powers_of_basis().as_slice()), commitment);

    // Compute new instance (commitment and inner products)
    let commitment_new = mul_f_trit(&committed_decomp_witness, &challenge);
    let inner_products_new = mul_trit_transpose_sym_trit(&inner_products_decomp, &challenge);
    Ok(Instance { commitment: commitment_new, inner_products: inner_products_new })
}