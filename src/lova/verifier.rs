use nimue::{Merlin, ProofError};

use crate::{check, check_eq};
use crate::lattice_arithmetic::balanced_decomposition::decompose_matrix;
use crate::lattice_arithmetic::challenge_set::ternary::{TernaryChallengeSet, Trit};
use crate::lattice_arithmetic::matrix::Matrix;
use crate::lattice_arithmetic::poly_ring::ConvertibleField;
use crate::lattice_arithmetic::traits::WithL2Norm;
use crate::lova::util::{CRS, SECPARAM};
use crate::nimue::merlin::SerMerlin;
use crate::nimue::traits::ChallengeFromRandomBytes;

pub fn verify_folding<F: ConvertibleField>(merlin: &mut Merlin, crs: &CRS<F>, t: &Matrix<F>) -> Result<Matrix<F>, ProofError> {
    // Check W = decompose(T)
    // TODO: or recompute ourselves? We'll need to do a pass on the full matrix anyway
    let w = merlin.next_matrix(crs.commitment_mat.nrows(), 2 * SECPARAM * crs.decomposition_length)?;
    check_eq!(w, decompose_matrix(&t, crs.decomposition_basis, crs.decomposition_length));

    let g = merlin.next_symmetric_matrix::<F>(2 * SECPARAM)?;

    let c = merlin.challenge_matrix::<Trit, TernaryChallengeSet<Trit>>(2 * SECPARAM * crs.decomposition_length, SECPARAM)?;

    let t_new = mul_F_Trit(&w, &c);
    Ok(t_new)
}

pub fn decide<F: ConvertibleField>(merlin: &mut Merlin, crs: &CRS<F>, t: &Matrix<F>, z: &Matrix<F>) -> Result<(), ProofError> {
    // Check W = decompose(T)
    // TODO: or recompute ourselves? We'll need to do a pass on the full matrix anyway
    let w = merlin.next_matrix(crs.commitment_mat.nrows(), 2 * SECPARAM * crs.decomposition_length)?;
    check_eq!(w, decompose_matrix(&t, crs.decomposition_basis, crs.decomposition_length));

    let g = merlin.next_symmetric_matrix(2 * SECPARAM)?;

    let c = merlin.challenge_matrix::<Trit, TernaryChallengeSet<Trit>>(2 * SECPARAM * crs.decomposition_length, SECPARAM)?;

    // Check ||Z|| <= beta_z
    let norm_bound_z = crs.norm_bound * (2 * SECPARAM) as f64; // TODO
    let norm_bound_z_sq = norm_bound_z.powi(2);
    for z_i in z.row_iter() {
        check!(z_i.transpose().l2_norm_squared() as f64 <= norm_bound_z_sq);
    }

    // Check AZ == T
    check_eq!(&crs.commitment_mat * z, w);

    // Check <Z_i, Z_i> == C_i * G * C^T_i
    let rhs = mul_Trit_F(c.transpose(), &mul_F_Trit_sym(&g, &c));
    debug_assert_eq!(rhs.nrows(), 1);
    debug_assert_eq!(rhs.ncols(), SECPARAM);
    for (z_i, rhs_i) in z.row_iter().zip(rhs.as_slice()) {
        let lhs = z_i.dot(&z_i);
        check_eq!(lhs, *rhs_i);
    }

    Ok(())
}