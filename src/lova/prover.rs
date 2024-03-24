use ark_std::iterable::Iterable;
use nimue::{Arthur, IOPatternError};

use crate::lattice_arithmetic::balanced_decomposition::decompose_matrix;
use crate::lattice_arithmetic::challenge_set::ternary::{mul_f_trit, TernaryChallengeSet, Trit};
use crate::lattice_arithmetic::matrix::Matrix;
use crate::lattice_arithmetic::poly_ring::ConvertibleField;
use crate::lattice_arithmetic::traits::WithL2Norm;
use crate::lova::util::{CRS, SECPARAM};
use crate::nimue::arthur::SerArthur;
use crate::nimue::traits::ChallengeFromRandomBytes;

pub fn prove_folding<F: ConvertibleField>(arthur: &mut Arthur, crs: &CRS<F>, t: &Matrix<F>, s: &Matrix<F>) -> Result<(Matrix<F>), IOPatternError> {
    let n = s.nrows();
    debug_assert_eq!(s.ncols(), 2 * SECPARAM);

    // d in {x in F | ||x|| <= decomposition_basis}^(n x 2*SEC*decomposition_length)
    let d = decompose_matrix(&s, crs.decomposition_basis, crs.decomposition_length);

    debug_assert_eq!(d.nrows(), n);
    debug_assert_eq!(d.ncols(), 2 * SECPARAM * crs.decomposition_length);

    // w in F^(2*SEC*decomposition_length x m)
    let w = &crs.commitment_mat * &d;

    // Add w to FS transcript
   arthur.absorb_matrix(&w)?;

    // g in {x in F | ||x|| <= norm_bound^2}^(2*SEC*decomposition_length x 2*SEC*decomposition_length)
    let g = &d.transpose() * &d; // TODO: compute only lower triangular part

    // Add g to FS transcript
    for g_i in g.row_iter() {
        for g_ij in g_i.iter() {
            arthur.absorb_canonical_serializable(g_ij)?;
        }
    }

    // Get challenge
    let c = arthur.challenge_matrix::<Trit, TernaryChallengeSet<Trit>>(2 * SECPARAM * crs.decomposition_length, SECPARAM)?;

    // Compute z
    let z = mul_f_trit(&d, &c);

    for z_i in z.iter() {
        debug_assert!(z_i.l2_norm_squared() as f64  <= crs.norm_bound);
    }

    Ok(z)
}