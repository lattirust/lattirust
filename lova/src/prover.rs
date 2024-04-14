use log::debug;
use nimue::{Arthur, IOPatternError};

use lattirust_arithmetic::balanced_decomposition::decompose_matrix;
use lattirust_arithmetic::challenge_set::ternary::{mul_f_trit, TernaryChallengeSet, Trit};
use lattirust_arithmetic::linear_algebra::inner_products::inner_products_mat;
use lattirust_arithmetic::linear_algebra::Matrix;
use lattirust_arithmetic::linear_algebra::SymmetricMatrix;
use lattirust_arithmetic::nimue::arthur::SerArthur;
use lattirust_arithmetic::nimue::traits::ChallengeFromRandomBytes;
use lattirust_arithmetic::poly_ring::{ConvertibleField, SignedRepresentative};
use relations::traits::Relation;

use crate::util::{BaseRelation, Instance, norm_l2_squared_columnwise, PublicParameters, Witness};

pub fn prove_folding<F: ConvertibleField>(
    arthur: &mut Arthur,
    pp: &PublicParameters<F>,
    instance_1: &Instance<F>,
    witness_1: &Witness<F>,
    instance_2: &Instance<F>,
    witness_2: &Witness<F>,
) -> Result<Matrix<F>, IOPatternError> {
    debug_assert_eq!(witness_1.ncols(), pp.security_parameter);
    debug_assert_eq!(witness_2.ncols(), pp.security_parameter);

    debug_assert!(BaseRelation::is_satisfied(pp, instance_1, witness_1));
    debug_assert!(BaseRelation::is_satisfied(pp, instance_2, witness_2));

    let mut witness = witness_1.clone();
    witness.extend(witness_2.column_iter());

    // Decompose witness matrix into wider matrix with small (column-wise) norms
    let decomp_witness =
        decompose_matrix(&witness, pp.decomposition_basis, pp.decomposition_length);
    debug_assert_eq!(decomp_witness.nrows(), witness.nrows());
    debug_assert_eq!(
        decomp_witness.ncols(),
        witness.ncols() * pp.decomposition_length
    );

    // Commit to decomposed witness matrix (mod q)
    let committed_decomp_witness = &pp.commitment_mat * &decomp_witness;

    // Add to FS transcript
    arthur.absorb_matrix::<F>(&committed_decomp_witness)?;
    arthur.ratchet()?;

    // Compute inner products over the integers
    let decomp_witness_int = &decomp_witness.map(|x| Into::<SignedRepresentative>::into(x).0);
    let inner_products_decomp: SymmetricMatrix<i128> =
        inner_products_mat(&decomp_witness_int).into();
    debug_assert_eq!(
        inner_products_decomp.size(),
        2 * pp.security_parameter * pp.decomposition_length
    );

    // Add to FS transcript
    arthur.absorb_symmetric_matrix::<i128>(&inner_products_decomp)?;
    arthur.ratchet()?;

    // Get challenge
    let challenge = arthur.challenge_matrix::<Trit, TernaryChallengeSet<Trit>>(
        2 * pp.security_parameter * pp.decomposition_length,
        pp.security_parameter,
    )?;
    arthur.ratchet()?;

    // Compute next witness
    // TODO: We don't actually have to do this mod q, but with the current implementation we're barely doing modular arithmetic, not sure if it makes sense to work over the integers instead here.
    let new_witness = mul_f_trit(&decomp_witness, &challenge);

    // Check that the new witness does not grow in norm
    let new_norms = norm_l2_squared_columnwise(&new_witness);
    debug!(
        "Columns of folded witness have norms (min, mean, max) = ({}, {}, {})",
        new_norms.iter().min().unwrap(),
        new_norms.iter().sum::<u128>() as f64 / new_norms.len() as f64,
        new_norms.iter().max().unwrap()
    );
    debug!("(Assuming uniform witness norms), expected norm should be should give (mean, max) = ({}, {})",
        (2 * pp.decomposition_length * pp.security_parameter * pp.decomposition_basis as usize) as f64 / 3.,
        (2 * pp.decomposition_length * pp.security_parameter * pp.decomposition_basis as usize) as f64);
    debug_assert!(new_norms
        .into_iter()
        .all(|l2_norm_sq| l2_norm_sq as f64 <= pp.norm_bound));

    Ok(new_witness)
}
