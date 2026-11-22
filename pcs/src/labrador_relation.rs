//! Convert Greyhound's short affine constraints to lattirust's native
//! LaBRADOR principal relation. This is not the LaBRADOR verifier itself.

use lattirust_arithmetic::{linear_algebra::Vector, ring::PolyRing};
use relations::principal_relation::{
    ConstantQuadraticConstraint, Index, Instance, QuadraticConstraint, Size,
};

use crate::greyhound::{Projection, ShortRelation};

#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum RelationError {
    #[error("invalid Greyhound row or LaBRADOR witness layout")]
    Shape,
    #[error("integer norm bound is too large for the LaBRADOR index")]
    Norm,
}

/// Converts all rows and preserves the squared norm bound exactly in the
/// integer range where the backend's f64 index and comparisons are exact
/// near the threshold. Unused padding contributes to the squared norm
/// and can be removed from an extracted witness.
pub fn convert<R: PolyRing>(
    relation: &ShortRelation<R>,
    words_per_vector: usize,
) -> Result<(Index<R>, Instance<R>), RelationError> {
    if relation.witness_len == 0 || words_per_vector == 0 {
        return Err(RelationError::Shape);
    }
    if relation.norm_squared == 0 || relation.norm_squared >= (1u64 << 52) {
        return Err(RelationError::Norm);
    }
    let rank = relation.witness_len.div_ceil(words_per_vector);
    let padded = rank
        .checked_mul(words_per_vector)
        .ok_or(RelationError::Shape)?;

    let mut full = Vec::new();
    let mut constant = Vec::new();
    for row in &relation.rows {
        let mut coeffs = vec![R::ZERO; padded];
        for &(i, coeff) in &row.terms {
            if i >= relation.witness_len {
                return Err(RelationError::Shape);
            }
            coeffs[i] += coeff;
        }
        let phi: Vec<Vector<R>> = coeffs
            .chunks_exact(words_per_vector)
            .map(|v| Vector::<R>::from(v.to_vec()))
            .collect();
        match row.projection {
            Projection::Ring => full.push(QuadraticConstraint::new_linear(phi, row.value)),
            Projection::ConstantCoefficient => {
                let c = row.value.coefficients();
                if R::from_scalar(c[0]) != row.value {
                    return Err(RelationError::Shape);
                }
                constant.push(ConstantQuadraticConstraint::new_linear(phi, c[0]));
            }
        }
    }
    let size = Size {
        num_witnesses: rank,
        witness_len: words_per_vector,
        norm_bound_sq: relation.norm_squared as f64,
        num_constraints: full.len(),
        num_constant_constraints: constant.len(),
    };
    Ok((
        Index::new(&size),
        Instance {
            quad_dot_prod_funcs: full,
            ct_quad_dot_prod_funcs: constant,
        },
    ))
}

#[cfg(test)]
mod tests {
    use lattirust_arithmetic::ring::Pow2CyclotomicPolyRing;
    use relations::{
        principal_relation::{PrincipalRelation, Witness},
        Relation,
    };

    use super::*;
    use crate::{field::Fq, greyhound::Row};

    type R = Pow2CyclotomicPolyRing<Fq, 8>;

    #[test]
    fn full_and_constant_rows_survive_conversion_and_padding() {
        let one = R::from(Fq::from(1u64));
        let two = one + one;
        let x = R::from_fn(|i| Fq::from(u64::from(i == 1)));
        let rows = vec![
            Row {
                terms: vec![(0, one), (2, one), (2, one)],
                value: one + two,
                projection: Projection::Ring,
            },
            Row {
                terms: vec![(1, one), (2, one)],
                value: one,
                projection: Projection::ConstantCoefficient,
            },
        ];
        let relation = ShortRelation {
            witness_len: 3,
            rows,
            norm_squared: 3,
            fold: vec![],
            transcript_digest: [0u8; 32],
        };
        let (index, instance) = convert(&relation, 2).unwrap();
        let witness = Witness::new(vec![
            Vector::<R>::from(vec![one, x]),
            Vector::<R>::from(vec![one, R::from(Fq::from(0u64))]),
        ]);
        assert!(PrincipalRelation::<R>::is_satisfied(
            &index, &instance, &witness
        ));
        let wrong = Witness::new(vec![
            Vector::<R>::from(vec![one, x]),
            Vector::<R>::from(vec![two, R::from(Fq::from(0u64))]),
        ]);
        assert!(!PrincipalRelation::<R>::is_satisfied(
            &index, &instance, &wrong
        ));
    }
}
