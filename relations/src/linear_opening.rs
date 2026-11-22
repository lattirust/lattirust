//! Linear constraints for split module-SIS opening equations.
//! Builds principal linear constraints for split module-SIS responses.
//! Transcript binding, Gaussian sampling, and response norm bounds belong to
//! the caller. This module has no protocol transcript or witness cache.
use crate::principal_relation::QuadraticConstraint;
use lattirust_arithmetic::linear_algebra::{Matrix, Vector};
use lattirust_arithmetic::ring::{PolyRing, Ring};
use std::ops::Range;

/// Invalid matrix, range, or plane dimensions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Shape,
}

pub struct ResponseBlock<'a, R: PolyRing> {
    pub key: &'a Matrix<R>,
    pub low: Range<usize>,
    pub high: Range<usize>,
    pub high_weight: R::BaseRing,
}

pub struct OpeningPlanes {
    pub announcement: Range<usize>,
    pub commitment: Range<usize>,
    pub width: usize,
}

fn place<R: PolyRing>(
    row: &mut [R],
    at: &Range<usize>,
    coefficients: impl Iterator<Item = R>,
) -> Result<(), Error> {
    let segment = row.get_mut(at.clone()).ok_or(Error::Shape)?;
    for (slot, value) in segment.iter_mut().zip(coefficients) {
        *slot += value;
    }
    Ok(())
}

/// Compile A(z_lo + B*z_hi) = decode(w) + c*decode(t).
/// The middle values w and t use signed radix-2 planes. Their corresponding
/// outer opening constraints have to use exactly these plane coordinates.
/// The caller must enforce response norm constraints on z_lo and z_hi separately.
pub fn compile<R: PolyRing>(
    count: usize,
    responses: &[ResponseBlock<'_, R>],
    planes: &OpeningPlanes,
    challenge: R,
) -> Result<Vec<QuadraticConstraint<R>>, Error> {
    let rank = responses.first().ok_or(Error::Shape)?.key.nrows();
    let width = planes.width;
    if rank == 0
        || width == 0
        || rank.checked_mul(width) != Some(planes.announcement.len())
        || planes.announcement.len() != planes.commitment.len()
        || responses.iter().any(|b| {
            b.key.nrows() != rank || b.key.ncols() != b.low.len() || b.low.len() != b.high.len()
        })
    {
        return Err(Error::Shape);
    }
    let mut equations = Vec::with_capacity(rank);
    for i in 0..rank {
        let mut row = vec![R::ZERO; count];
        for block in responses {
            place(
                &mut row,
                &block.low,
                (0..block.key.ncols()).map(|j| block.key[(i, j)]),
            )?;
            place(
                &mut row,
                &block.high,
                (0..block.key.ncols()).map(|j| block.key[(i, j)] * block.high_weight),
            )?;
        }
        let mut radix = R::BaseRing::ONE;
        for j in 0..width {
            let offset = j * rank + i;
            let negative = if j + 1 == width { radix } else { -radix };
            for (range, factor) in [
                (&planes.announcement, R::ONE),
                (&planes.commitment, challenge),
            ] {
                let position = range.start.checked_add(offset).ok_or(Error::Shape)?;
                *row.get_mut(position).ok_or(Error::Shape)? += factor * negative;
            }
            radix += radix;
        }
        equations.push(QuadraticConstraint::new_homogeneous_linear(vec![
            Vector::from_vec(row),
        ]));
    }
    Ok(equations)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::principal_relation::Witness;
    use lattirust_arithmetic::ring::{Pow2CyclotomicPolyRing, Zq1};
    type F = Zq1<3329>;
    type P = Pow2CyclotomicPolyRing<F, 8>;

    fn signed_scalar(value: i64) -> P {
        let x = F::try_from(value.unsigned_abs()).unwrap();
        P::from_scalar(if value < 0 { -x } else { x })
    }

    #[test]
    fn split_opening_equation_matches_independent_scalar_identity() {
        let key = Matrix::from_fn(1, 2, |_, j| signed_scalar(7 + j as i64));
        let block = ResponseBlock {
            key: &key,
            low: 0..2,
            high: 2..4,
            high_weight: F::try_from(17u64).unwrap(),
        };
        let planes = OpeningPlanes {
            announcement: 4..16,
            commitment: 16..28,
            width: 12,
        };
        let constraints = compile(28, &[block], &planes, signed_scalar(7)).unwrap();
        // 7*(3+17*1) + 8*(7-17*2) = -76 = -181 + 7*15.
        let mut words = vec![
            signed_scalar(3),
            signed_scalar(7),
            signed_scalar(1),
            signed_scalar(-2),
        ];
        for value in [-181i64, 15] {
            let encoded = value.rem_euclid(1 << 12) as u64;
            words.extend((0..12).map(|bit| signed_scalar(((encoded >> bit) & 1) as i64)));
        }
        let eval =
            |words: &[P]| constraints[0].eval(&Witness::new(vec![Vector::from_slice(words)]));
        assert_eq!(eval(&words), P::ZERO);
        for position in [0, 2, 4, 16] {
            words[position] += P::ONE;
            assert_ne!(eval(&words), P::ZERO);
            words[position] -= P::ONE;
        }
    }

    #[test]
    fn rejects_inconsistent_ranges_and_out_of_bounds_coordinates() {
        let key = Matrix::from_fn(1, 2, |_, _| P::ONE);
        let invalid = ResponseBlock {
            key: &key,
            low: 0..1,
            high: 1..3,
            high_weight: F::ONE,
        };
        let planes = OpeningPlanes {
            announcement: 3..5,
            commitment: 5..7,
            width: 2,
        };
        assert!(matches!(
            compile(7, &[invalid], &planes, P::ONE),
            Err(Error::Shape)
        ));
        let valid = ResponseBlock {
            key: &key,
            low: 0..2,
            high: 2..4,
            high_weight: F::ONE,
        };
        assert!(matches!(
            compile(4, &[valid], &planes, P::ONE),
            Err(Error::Shape)
        ));
        assert!(matches!(
            compile::<P>(7, &[], &planes, P::ONE),
            Err(Error::Shape)
        ));
    }
}
