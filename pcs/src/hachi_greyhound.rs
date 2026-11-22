//! Hachi's final evaluation over the quartic field, authenticated through
//! Greyhound. Use the generic trace-dual packing of §3.1 without peeling
//! coordinates. Do not use the fixed-basis packing of extension-valued
//! partial evaluations in §4.5: it does not bind the original evaluation.

use crate::field::fq;
use ark_ff::{AdditiveGroup, Field};
use lattirust_arithmetic::ring::{PolyRing, Pow2CyclotomicPolyRing};

use crate::field::{Fq, Fq4};
use crate::greyhound::{Claim, Evaluation, GreyhoundKey, Projection};
use crate::hachi::FinalOpening;

pub type R64 = Pow2CyclotomicPolyRing<Fq, 64>;
const D: usize = 64;
const K: usize = 4;
const WIDTH: usize = D / K;

/// Parameters and indices are invalid, or a claimed evaluation is inconsistent.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum AdapterError {
    #[error("wrong Hachi successor table or evaluation dimensions")]
    Dimensions,
    #[error("inconsistent Hachi successor evaluation")]
    Evaluation,
}

/// The basis Z_i of the real quartic subfield. The conjugations defining the
/// trace fix these elements. For i>0, Z_i = X^(8i) + X^(-8i).
fn basis(i: usize) -> R64 {
    if i == 0 {
        return R64::from(Fq::ONE);
    }
    let n = (D / (2 * K)) * i;
    R64::from_fn(|j| match j {
        j if j == n => Fq::ONE,
        j if j == D - n => -Fq::ONE,
        _ => Fq::ZERO,
    })
}

/// Embed the arkworks tower in the precisely chosen real subfield of R64.
pub(crate) fn embed(a: Fq4) -> R64 {
    let s = basis(2); // s² = 2
    const A: Fq = fq(1492342050);
    const I: Fq = fq(3311696422);
    // T = A (Z1 + I Z3). I²=-1 and 2(1+I)A²=1, hence T²=s.
    let t = R64::from(A) * (basis(1) + R64::from(I) * basis(3));
    (R64::from(a.c0.c0) + s * R64::from(a.c0.c1))
        + t * (R64::from(a.c1.c0) + s * R64::from(a.c1.c1))
}

fn monomial(i: usize) -> R64 {
    R64::from_fn(|j| if i == j { Fq::ONE } else { Fq::ZERO })
}

/// Eq. (8): shift first-half field elements by X^i and second-half elements
/// by X^(D/2+i). The input order is the Boolean suffix order, MSB first.
fn psi(input: &[R64]) -> R64 {
    debug_assert_eq!(input.len(), WIDTH);
    let mut out = R64::from(Fq::ZERO);
    let half = WIDTH / 2;
    for (i, item) in input.iter().enumerate() {
        let shift = if i < half { i } else { D / 2 + i - half };
        out += *item * monomial(shift);
    }
    out
}

/// Pack the exact successor coefficients to which the Greyhound key commits.
/// Each base coefficient is a base scalar in the quartic field. A ring word
/// packs D/K=16 source coefficients, four times fewer than in the §4.5
/// specialization; no extension-valued partial evaluations are combined.
pub fn pack_successor(coefficients: &[Fq]) -> Result<Vec<R64>, AdapterError> {
    if coefficients.len() < WIDTH || !coefficients.len().is_power_of_two() {
        return Err(AdapterError::Dimensions);
    }
    Ok(coefficients
        .chunks_exact(WIDTH)
        .map(|word| {
            let mut out = [Fq::ZERO; D];
            for (i, &value) in word.iter().enumerate() {
                let shift = if i < WIDTH / 2 {
                    i
                } else {
                    D / 2 + i - WIDTH / 2
                };
                out[shift] = value;
            }
            R64::from(out)
        })
        .collect())
}

fn conjugate(input: R64) -> R64 {
    let a = input.coefficients();
    R64::from_fn(|i| if i == 0 { a[0] } else { -a[D - i] })
}

/// Normalize the subgroup trace rather than absorbing D/K into every RHS.
fn trace(input: R64) -> R64 {
    let c = input.coefficients();
    let step = D / (2 * K);
    let coeff: [Fq; D] = std::array::from_fn(|i| {
        if i == 0 {
            c[0]
        } else if i % step == 0 {
            // The cyclic subgroup projects onto multiples of D/(2K).
            // Average the projected polynomial with its conjugate.
            (c[i] - c[D - i]) * fq(2147483599)
        } else {
            Fq::ZERO
        }
    });
    R64::from(coeff)
}

fn eq_weights(point: &[Fq4]) -> Vec<Fq4> {
    let mut out = vec![Fq4::ONE];
    for x in point {
        out = out
            .into_iter()
            .flat_map(|v| [v * (Fq4::ONE - x), v * x])
            .collect();
    }
    out
}

/// Checks the trace of the exact Greyhound ring target. The weights cover
/// the entire original point, and the trace is compared directly to the
/// claimed quartic-field value. This is only a reduction: authenticate the
/// returned ring evaluation under the same prebound commitment.
pub fn reduce_evaluation(
    point: &[Fq4],
    value: Fq4,
    ring_target: R64,
    row_log: usize,
) -> Result<Evaluation<R64>, AdapterError> {
    let (outer, inner) = weights(point, row_log)?;
    if trace(ring_target) != embed(value) {
        return Err(AdapterError::Evaluation);
    }
    Ok(Evaluation {
        outer,
        inner,
        value: ring_target,
        projection: Projection::Ring,
    })
}

fn weights(point: &[Fq4], row_log: usize) -> Result<(Vec<R64>, Vec<R64>), AdapterError> {
    const LOG_WIDTH: usize = 4;
    if point.len() < LOG_WIDTH
        || point.len() >= usize::BITS as usize
        || row_log > point.len() - LOG_WIDTH
    {
        return Err(AdapterError::Dimensions);
    }
    let outer = eq_weights(&point[..row_log])
        .into_iter()
        .map(embed)
        .collect();
    let lower = eq_weights(&point[row_log..]);
    let inner = lower
        .chunks_exact(WIDTH)
        .map(|chunk| {
            let values: [R64; WIDTH] = std::array::from_fn(|i| embed(chunk[i]));
            conjugate(psi(&values))
        })
        .collect();
    Ok((outer, inner))
}

/// The prover's ring evaluation of the packed successor. Its authenticity
/// is verified using Greyhound, not by accepting this value alone.
#[derive(Clone, Debug, PartialEq)]
pub struct AdapterProof {
    pub ring_target: R64,
}

/// Reference prover for the trace adapter, without the Greyhound proof.
pub fn prepare_adapter(
    coefficients: &[Fq],
    point: &[Fq4],
    row_log: usize,
) -> Result<(Fq4, AdapterProof), AdapterError> {
    let domain = u32::try_from(point.len())
        .ok()
        .and_then(|n| 1usize.checked_shl(n))
        .ok_or(AdapterError::Dimensions)?;
    if domain != coefficients.len() {
        return Err(AdapterError::Dimensions);
    }
    let (outer, inner) = weights(point, row_log)?;
    let value = eq_weights(point)
        .into_iter()
        .zip(coefficients)
        .map(|(w, &c)| w * Fq4::from_base_prime_field(c))
        .sum();
    let packed = pack_successor(coefficients)?;
    let mut ring_target = R64::from(Fq::ZERO);
    for (r, a) in packed.chunks_exact(inner.len()).zip(outer) {
        ring_target += a * r.iter().zip(&inner).map(|(&b, &c)| b * c).sum::<R64>();
    }
    if trace(ring_target) != embed(value) {
        return Err(AdapterError::Evaluation);
    }
    Ok((value, AdapterProof { ring_target }))
}

/// Convert the authenticated Hachi endpoint into its still-unproved Greyhound claim.
pub fn reduce_final(
    opening: &FinalOpening<Fq4, Vec<R64>>,
    key: &GreyhoundKey<R64>,
    adapter: &AdapterProof,
) -> Result<Claim<R64>, AdapterError> {
    if !key.rows.is_power_of_two() || key.input_depth != 1 {
        return Err(AdapterError::Dimensions);
    }
    let evaluation = reduce_evaluation(
        opening.point(),
        opening.value(),
        adapter.ring_target,
        key.rows.ilog2() as usize,
    )?;
    Ok(Claim {
        commitment: opening.commitment().clone(),
        evaluations: vec![evaluation],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Fq2;
    use crate::greyhound::{self, OpeningMessage};

    fn e(x: u64, y: u64, z: u64, w: u64) -> Fq4 {
        Fq4::new(Fq2::new(x.into(), y.into()), Fq2::new(z.into(), w.into()))
    }

    #[test]
    fn embedding_preserves_multiplication() {
        for i in 0..10 {
            let a = e(i, i + 13, 12, 9);
            let b = e(24, 0, i * 7, i + 3);
            assert_eq!(embed(a * b), embed(a) * embed(b));
            assert_eq!(trace(embed(a)), embed(a));
        }
    }

    #[test]
    fn unpeeled_trace_matches_direct_mle_and_binds_target() {
        let point = [e(4, 9, 12, 1); 8];
        let table: Vec<Fq> = (0..256).map(|i| Fq::from((i * i + 3) as u64)).collect();
        let expected: Fq4 = eq_weights(&point)
            .into_iter()
            .zip(&table)
            .map(|(w, &t)| w * Fq4::from_base_prime_field(t))
            .sum();
        for split in 0..=4 {
            let (value, p) = prepare_adapter(&table, &point, split).unwrap();
            assert_eq!(expected, value);
            assert!(reduce_evaluation(&point, value, p.ring_target, split).is_ok());
            assert!(reduce_evaluation(&point, value + Fq4::ONE, p.ring_target, split).is_err());
            assert!(
                reduce_evaluation(&point, value, p.ring_target + R64::from(Fq::ONE), split)
                    .is_err()
            );
        }
    }

    #[test]
    fn same_commitment_closes_both_protocol_reductions() {
        use crate::greyhound::{check_disclosed_witness, Error};
        use crate::hachi::{prove_linear_reference, reduce_linear_relation, LinearStatement};
        use ark_ff::PrimeField;

        // Deliberately small, nonsecurity parameters and a disclosed witness.
        // The test checks composition and rejection, not succinctness.
        fn verify_final(
            opening: &FinalOpening<Fq4, Vec<R64>>,
            key: &GreyhoundKey<R64>,
            adapter: &AdapterProof,
            message: &OpeningMessage<R64>,
            witness: &[R64],
        ) -> Result<(), Error> {
            let claim = reduce_final(opening, key, adapter).map_err(|_| Error::Unsatisfied)?;
            let relation = greyhound::reduce(key, &claim, message, opening.context())?;
            check_disclosed_witness(&relation, witness)
                .then_some(())
                .ok_or(Error::Unsatisfied)
        }
        let scalar = |x: u64| R64::from(Fq::from(x));
        let decompose = |value: R64, depth| {
            let c = value.coefficients();
            Ok((0..depth)
                .map(|j| {
                    R64::from_fn(|i| {
                        let v = c[i].into_bigint();
                        Fq::from((v.as_ref()[0] >> (4 * j)) & 15)
                    })
                })
                .collect())
        };
        let key = GreyhoundKey {
            rows: 2,
            cols: 8,
            radix: 16,
            depth: 8,
            input_depth: 1,
            challenge_weight: 3,
            norm_squared: 2_000_000,
            a: vec![(3..=10).map(scalar).collect()],
            b: vec![(1..=16).map(scalar).collect()],
            d: vec![(17..=32).map(scalar).collect()],
        };
        let statement = LinearStatement {
            ring_degree: 8,
            base: 4,
            matrix: vec![vec![vec![Fq::ZERO, Fq::ONE]]],
            rhs: vec![vec![-Fq::ONE]],
        };
        let mut coefficients = vec![Fq::ZERO; 256];
        coefficients[7] = Fq::ONE; // X^7 * X = -1 + (X^8 + 1)
        coefficients[8] = Fq::ONE; // Quotient in its first radix digit.
        let table: Vec<_> = pack_successor(&coefficients)
            .unwrap()
            .chunks_exact(8)
            .map(|row| row.to_vec())
            .collect();
        let root = greyhound::commit(&key, &table, &decompose).unwrap();
        let ctx = b"same commitment integration fixture";
        let hachi =
            prove_linear_reference::<Fq4, _>(ctx, &statement, &coefficients, root.clone()).unwrap();
        let request = reduce_linear_relation(ctx, &statement, &hachi).unwrap();
        let (value, adapter) = prepare_adapter(&coefficients, request.point(), 1).unwrap();
        assert_eq!(request.value(), value);
        let evaluation = reduce_evaluation(request.point(), value, adapter.ring_target, 1).unwrap();
        let claim = Claim {
            commitment: root,
            evaluations: vec![evaluation],
        };
        let proof = greyhound::prepare_opening(&key, &claim, &table, &decompose, request.context())
            .unwrap();
        assert!(verify_final(
            &request,
            &key,
            &adapter,
            &proof.message,
            &proof.short_witness
        )
        .is_ok());

        let mut wrong_message = proof.message.clone();
        wrong_message.auxiliary_commitments[0][0] += scalar(1);
        assert!(verify_final(
            &request,
            &key,
            &adapter,
            &wrong_message,
            &proof.short_witness
        )
        .is_err());
        let mut wrong_adapter = adapter;
        wrong_adapter.ring_target += scalar(1);
        assert!(verify_final(
            &request,
            &key,
            &wrong_adapter,
            &proof.message,
            &proof.short_witness
        )
        .is_err());
    }
}
