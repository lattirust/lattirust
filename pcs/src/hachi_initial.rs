//! Hachi §3.1: base-field monomial coefficients to a ring-valued opening.
//!
//! Use the generic trace map, embedding each base coefficient as one field
//! element. This deliberately omits the §3.2 optimization combining four
//! extension-valued partial evaluations with a fixed basis. As flagged by
//! the ArkLib audit, that combination does not authenticate those partials.
//! This choice multiplies the packed word count by four and requires keys and
//! commitments produced for this layout.

use crate::field::fq;
use ark_ff::{AdditiveGroup, Field};
use lattirust_arithmetic::ring::{PolyRing, Pow2CyclotomicPolyRing};

use crate::field::{Fq, Fq4};
use crate::greyhound::{self, Claim, GreyhoundKey, OpeningMessage, ShortRelation};
use crate::hachi::{self, FoldParams, FoldProof, FoldStatement, HachiError};
use crate::hachi_greyhound::embed;
use crate::hachi_greyhound::{self, reduce_final, AdapterError, AdapterProof, R64};

type R<const N: usize> = Pow2CyclotomicPolyRing<Fq, N>;

#[derive(Copy, Clone, Debug, Eq, PartialEq, thiserror::Error)]
pub enum Error {
    #[error("invalid Hachi initial domain or split")]
    Dimensions,
    #[error("invalid Hachi initial ring evaluation trace")]
    Evaluation,
}

fn dimension<const N: usize>(len: usize) -> Result<usize, Error> {
    if N < 8 || !N.is_power_of_two() || !len.is_power_of_two() || len < N / 4 {
        return Err(Error::Dimensions);
    }
    Ok(len / (N / 4))
}

fn monomial_weights(point: &[Fq4]) -> Vec<Fq4> {
    let mut weights = vec![Fq4::ONE];
    for &x in point {
        weights = weights.into_iter().flat_map(|w| [w, w * x]).collect();
    }
    weights
}

fn invariant<const N: usize>(a: Fq4) -> R<N> {
    let standard = embed(a).coefficients();
    R::from_fn(|j| {
        if j % (N / 8) == 0 {
            standard[8 * (j / (N / 8))]
        } else {
            Fq::ZERO
        }
    })
}

/// Eq. (8) using shifts in the negacyclic ring; shares the same coefficient
/// ordering for packing monomials and their evaluation weights.
fn psi<const N: usize>(input: impl IntoIterator<Item = R<N>>) -> R<N> {
    let mut out = [Fq::ZERO; N];
    for (i, rad) in input.into_iter().enumerate() {
        let shift = if i < N / 8 { i } else { N / 2 + i - N / 8 };
        for (j, &coefficient) in rad.coefficients().iter().enumerate() {
            let u = j + shift;
            if u < N {
                out[u] += coefficient;
            } else {
                out[u - N] -= coefficient;
            }
        }
    }
    R::from(out)
}

/// Input uses MSB-first binary indexing of the monomial coefficient table.
/// The last log2(N)-2 bits are packed by the trace-dual map with one base
/// scalar per subfield element; the preceding bits index ring words.
pub fn pack_monomials<const N: usize>(coefficients: &[Fq]) -> Result<Vec<R<N>>, Error> {
    dimension::<N>(coefficients.len())?;
    Ok(coefficients
        .chunks_exact(N / 4)
        .map(|word| {
            // psi([embed(c_i)]_i), specialized to scalar c_i.
            let mut out = [Fq::ZERO; N];
            out[..N / 8].copy_from_slice(&word[..N / 8]);
            out[N / 2..N / 2 + N / 8].copy_from_slice(&word[N / 8..]);
            R::from(out)
        })
        .collect())
}

/// The ring value whose normalized trace authenticates the field value.
/// There are no separate extension-valued partial evaluation claims.
#[derive(Clone, Debug, PartialEq)]
pub struct InitialMessage<const N: usize> {
    pub ring_value: R<N>,
}

fn trace<const N: usize>(a: R<N>) -> R<N> {
    let c = a.coefficients();
    R::from_fn(|i| {
        if i == 0 {
            c[0]
        } else if i % (N / 8) == 0 {
            (c[i] - c[N - i]) * fq(2147483599)
        } else {
            Fq::ZERO
        }
    })
}

fn split<const N: usize>(point: &[Fq4]) -> Result<(usize, usize), Error> {
    let domain = u32::try_from(point.len())
        .ok()
        .and_then(|n| 1usize.checked_shl(n))
        .ok_or(Error::Dimensions)?;
    let words = dimension::<N>(domain)?;
    Ok((point.len() - (N.ilog2() as usize - 2), words))
}

/// Reference prover for §3.1; this does not compute an opening proof.
pub fn prepare<const N: usize>(
    coefficients: &[Fq],
    point: &[Fq4],
) -> Result<(Fq4, InitialMessage<N>), Error> {
    let (end, words) = split::<N>(point)?;
    if (N / 4).checked_mul(words) != Some(coefficients.len()) {
        return Err(Error::Dimensions);
    }
    let packed = pack_monomials::<N>(coefficients)?;
    let value = monomial_weights(point)
        .into_iter()
        .zip(coefficients)
        .map(|(w, &c)| w * Fq4::from_base_prime_field(c))
        .sum();
    let ring_value = monomial_weights(&point[..end])
        .into_iter()
        .zip(packed)
        .map(|(w, f)| invariant::<N>(w) * f)
        .sum();
    Ok((value, InitialMessage { ring_value }))
}

/// Checks the normalized trace, then returns the full ring-valued claim
/// that still must open against the same provided Hachi commitment.
pub fn reduce<const N: usize>(
    commitment: &[R<N>],
    point: &[Fq4],
    value: Fq4,
    message: &InitialMessage<N>,
    row_log: usize,
) -> Result<FoldStatement<Fq>, Error> {
    let (end, words) = split::<N>(point)?;
    if row_log > words.ilog2() as usize || commitment.is_empty() {
        return Err(Error::Dimensions);
    }
    let v = psi(monomial_weights(&point[end..])
        .into_iter()
        .map(invariant::<N>));
    let c = v.coefficients();
    let conjugate = R::from_fn(|i| if i == 0 { c[0] } else { -c[N - i] });
    if trace(message.ring_value * conjugate) != invariant::<N>(value) {
        return Err(Error::Evaluation);
    }
    let weights = |p: &[Fq4]| {
        monomial_weights(p)
            .into_iter()
            .map(|x| invariant::<N>(x).coefficients())
            .collect()
    };
    Ok(FoldStatement {
        outer_commitment: commitment.iter().map(|x| x.coefficients()).collect(),
        value: message.ring_value.coefficients(),
        left: weights(&point[row_log..end]),
        right: weights(&point[..row_log]),
    })
}

/// Public Hachi messages. The result remains a Greyhound opening obligation.
#[derive(Clone, Debug, PartialEq)]
pub struct HachiMessage<const N: usize> {
    pub initial: InitialMessage<N>,
    pub fold: FoldProof<Fq4, Vec<R64>>,
    pub adapter: AdapterProof,
}

#[derive(Debug)]
pub enum ReductionError {
    Initial(Error),
    Fold(HachiError),
    Adapter(AdapterError),
}

/// Check only the Hachi reduction and return the exact successor claim and context.
pub fn reduce_hachi<const N: usize>(
    hachi_key: &FoldParams<Fq>,
    greyhound_key: &GreyhoundKey<R64>,
    commitment: &[R<N>],
    row_log: usize,
    point: &[Fq4],
    value: Fq4,
    context: &[u8],
    message: &HachiMessage<N>,
) -> Result<(Claim<R64>, [u8; 32]), ReductionError> {
    if hachi_key.ring_degree != N {
        return Err(ReductionError::Initial(Error::Dimensions));
    }
    let bound_context =
        bind_statement(commitment, point, value, &message.initial, row_log, context)
            .map_err(ReductionError::Fold)?;
    let statement = reduce(commitment, point, value, &message.initial, row_log)
        .map_err(ReductionError::Initial)?;
    let open = hachi::reduce_ring_opening(&bound_context, hachi_key, &statement, &message.fold)
        .map_err(ReductionError::Fold)?;
    let claim =
        reduce_final(&open, greyhound_key, &message.adapter).map_err(ReductionError::Adapter)?;
    Ok((claim, *open.context()))
}

fn bind_statement<const N: usize>(
    commitment: &[R<N>],
    point: &[Fq4],
    value: Fq4,
    message: &InitialMessage<N>,
    row_log: usize,
    context: &[u8],
) -> Result<[u8; 32], HachiError> {
    // Bind the original field statement and initial reduction messages,
    // not only the induced ring evaluation, before sampling the fold.
    use ark_serialize::CanonicalSerialize;
    use sha3::{
        digest::{ExtendableOutput, Update, XofReader},
        Shake256,
    };
    let mut fs = Shake256::default();
    fs.update(b"lattirust/hachi/base-field-statement/no-peeling/v1");
    fs.update(&(N as u64).to_le_bytes());
    fs.update(&(row_log as u64).to_le_bytes());
    fs.update(&(context.len() as u64).to_le_bytes());
    fs.update(context);
    let mut canonical = Vec::new();
    (|| {
        commitment.serialize_compressed(&mut canonical)?;
        point.serialize_compressed(&mut canonical)?;
        value.serialize_compressed(&mut canonical)?;
        message.ring_value.serialize_compressed(&mut canonical)
    })()
    .map_err(|_| HachiError::Serialization)?;
    fs.update(&canonical);
    let mut bound_context = [0u8; 32];
    fs.finalize_xof().read(&mut bound_context);
    Ok(bound_context)
}

/// Openings waiting on a separate short-relation prover. The short witness
/// is private prover state and must not be serialized as a succinct proof.
pub struct PreparedOpening<const N: usize> {
    pub initial: InitialMessage<N>,
    pub fold: FoldProof<Fq4, Vec<R64>>,
    pub adapter: AdapterProof,
    pub greyhound: OpeningMessage<R64>,
    pub short_relation: ShortRelation<R64>,
    pub short_witness: Vec<R64>,
}

#[derive(Debug)]
pub enum ProverError {
    Initial(Error),
    Fold(HachiError),
    Adapter(AdapterError),
    Greyhound(greyhound::Error),
}

/// Compute the transparent Hachi commitment for base-field monomial
/// coefficients, using verifier-selected matrices.
pub fn commit<const N: usize>(
    key: &FoldParams<Fq>,
    coefficients: &[Fq],
    row_log: usize,
) -> Result<Vec<R<N>>, ProverError> {
    let input = ring_rows::<N>(key, coefficients, row_log)?;
    Ok(hachi::commit_ring_reference(key, &input)
        .map_err(ProverError::Fold)?
        .into_iter()
        .map(R::from)
        .collect())
}

fn ring_rows<const N: usize>(
    key: &FoldParams<Fq>,
    coefficients: &[Fq],
    row_log: usize,
) -> Result<Vec<Vec<Vec<Fq>>>, ProverError> {
    let packed = pack_monomials::<N>(coefficients).map_err(ProverError::Initial)?;
    if key.ring_degree != N || row_log > packed.len().ilog2() as usize {
        return Err(ProverError::Initial(Error::Dimensions));
    }
    Ok(packed
        .chunks_exact(packed.len() >> row_log)
        .map(|row| row.iter().map(PolyRing::coefficients).collect())
        .collect())
}

/// A correctness-first prover for every stage before the succinct backend.
/// Sumcheck is intentionally reference-only; use a specialized prover for
/// large instances.
pub fn prepare_hachi_reference<const N: usize>(
    hachi_key: &FoldParams<Fq>,
    greyhound_key: &GreyhoundKey<R64>,
    commitment: &[R<N>],
    row_log: usize,
    point: &[Fq4],
    coefficients: &[Fq],
    context: &[u8],
    decompose: &impl Fn(R64, usize) -> Result<Vec<R64>, greyhound::Error>,
) -> Result<(Fq4, HachiMessage<N>, Claim<R64>, Vec<Vec<R64>>, [u8; 32]), ProverError> {
    let input = ring_rows::<N>(hachi_key, coefficients, row_log)?;
    let (value, initial) = prepare::<N>(coefficients, point).map_err(ProverError::Initial)?;
    let statement =
        reduce(commitment, point, value, &initial, row_log).map_err(ProverError::Initial)?;
    let ctx = bind_statement(commitment, point, value, &initial, row_log, context)
        .map_err(ProverError::Fold)?;
    let (v, digits) =
        hachi::prepare_fold(&ctx, hachi_key, &statement, &input).map_err(ProverError::Fold)?;
    let packed = hachi_greyhound::pack_successor(&digits).map_err(ProverError::Adapter)?;
    if greyhound_key.input_depth != 1
        || greyhound_key.cols.checked_mul(greyhound_key.rows) != Some(packed.len())
        || !greyhound_key.rows.is_power_of_two()
    {
        return Err(ProverError::Adapter(AdapterError::Dimensions));
    }
    let table: Vec<_> = packed
        .chunks_exact(greyhound_key.cols)
        .map(|r| r.to_vec())
        .collect();
    let t = greyhound::commit(greyhound_key, &table, decompose).map_err(ProverError::Greyhound)?;
    let fold =
        hachi::prove_ring_reference::<Fq4, _>(&ctx, hachi_key, &statement, v, &digits, t.clone())
            .map_err(ProverError::Fold)?;
    let request = hachi::reduce_ring_opening(&ctx, hachi_key, &statement, &fold)
        .map_err(ProverError::Fold)?;
    let rows = greyhound_key.rows.ilog2() as usize;
    let (opened, adapter) = hachi_greyhound::prepare_adapter(&digits, request.point(), rows)
        .map_err(ProverError::Adapter)?;
    let evaluation = hachi_greyhound::reduce_evaluation(
        request.point(),
        request.value(),
        adapter.ring_target,
        rows,
    )
    .map_err(ProverError::Adapter)?;
    debug_assert_eq!(opened, request.value());
    let claim = Claim {
        commitment: t,
        evaluations: vec![evaluation],
    };
    Ok((
        value,
        HachiMessage {
            initial,
            fold,
            adapter,
        },
        claim,
        table,
        *request.context(),
    ))
}

/// Reference composition retained for callers inspecting the private reduced witness.
pub fn prepare_opening_reference<const N: usize>(
    hachi_key: &FoldParams<Fq>,
    greyhound_key: &GreyhoundKey<R64>,
    commitment: &[R<N>],
    row_log: usize,
    point: &[Fq4],
    coefficients: &[Fq],
    context: &[u8],
    decompose: &impl Fn(R64, usize) -> Result<Vec<R64>, greyhound::Error>,
) -> Result<(Fq4, PreparedOpening<N>), ProverError> {
    let (value, message, claim, table, next_context) = prepare_hachi_reference(
        hachi_key,
        greyhound_key,
        commitment,
        row_log,
        point,
        coefficients,
        context,
        decompose,
    )?;
    let closing =
        greyhound::prepare_opening(greyhound_key, &claim, &table, decompose, &next_context)
            .map_err(ProverError::Greyhound)?;
    Ok((
        value,
        PreparedOpening {
            initial: message.initial,
            fold: message.fold,
            adapter: message.adapter,
            greyhound: closing.message,
            short_relation: closing.relation,
            short_witness: closing.short_witness,
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Fq2;

    #[test]
    fn initial_adapter_matches_direct_monomials() {
        const N: usize = 16;
        let x = Fq4::new(
            Fq2::new(4u64.into(), 2u64.into()),
            Fq2::new(8u64.into(), 3u64.into()),
        );
        let point = [x, x + Fq4::ONE, x, x, x, x];
        let coefficients: Vec<_> = (0..64).map(|i| Fq::from((i * i + 2) as u64)).collect();
        let expected: Fq4 = monomial_weights(&point)
            .into_iter()
            .zip(&coefficients)
            .map(|(w, c)| w * Fq4::from_base_prime_field(*c))
            .sum();
        let (value, message) = prepare::<N>(&coefficients, &point).unwrap();
        assert_eq!(value, expected);
        let commitment = [R::from(Fq::ONE)];
        let packed = pack_monomials::<N>(&coefficients).unwrap();
        assert_eq!(packed.len(), 16);
        for row_log in 0..=4 {
            let s = reduce(&commitment, &point, value, &message, row_log).unwrap();
            let left: Vec<R<N>> = s.left.iter().cloned().map(R::from).collect();
            let right: Vec<R<N>> = s.right.iter().cloned().map(R::from).collect();
            let opened: R<N> = packed
                .chunks_exact(left.len())
                .zip(right)
                .map(|(row, b)| b * row.iter().zip(&left).map(|(&f, &a)| f * a).sum::<R<N>>())
                .sum();
            assert_eq!(opened, message.ring_value);
        }
        let mut wrong = message.clone();
        wrong.ring_value += R::from(Fq::ONE);
        assert!(reduce(&commitment, &point, value, &wrong, 1).is_err());
    }
    #[test]
    fn field_opening_rejects_wrong_outer_claim_and_context() {
        use ark_ff::PrimeField;
        const N: usize = 8;
        let s = |v: u64| R64::from(Fq::from(v));
        let key = FoldParams {
            ring_degree: N,
            base: 16,
            folded_digits: 2,
            challenge_weight: 1,
            a: vec![(1..=32).map(|i: u64| vec![Fq::from(i)]).collect()],
            b: vec![(9..=16).map(|i: u64| vec![Fq::from(i)]).collect()],
            d: vec![(17..=24).map(|i: u64| vec![Fq::from(i)]).collect()],
        };
        let gh = GreyhoundKey {
            rows: 4,
            cols: 16,
            radix: 16,
            depth: 8,
            input_depth: 1,
            challenge_weight: 1,
            norm_squared: 10_000_000,
            a: vec![(1..=16).map(s).collect()],
            b: vec![(5..=36).map(s).collect()],
            d: vec![(37..=68).map(s).collect()],
        };
        let d = |value: R64, depth| {
            let a = value.coefficients();
            Ok((0..depth)
                .map(|j| {
                    R64::from_fn(|i| Fq::from((a[i].into_bigint().as_ref()[0] >> (4 * j)) & 15))
                })
                .collect())
        };
        let coefficients: Vec<_> = (0..8).map(|i| Fq::from(i as u64)).collect();
        let x = Fq4::new(
            Fq2::new(4u64.into(), 2u64.into()),
            Fq2::new(8u64.into(), 3u64.into()),
        );
        let point = [x, x + Fq4::ONE, x];
        let commitment = commit::<N>(&key, &coefficients, 0).unwrap();
        let context = b"full field opening correctness fixture";
        let (y, prepared) = prepare_opening_reference(
            &key,
            &gh,
            &commitment,
            0,
            &point,
            &coefficients,
            context,
            &d,
        )
        .unwrap();
        let witness = prepared.short_witness.clone();
        let opening = HachiMessage {
            initial: prepared.initial,
            fold: prepared.fold,
            adapter: prepared.adapter,
        };
        let verify = |key,
                      gh,
                      commitment: &Vec<R<N>>,
                      row_log,
                      point: &[Fq4],
                      value,
                      context: &[u8],
                      opening: &HachiMessage<N>| {
            let (claim, next) =
                reduce_hachi(key, gh, commitment, row_log, point, value, context, opening)
                    .map_err(|_| ())?;
            let relation =
                greyhound::reduce(gh, &claim, &prepared.greyhound, &next).map_err(|_| ())?;
            greyhound::check_disclosed_witness(&relation, &witness)
                .then_some(())
                .ok_or(())
        };
        assert!(verify(&key, &gh, &commitment, 0, &point, y, context, &opening,).is_ok());
        assert!(verify(
            &key,
            &gh,
            &commitment,
            0,
            &point,
            y + Fq4::ONE,
            context,
            &opening,
        )
        .is_err());
        assert!(verify(&key, &gh, &commitment, 0, &point, y, b"wrong", &opening,).is_err());
    }
}
