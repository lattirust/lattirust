//! Greyhound's verifier-side reduction to a short affine relation.
//!
//! The reduction is not itself an argument. A verifier must run a LaBRADOR
//! argument for the resulting short relation for the returned statement.
//! This uses an explicit SHAKE256 transcript, not the C implementation's wire
//! format. Parameters and keys must be fixed by the verifier; this module does
//! not calibrate Module-SIS hardness, extraction bounds, or short differences.
use lattirust_arithmetic::ring::{PolyRing, Ring};
use lattirust_arithmetic::traits::WithL2Norm;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};
pub mod grouped;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Error {
    Shape,
    Parameters,
    Serialization,
    Sampling,
    Unsatisfied,
}

/// `A` has `cols*input_depth` columns, `B` has
/// `rows*rank(A)*depth` columns, and `D` has `rows*depth` columns.
/// Set input_depth=1 for a pre-shortened table (e.g. Hachi §4.5).
#[derive(Clone, Debug, PartialEq)]
pub struct GreyhoundKey<R> {
    pub rows: usize,
    pub cols: usize,
    pub radix: u32,
    pub depth: usize,
    pub input_depth: usize,
    pub challenge_weight: usize,
    pub norm_squared: u64,
    pub a: Vec<Vec<R>>,
    pub b: Vec<Vec<R>>,
    pub d: Vec<Vec<R>>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Projection {
    Ring,
    ConstantCoefficient,
}

/// Value of outer^T W inner, optionally projected to the constant coefficient.
/// A projected value must be encoded as a constant polynomial.
#[derive(Clone, Debug, PartialEq)]
pub struct Evaluation<R> {
    pub outer: Vec<R>,
    pub inner: Vec<R>,
    pub value: R,
    pub projection: Projection,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Claim<R> {
    pub commitment: Vec<R>,
    pub evaluations: Vec<Evaluation<R>>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct OpeningMessage<R> {
    /// One D commitment per evaluation, all bound before the common fold.
    pub auxiliary_commitments: Vec<Vec<R>>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Row<R> {
    pub terms: Vec<(usize, R)>,
    pub value: R,
    pub projection: Projection,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ShortRelation<R> {
    pub witness_len: usize,
    pub rows: Vec<Row<R>>,
    pub norm_squared: u64,
    pub fold: Vec<R>,
    /// Domain-separated digest binding parameters, public claim and message.
    pub transcript_digest: [u8; 32],
}

/// Deliberately non-succinct check for fixtures and debugging, not a PCS proof.
pub fn check_disclosed_witness<R: PolyRing>(relation: &ShortRelation<R>, witness: &[R]) -> bool {
    witness.len() == relation.witness_len
        && witness.l2_norm_squared() <= relation.norm_squared.into()
        && relation.rows.iter().all(|row| {
            let sum = row.terms.iter().try_fold(R::ZERO, |sum, &(i, coef)| {
                Some(sum + coef * *witness.get(i)?)
            });
            sum.is_some_and(|sum| project(sum, row.projection) == row.value)
        })
}

fn project<R: PolyRing>(x: R, projection: Projection) -> R {
    match projection {
        Projection::Ring => x,
        Projection::ConstantCoefficient => R::from_scalar(
            x.coefficients()
                .first()
                .copied()
                .unwrap_or(R::BaseRing::ZERO),
        ),
    }
}

fn number(t: &mut Shake256, n: usize) -> Result<(), Error> {
    t.update(
        &u64::try_from(n)
            .map_err(|_| Error::Parameters)?
            .to_le_bytes(),
    );
    Ok(())
}

fn bytes(t: &mut Shake256, s: &[u8]) -> Result<(), Error> {
    number(t, s.len())?;
    t.update(s);
    Ok(())
}

fn values<R: Ring>(t: &mut Shake256, vals: &[R]) -> Result<(), Error> {
    number(t, vals.len())?;
    let mut buf = Vec::new();
    for x in vals {
        buf.clear();
        x.serialize_compressed(&mut buf)
            .map_err(|_| Error::Serialization)?;
        bytes(t, &buf)?;
    }
    Ok(())
}

fn matrix<R: Ring>(t: &mut Shake256, rows: &[Vec<R>]) -> Result<(), Error> {
    number(t, rows.len())?;
    for row in rows {
        values(t, row)?;
    }
    Ok(())
}

fn checked_mul(a: usize, b: usize) -> Result<usize, Error> {
    a.checked_mul(b).ok_or(Error::Parameters)
}

fn radix_powers<R: PolyRing>(radix: u32, count: usize) -> Result<Vec<R>, Error> {
    let radix = R::BaseRing::try_from(radix).map_err(|_| Error::Parameters)?;
    let mut r = R::BaseRing::ONE;
    Ok((0..count)
        .map(|_| {
            let out = r;
            r *= radix;
            R::from_scalar(out)
        })
        .collect())
}

// Uniform over [0,n). Retain a bounded, explicitly failing rejection path.
fn uniform(reader: &mut impl XofReader, n: u64) -> Result<usize, Error> {
    let threshold = n.wrapping_neg() % n;
    let mut bytes = [0u8; 32];
    reader.read(&mut bytes);
    for chunk in bytes.chunks_exact(8) {
        let x = u64::from_le_bytes(chunk.try_into().map_err(|_| Error::Sampling)?);
        if x >= threshold {
            return usize::try_from(x % n).map_err(|_| Error::Sampling);
        }
    }
    Err(Error::Sampling)
}

fn sparse<R: PolyRing>(reader: &mut impl XofReader, weight: usize) -> Result<R, Error> {
    let dim = R::dimension();
    let mut pool: Vec<_> = (0..dim).collect();
    let mut coefficients = vec![R::BaseRing::ZERO; dim];
    for i in 0..weight {
        let n = u64::try_from(dim - i).map_err(|_| Error::Parameters)?;
        let j = i + uniform(reader, n)?;
        pool.swap(i, j);
        let mut sign = [0u8; 1];
        reader.read(&mut sign);
        coefficients[pool[i]] = if sign[0] & 1 == 0 {
            R::BaseRing::ONE
        } else {
            -R::BaseRing::ONE
        };
    }
    R::try_from_coefficients(&coefficients).ok_or(Error::Parameters)
}

/// Construct the full short relation for Fig. 1 (one evaluation) or batched
/// separable evaluations (one auxiliary commitment per query and common fold).
/// The expected claim is supplied by the application, never read from the proof.
/// Check public geometry without hashing, sampling or allocating relation rows.
pub(crate) fn validate_claim<R: PolyRing>(
    key: &GreyhoundKey<R>,
    expected: &Claim<R>,
) -> Result<(), Error> {
    let k = key;
    if k.rows == 0
        || k.cols == 0
        || k.radix < 2
        || k.depth == 0
        || k.input_depth == 0
        || k.norm_squared == 0
        || k.challenge_weight == 0
        || k.challenge_weight > R::dimension()
        || !R::dimension().is_power_of_two()
        || k.a.is_empty()
        || k.b.is_empty()
        || k.d.is_empty()
        || expected.evaluations.is_empty()
    {
        return Err(Error::Parameters);
    }
    let z_len = checked_mul(k.cols, k.input_depth)?;
    let t_block = checked_mul(k.a.len(), k.depth)?;
    let t_len = checked_mul(k.rows, t_block)?;
    let w_len = checked_mul(k.rows, k.depth)?;
    let w_total = checked_mul(w_len, expected.evaluations.len())?;
    let w_start = z_len.checked_add(t_len).ok_or(Error::Parameters)?;
    let _witness_len = w_start.checked_add(w_total).ok_or(Error::Parameters)?;
    if k.a.iter().any(|r| r.len() != z_len)
        || k.b.iter().any(|r| r.len() != t_len)
        || k.d.iter().any(|r| r.len() != w_len)
        || expected.commitment.len() != k.b.len()
        || expected.evaluations.iter().any(|e| {
            e.outer.len() != k.rows
                || e.inner.len() != k.cols
                || project(e.value, e.projection) != e.value
        })
    {
        return Err(Error::Shape);
    }
    Ok(())
}

pub fn reduce<R: PolyRing>(
    key: &GreyhoundKey<R>,
    expected: &Claim<R>,
    message: &OpeningMessage<R>,
    context: &[u8],
) -> Result<ShortRelation<R>, Error> {
    let k = key;
    validate_claim(k, expected)?;
    let z_len = checked_mul(k.cols, k.input_depth)?;
    let t_block = checked_mul(k.a.len(), k.depth)?;
    let t_len = checked_mul(k.rows, t_block)?;
    let w_len = checked_mul(k.rows, k.depth)?;
    let w_start = z_len.checked_add(t_len).ok_or(Error::Parameters)?;
    let witness_len = w_start
        .checked_add(checked_mul(w_len, expected.evaluations.len())?)
        .ok_or(Error::Parameters)?;
    if message.auxiliary_commitments.len() != expected.evaluations.len()
        || message
            .auxiliary_commitments
            .iter()
            .any(|v| v.len() != k.d.len())
    {
        return Err(Error::Shape);
    }
    let mut transcript = Shake256::default();
    bytes(&mut transcript, b"lattirust/greyhound/short-relation/v2")?;
    bytes(&mut transcript, context)?;
    bytes(&mut transcript, &R::modulus().to_bytes_le())?;
    for n in [
        R::dimension(),
        k.rows,
        k.cols,
        k.depth,
        k.input_depth,
        k.challenge_weight,
    ] {
        number(&mut transcript, n)?;
    }
    transcript.update(&k.radix.to_le_bytes());
    transcript.update(&k.norm_squared.to_le_bytes());
    matrix(&mut transcript, &k.a)?;
    matrix(&mut transcript, &k.b)?;
    matrix(&mut transcript, &k.d)?;
    values(&mut transcript, &expected.commitment)?;
    number(&mut transcript, expected.evaluations.len())?;
    for e in &expected.evaluations {
        transcript.update(&[match e.projection {
            Projection::Ring => 0,
            Projection::ConstantCoefficient => 1,
        }]);
        values(&mut transcript, &e.outer)?;
        values(&mut transcript, &e.inner)?;
        values(&mut transcript, &[e.value])?;
    }
    matrix(&mut transcript, &message.auxiliary_commitments)?;
    let mut reader = transcript.finalize_xof();
    let mut transcript_digest = [0u8; 32];
    reader.read(&mut transcript_digest);
    let fold = (0..k.rows)
        .map(|_| sparse::<R>(&mut reader, k.challenge_weight))
        .collect::<Result<Vec<_>, _>>()?;

    let g = radix_powers::<R>(k.radix, k.depth)?;
    let g_in = radix_powers::<R>(k.radix, k.input_depth)?;
    let mut rows = Vec::new();
    // Unknowns ordered [z, flattened t_hat, w_hat for each query].
    // B t_hat = u.
    for (b, &value) in k.b.iter().zip(&expected.commitment) {
        rows.push(Row {
            terms: b.iter().enumerate().map(|(i, &x)| (z_len + i, x)).collect(),
            value,
            projection: Projection::Ring,
        });
    }
    // A z = sum_i c_i Gadget(t_hat_i).
    for (a_row, a) in k.a.iter().enumerate() {
        let mut terms: Vec<_> = a.iter().copied().enumerate().collect();
        for (i, &c) in fold.iter().enumerate() {
            for (digit, &scale) in g.iter().enumerate() {
                terms.push((z_len + i * t_block + a_row * k.depth + digit, -(c * scale)));
            }
        }
        rows.push(Row {
            terms,
            value: R::ZERO,
            projection: Projection::Ring,
        });
    }
    for (query, (e, v)) in expected
        .evaluations
        .iter()
        .zip(&message.auxiliary_commitments)
        .enumerate()
    {
        let offset = w_start + query * w_len;
        // D w_hat = v.
        for (d, &value) in k.d.iter().zip(v) {
            rows.push(Row {
                terms: d
                    .iter()
                    .enumerate()
                    .map(|(i, &x)| (offset + i, x))
                    .collect(),
                value,
                projection: Projection::Ring,
            });
        }
        // Outer-weighted evaluation of un-folded w.
        let weighted = |coeffs: &[R]| {
            coeffs
                .iter()
                .enumerate()
                .flat_map(|(i, &w)| {
                    g.iter()
                        .enumerate()
                        .map(move |(j, &scale)| (offset + i * k.depth + j, w * scale))
                })
                .collect::<Vec<_>>()
        };
        rows.push(Row {
            terms: weighted(&e.outer),
            value: e.value,
            projection: e.projection,
        });
        // Inner-weighted folded z = sum_i c_i Gadget(w_hat_i).
        let mut terms = Vec::with_capacity(z_len);
        for (i, &weight) in e.inner.iter().enumerate() {
            for (j, &scale) in g_in.iter().enumerate() {
                terms.push((i * k.input_depth + j, weight * scale));
            }
        }
        terms.extend(weighted(&fold).into_iter().map(|(i, x)| (i, -x)));
        rows.push(Row {
            terms,
            value: R::ZERO,
            projection: Projection::Ring,
        });
    }
    Ok(ShortRelation {
        witness_len,
        rows,
        norm_squared: k.norm_squared,
        fold,
        transcript_digest,
    })
}

/// Private fold output for the next reduction, not a succinct opening proof.
pub struct PreparedOpening<R> {
    pub message: OpeningMessage<R>,
    pub relation: ShortRelation<R>,
    pub short_witness: Vec<R>,
}

fn dot<R: Ring>(a: &[R], b: &[R]) -> R {
    debug_assert_eq!(a.len(), b.len());
    a.iter().zip(b).map(|(&x, &y)| x * y).sum()
}

fn digits<R: PolyRing>(
    radix: u32,
    depth: usize,
    value: R,
    decompose: &impl Fn(R, usize) -> Result<Vec<R>, Error>,
) -> Result<Vec<R>, Error> {
    let digits = decompose(value, depth)?;
    if digits.len() != depth || dot(&digits, &radix_powers::<R>(radix, depth)?) != value {
        return Err(Error::Unsatisfied);
    }
    Ok(digits)
}

fn inner_commitments<R: PolyRing>(
    key: &GreyhoundKey<R>,
    table: &[Vec<R>],
    decompose: &impl Fn(R, usize) -> Result<Vec<R>, Error>,
) -> Result<(Vec<Vec<R>>, Vec<R>), Error> {
    if key.rows == 0
        || key.cols == 0
        || key.radix < 2
        || key.depth == 0
        || key.input_depth == 0
        || key.a.is_empty()
        || key.b.is_empty()
        || table.len() != key.rows
        || table.iter().any(|r| r.len() != key.cols)
        || key
            .a
            .iter()
            .any(|r| r.len() != key.cols.checked_mul(key.input_depth).unwrap_or(0))
        || key.b.iter().any(|r| {
            r.len()
                != key
                    .rows
                    .checked_mul(key.a.len())
                    .and_then(|n| n.checked_mul(key.depth))
                    .unwrap_or(0)
        })
    {
        return Err(Error::Shape);
    }
    let mut source = Vec::new();
    let mut t_hat = Vec::new();
    for row in table {
        let mut s = Vec::new();
        for &x in row {
            if key.input_depth == 1 {
                s.push(x);
            } else {
                s.extend(digits(key.radix, key.input_depth, x, decompose)?);
            }
        }
        for a in &key.a {
            t_hat.extend(digits(key.radix, key.depth, dot(a, &s), decompose)?);
        }
        source.push(s);
    }
    Ok((source, t_hat))
}

/// Recompute the standard inner/outer commitment. `decompose` defines the
/// parameter set's deterministic digit convention and must return short digits.
pub fn commit<R: PolyRing>(
    key: &GreyhoundKey<R>,
    table: &[Vec<R>],
    decompose: &impl Fn(R, usize) -> Result<Vec<R>, Error>,
) -> Result<Vec<R>, Error> {
    let (_, t) = inner_commitments(key, table, decompose)?;
    Ok(key.b.iter().map(|b| dot(b, &t)).collect())
}

/// Produce the honest short witness for the next reduction.
/// Verifier acceptance still requires an argument for the reduced relation.
pub fn prepare_opening<R: PolyRing>(
    key: &GreyhoundKey<R>,
    claim: &Claim<R>,
    table: &[Vec<R>],
    decompose: &impl Fn(R, usize) -> Result<Vec<R>, Error>,
    context: &[u8],
) -> Result<PreparedOpening<R>, Error> {
    let (source, t) = inner_commitments(key, table, decompose)?;
    let w_len = checked_mul(key.rows, key.depth)?;
    if key.d.is_empty()
        || key.d.iter().any(|d| d.len() != w_len)
        || claim
            .evaluations
            .iter()
            .any(|e| e.inner.len() != key.cols || e.outer.len() != key.rows)
    {
        return Err(Error::Shape);
    }
    let mut message = OpeningMessage {
        auxiliary_commitments: Vec::new(),
    };
    let mut witnesses = Vec::new();
    for e in &claim.evaluations {
        let mut w = Vec::new();
        for row in table {
            w.extend(digits(key.radix, key.depth, dot(&e.inner, row), decompose)?);
        }
        message
            .auxiliary_commitments
            .push(key.d.iter().map(|d| dot(d, &w)).collect());
        witnesses.extend(w);
    }
    let relation = reduce(key, claim, &message, context)?;
    let mut short_witness = vec![R::ZERO; checked_mul(key.cols, key.input_depth)?];
    for (src, c) in source.iter().zip(&relation.fold) {
        for (acc, x) in short_witness.iter_mut().zip(src) {
            *acc += *x * *c;
        }
    }
    short_witness.extend(t);
    short_witness.extend(witnesses);
    if !check_disclosed_witness(&relation, &short_witness) {
        return Err(Error::Unsatisfied);
    }
    Ok(PreparedOpening {
        message,
        relation,
        short_witness,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{BigInteger, PrimeField};
    use lattirust_arithmetic::ring::{f_p::Fq, Pow2CyclotomicPolyRing};

    type F = Fq<4294967197>;
    type R = Pow2CyclotomicPolyRing<F, 8>;

    fn r(x: u32) -> R {
        R::from_scalar(F::from(x))
    }

    fn base16(value: R, depth: usize) -> Result<Vec<R>, Error> {
        let c = value.coefficients();
        Ok((0..depth)
            .map(|d| {
                R::from_fn(|i| {
                    let raw = c
                        .get(i)
                        .copied()
                        .unwrap_or(F::ZERO)
                        .into_bigint()
                        .to_bytes_le();
                    F::from(u32::from((raw[d / 2] >> (4 * (d % 2))) & 15))
                })
            })
            .collect())
    }

    fn key() -> GreyhoundKey<R> {
        GreyhoundKey {
            rows: 2,
            cols: 2,
            radix: 16,
            depth: 8,
            input_depth: 1,
            challenge_weight: 2,
            norm_squared: 1_000_000,
            a: vec![vec![r(2), r(3)], vec![r(4), r(5)]],
            b: vec![(0..32).map(|j| r(j + 1)).collect()],
            d: vec![(0..16).map(|j| r(j + 2)).collect()],
        }
    }

    // Test-only backend; it discloses the entire folded witness.
    fn verify(
        key: &GreyhoundKey<R>,
        claim: &Claim<R>,
        message: &OpeningMessage<R>,
        context: &[u8],
        witness: &[R],
    ) -> Result<(), Error> {
        let relation = reduce(key, claim, message, context)?;
        check_disclosed_witness(&relation, witness)
            .then_some(())
            .ok_or(Error::Unsatisfied)
    }

    #[test]
    fn one_and_two_separable_claims_are_bound_and_reduced() {
        let table = vec![vec![r(1), r(2)], vec![r(3), r(4)]];
        let key = key();
        let commitment = commit(&key, &table, &base16).unwrap();
        let query = Evaluation {
            outer: vec![r(2), r(1)],
            inner: vec![r(1), r(3)],
            value: r(29),
            projection: Projection::Ring,
        };
        let mut claim = Claim {
            commitment,
            evaluations: vec![query],
        };
        let open = prepare_opening(&key, &claim, &table, &base16, b"proof/session/v1").unwrap();
        assert_eq!(open.relation.rows.len(), 6);
        assert!(check_disclosed_witness(&open.relation, &open.short_witness));
        claim.evaluations.push(Evaluation {
            outer: vec![r(1), r(2)],
            inner: vec![r(1), r(1)],
            value: r(17),
            projection: Projection::ConstantCoefficient,
        });
        let open = prepare_opening(&key, &claim, &table, &base16, b"proof/session/v1").unwrap();
        assert_eq!(open.relation.rows.len(), 9);
        assert!(verify(
            &key,
            &claim,
            &open.message,
            b"proof/session/v1",
            &open.short_witness
        )
        .is_ok());
        let mut invalid_witness = open.short_witness.clone();
        invalid_witness[0] += r(1);
        assert!(matches!(
            verify(
                &key,
                &claim,
                &open.message,
                b"proof/session/v1",
                &invalid_witness
            ),
            Err(Error::Unsatisfied)
        ));
        assert!(check_disclosed_witness(&open.relation, &open.short_witness));
        let changed = reduce(&key, &claim, &open.message, b"proof/session/v2").unwrap();
        assert_ne!(open.relation.transcript_digest, changed.transcript_digest);
        for c in &changed.fold {
            assert_eq!(
                c.coefficients().iter().filter(|&&x| x != F::ZERO).count(),
                2
            );
        }
        let mut v = open.message.clone();
        v.auxiliary_commitments[1][0] += r(1);
        let changed = reduce(&key, &claim, &v, b"proof/session/v1").unwrap();
        assert_ne!(open.relation.transcript_digest, changed.transcript_digest);
        assert!(!check_disclosed_witness(&changed, &open.short_witness));
        claim.evaluations[1].value += r(1);
        let changed = reduce(&key, &claim, &open.message, b"proof/session/v1").unwrap();
        assert_ne!(open.relation.transcript_digest, changed.transcript_digest);
        assert!(!check_disclosed_witness(&changed, &open.short_witness));
        assert!(matches!(
            prepare_opening(&key, &claim, &table, &base16, b"proof/session/v1"),
            Err(Error::Unsatisfied)
        ));
    }

    #[test]
    fn canonical_digits_and_constant_projection_are_distinct_operations() {
        let mut key = key();
        key.input_depth = 8;
        key.a = vec![
            (0..16).map(|j| r(j + 1)).collect(),
            (0..16).map(|j| r(j + 2)).collect(),
        ];
        let x = R::from_fn(|i| if i == 1 { F::from(1u32) } else { F::ZERO });
        let table = vec![vec![r(19) + x, r(2)], vec![r(3), r(4)]];
        let commitment = commit(&key, &table, &base16).unwrap();
        let mut claim = Claim {
            commitment,
            evaluations: vec![Evaluation {
                outer: vec![r(1), r(0)],
                inner: vec![r(1), r(0)],
                value: r(19),
                projection: Projection::ConstantCoefficient,
            }],
        };
        let open = prepare_opening(&key, &claim, &table, &base16, b"projection-test").unwrap();
        assert!(check_disclosed_witness(&open.relation, &open.short_witness));
        claim.evaluations[0].projection = Projection::Ring;
        assert!(matches!(
            prepare_opening(&key, &claim, &table, &base16, b"projection-test"),
            Err(Error::Unsatisfied)
        ));
        let mut short = open.message;
        short.auxiliary_commitments[0].clear();
        assert!(matches!(
            reduce(&key, &claim, &short, b"projection-test"),
            Err(Error::Shape)
        ));
    }
}
