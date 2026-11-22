//! PCS relation handoffs through LattiRust's existing reduction framework.
//!
//! Hachi and Greyhound retain their domain-separated SHAKE transcripts.
//! Merlin/Arthur transport fixed-shape public messages; they do not replace
//! those challenge derivations. A returned relation is not an accepted proof.
use crate::{
    field::{Fq, Fq4},
    greyhound::{
        self, Claim, Evaluation, GreyhoundKey, OpeningMessage, Projection, Row, ShortRelation,
    },
    hachi::{self, FoldParams, FoldProof, LinearProof, Sumcheck},
    hachi_greyhound::{AdapterProof, R64},
    hachi_initial::{self, HachiMessage, InitialMessage},
    labrador_relation,
};
use anyhow::{ensure, Result};
use ark_ff::{AdditiveGroup, Field, PrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use lattirust_arithmetic::{
    linear_algebra::Vector,
    ring::{PolyRing, Pow2CyclotomicPolyRing, Ring},
};
use nimue::{Arthur, ByteReader, ByteWriter, IOPattern, Merlin, ProofError, ProofResult};
use relations::{
    principal_relation::{Index, Instance, PrincipalRelation, Witness},
    reduction::Reduction,
    Relation,
};
use std::{fmt, marker::PhantomData, sync::Arc};

fn invalid<E>(_: E) -> ProofError {
    ProofError::InvalidProof
}
fn mul(a: usize, b: usize) -> ProofResult<usize> {
    a.checked_mul(b).ok_or(ProofError::InvalidProof)
}
fn add(a: usize, b: usize) -> ProofResult<usize> {
    a.checked_add(b).ok_or(ProofError::InvalidProof)
}

/// Unsigned radix digits of canonical coefficient representatives. This is
/// the convention used by the framework prover; lower-level APIs still accept
/// an explicit decomposition function for other supported conventions.
pub fn canonical_digits<R: PolyRing>(
    value: R,
    radix: u32,
    depth: usize,
) -> Result<Vec<R>, greyhound::Error>
where
    R::BaseRing: PrimeField,
{
    if radix < 2 || depth == 0 {
        return Err(greyhound::Error::Parameters);
    }
    use ark_ff::BigInteger;
    if radix.is_power_of_two() {
        let mut coefficients: Vec<_> = value
            .coefficients()
            .into_iter()
            .map(|x| x.into_bigint())
            .collect();
        let bits = radix.ilog2();
        let mask = u64::from(radix - 1);
        let mut out = Vec::with_capacity(depth);
        for _ in 0..depth {
            let digits: Vec<_> = coefficients
                .iter_mut()
                .map(|x| {
                    let digit = R::BaseRing::from((x.as_ref()[0] & mask) as u32);
                    x.divn(bits);
                    digit
                })
                .collect();
            out.push(R::try_from_coefficients(&digits).ok_or(greyhound::Error::Shape)?);
        }
        if coefficients
            .iter()
            .any(|x| x.as_ref().iter().any(|&limb| limb != 0))
        {
            return Err(greyhound::Error::Unsatisfied);
        }
        return Ok(out);
    }
    let mut coefficients: Vec<_> = value
        .coefficients()
        .into_iter()
        .map(|x| num_bigint::BigUint::from_bytes_le(&x.into_bigint().to_bytes_le()))
        .collect();
    let radix_big = num_bigint::BigUint::from(radix);
    let mut out = Vec::with_capacity(depth);
    for _ in 0..depth {
        let digits: Vec<_> = coefficients
            .iter_mut()
            .map(|x| {
                let digit = (&*x % &radix_big)
                    .to_u32_digits()
                    .first()
                    .copied()
                    .unwrap_or(0);
                *x /= &radix_big;
                R::BaseRing::from(digit)
            })
            .collect();
        out.push(R::try_from_coefficients(&digits).ok_or(greyhound::Error::Shape)?);
    }
    if coefficients.iter().any(|x| x.bits() != 0) {
        return Err(greyhound::Error::Unsatisfied);
    }
    Ok(out)
}

#[derive(Clone, Debug, PartialEq)]
pub struct GreyhoundClaim<R> {
    pub claim: Claim<R>,
    pub context: Vec<u8>,
}
impl<R> fmt::Display for GreyhoundClaim<R> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Greyhound opening: {} evaluations",
            self.claim.evaluations.len()
        )
    }
}
/// Opening of a table under the canonical-digit Greyhound commitment.
pub struct GreyhoundOpening<R>(PhantomData<R>);
impl<R: PolyRing> Relation for GreyhoundOpening<R>
where
    R::BaseRing: PrimeField,
{
    type Size = Arc<GreyhoundKey<R>>;
    type Index = Arc<GreyhoundKey<R>>;
    type Instance = GreyhoundClaim<R>;
    type Witness = Vec<Vec<R>>;
    fn is_well_defined(i: &Self::Index, x: &Self::Instance, w: Option<&Self::Witness>) -> bool {
        Self::is_well_defined_err(i, x, w).is_ok()
    }
    fn is_well_defined_err(
        i: &Self::Index,
        x: &Self::Instance,
        w: Option<&Self::Witness>,
    ) -> Result<()> {
        if let Some(w) = w {
            ensure!(
                w.len() == i.rows && w.iter().all(|r| r.len() == i.cols),
                "table shape"
            );
        }
        greyhound::validate_claim(i, &x.claim).map_err(|e| anyhow::anyhow!("{e:?}"))?;
        Ok(())
    }
    fn is_satisfied(i: &Self::Index, x: &Self::Instance, w: &Self::Witness) -> bool {
        Self::is_satisfied_err(i, x, w).is_ok()
    }
    fn is_satisfied_err(i: &Self::Index, x: &Self::Instance, w: &Self::Witness) -> Result<()> {
        Self::is_well_defined_err(i, x, Some(w))?;
        let commitment = greyhound::commit(i, w, &|v, d| canonical_digits(v, i.radix, d))
            .map_err(|e| anyhow::anyhow!("{e:?}"))?;
        ensure!(commitment == x.claim.commitment, "commitment mismatch");
        for e in &x.claim.evaluations {
            let value: R = w
                .iter()
                .zip(&e.outer)
                .map(|(row, &a)| a * row.iter().zip(&e.inner).map(|(&b, &c)| b * c).sum::<R>())
                .sum();
            let value = match e.projection {
                Projection::Ring => value,
                Projection::ConstantCoefficient => R::from_scalar(value.coefficients()[0]),
            };
            ensure!(value == e.value, "evaluation mismatch");
        }
        Ok(())
    }
    fn generate_satisfied_instance(i: &Self::Size) -> (Self::Index, Self::Instance, Self::Witness) {
        let w = vec![vec![R::ZERO; i.cols]; i.rows];
        let commitment = greyhound::commit(i, &w, &|v, d| canonical_digits(v, i.radix, d))
            .expect("valid fixture parameters");
        let x = GreyhoundClaim {
            claim: Claim {
                commitment,
                evaluations: vec![Evaluation {
                    outer: vec![R::ONE; i.rows],
                    inner: vec![R::ONE; i.cols],
                    value: R::ZERO,
                    projection: Projection::Ring,
                }],
            },
            context: vec![],
        };
        Self::is_satisfied_err(i, &x, &w).expect("valid fixture parameters");
        (i.clone(), x, w)
    }
    fn generate_unsatisfied_instance(
        i: &Self::Size,
    ) -> (Self::Index, Self::Instance, Self::Witness) {
        let (i, mut x, w) = Self::generate_satisfied_instance(i);
        x.claim.evaluations[0].value = R::ONE;
        (i, x, w)
    }
}

/// Short affine equations and an integer squared-norm bound. The instance
/// retains the fold and transcript digest needed by subsequent arguments.
pub struct ShortAffine<R>(PhantomData<R>);
impl<R: PolyRing> Relation for ShortAffine<R> {
    type Size = ShortRelation<R>;
    type Index = ();
    type Instance = ShortRelation<R>;
    type Witness = Vec<R>;
    fn is_well_defined(i: &(), x: &Self::Instance, w: Option<&Self::Witness>) -> bool {
        Self::is_well_defined_err(i, x, w).is_ok()
    }
    fn is_well_defined_err(_: &(), x: &Self::Instance, w: Option<&Self::Witness>) -> Result<()> {
        ensure!(
            x.witness_len > 0 && x.norm_squared > 0,
            "short relation shape"
        );
        ensure!(
            w.is_none_or(|w| w.len() == x.witness_len),
            "short witness shape"
        );
        for row in &x.rows {
            ensure!(
                row.terms.iter().all(|(i, _)| *i < x.witness_len),
                "row index"
            );
            if row.projection == Projection::ConstantCoefficient {
                ensure!(
                    row.value == R::from_scalar(row.value.coefficients()[0]),
                    "constant row target"
                );
            }
        }
        Ok(())
    }
    fn is_satisfied(i: &(), x: &Self::Instance, w: &Self::Witness) -> bool {
        Self::is_satisfied_err(i, x, w).is_ok()
    }
    fn is_satisfied_err(i: &(), x: &Self::Instance, w: &Self::Witness) -> Result<()> {
        Self::is_well_defined_err(i, x, Some(w))?;
        ensure!(
            greyhound::check_disclosed_witness(x, w),
            "short equations or integer norm"
        );
        Ok(())
    }
    fn generate_satisfied_instance(size: &Self::Size) -> ((), Self::Instance, Self::Witness) {
        let mut x = size.clone();
        for row in &mut x.rows {
            row.value = R::ZERO;
        }
        let w = vec![R::ZERO; x.witness_len];
        Self::is_satisfied_err(&(), &x, &w).expect("valid fixture dimensions");
        ((), x, w)
    }
    fn generate_unsatisfied_instance(size: &Self::Size) -> ((), Self::Instance, Self::Witness) {
        let (i, mut x, w) = Self::generate_satisfied_instance(size);
        x.rows.push(Row {
            terms: vec![],
            value: R::ONE,
            projection: Projection::Ring,
        });
        (i, x, w)
    }
}
impl<R> fmt::Display for ShortRelation<R> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Short affine relation: {} words, {} rows",
            self.witness_len,
            self.rows.len()
        )
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct HachiIndex<const N: usize> {
    pub key: Arc<FoldParams<Fq>>,
    pub greyhound: Arc<GreyhoundKey<R64>>,
    pub coefficients: usize,
    pub row_log: usize,
}
#[derive(Clone, Debug, PartialEq)]
pub struct HachiClaim<const N: usize> {
    pub commitment: Vec<Pow2CyclotomicPolyRing<Fq, N>>,
    pub point: Vec<Fq4>,
    pub value: Fq4,
    pub context: Vec<u8>,
}
pub struct HachiOpening<const N: usize>;
impl<const N: usize> Relation for HachiOpening<N> {
    type Size = HachiIndex<N>;
    type Index = HachiIndex<N>;
    type Instance = HachiClaim<N>;
    type Witness = Vec<Fq>;
    fn is_well_defined(i: &Self::Index, x: &Self::Instance, w: Option<&Self::Witness>) -> bool {
        Self::is_well_defined_err(i, x, w).is_ok()
    }
    fn is_well_defined_err(
        i: &Self::Index,
        x: &Self::Instance,
        w: Option<&Self::Witness>,
    ) -> Result<()> {
        ensure!(
            N >= 8 && N.is_power_of_two() && i.key.ring_degree == N,
            "ring degree"
        );
        ensure!(
            i.coefficients.is_power_of_two() && i.coefficients >= N / 4,
            "coefficient count"
        );
        ensure!(
            x.point.len() == i.coefficients.ilog2() as usize && x.commitment.len() == i.key.b.len(),
            "opening shape"
        );
        ensure!(w.is_none_or(|w| w.len() == i.coefficients), "witness shape");
        hachi_layout(i).map_err(|_| anyhow::anyhow!("incompatible Hachi/Greyhound geometry"))?;
        Ok(())
    }
    fn is_satisfied(i: &Self::Index, x: &Self::Instance, w: &Self::Witness) -> bool {
        Self::is_satisfied_err(i, x, w).is_ok()
    }
    fn is_satisfied_err(i: &Self::Index, x: &Self::Instance, w: &Self::Witness) -> Result<()> {
        Self::is_well_defined_err(i, x, Some(w))?;
        let commitment = hachi_initial::commit::<N>(&i.key, w, i.row_log)
            .map_err(|e| anyhow::anyhow!("{e:?}"))?;
        let (value, _) =
            hachi_initial::prepare::<N>(w, &x.point).map_err(|e| anyhow::anyhow!("{e:?}"))?;
        ensure!(
            commitment == x.commitment && value == x.value,
            "Hachi commitment or evaluation"
        );
        Ok(())
    }
    fn generate_satisfied_instance(i: &Self::Size) -> (Self::Index, Self::Instance, Self::Witness) {
        let w = vec![Fq::from(0u64); i.coefficients];
        let x = HachiClaim {
            commitment: hachi_initial::commit::<N>(&i.key, &w, i.row_log)
                .expect("valid fixture parameters"),
            point: vec![Fq4::ZERO; i.coefficients.ilog2() as usize],
            value: Fq4::ZERO,
            context: vec![],
        };
        Self::is_satisfied_err(i, &x, &w).expect("valid fixture parameters");
        (i.clone(), x, w)
    }
    fn generate_unsatisfied_instance(
        i: &Self::Size,
    ) -> (Self::Index, Self::Instance, Self::Witness) {
        let (i, mut x, w) = Self::generate_satisfied_instance(i);
        x.value += Fq4::ONE;
        (i, x, w)
    }
}

// Fixed-shape field transport: no proof-controlled vector lengths are decoded.
fn put<T: CanonicalSerialize>(out: &mut Vec<u8>, v: &T) -> ProofResult<()> {
    v.serialize_compressed(out).map_err(invalid)
}
fn get<T: CanonicalDeserialize>(input: &mut &[u8]) -> ProofResult<T> {
    T::deserialize_compressed(input).map_err(invalid)
}
fn send(merlin: &mut Merlin, bytes: &[u8], expected: usize) -> ProofResult<()> {
    if bytes.len() != expected {
        return Err(ProofError::InvalidProof);
    }
    merlin.add_bytes(bytes)?;
    Ok(())
}
fn receive(arthur: &mut Arthur, size: usize) -> ProofResult<Vec<u8>> {
    let mut out = vec![0; size];
    arthur.fill_next_bytes(&mut out)?;
    Ok(out)
}
fn frame(label: &str, size: ProofResult<usize>) -> IOPattern {
    IOPattern::new(label).absorb(size.ok().filter(|&n| n > 0).unwrap_or(1), "messages")
}

fn hachi_layout<const N: usize>(i: &HachiIndex<N>) -> ProofResult<(usize, usize, usize)> {
    if N < 8
        || !N.is_power_of_two()
        || !i.coefficients.is_power_of_two()
        || i.coefficients < N / 4
        || i.row_log >= usize::BITS as usize
    {
        return Err(ProofError::InvalidProof);
    }
    let words = i.coefficients / (N / 4);
    let right = 1usize
        .checked_shl(i.row_log as u32)
        .ok_or(ProofError::InvalidProof)?;
    if right > words || i.key.ring_degree != N {
        return Err(ProofError::InvalidProof);
    }
    let (rounds, small) = hachi::opening_geometry(&i.key, words / right, right).map_err(invalid)?;
    let packed = 1usize
        .checked_shl(u32::try_from(rounds).map_err(invalid)?)
        .ok_or(ProofError::InvalidProof)?
        / 16;
    if i.greyhound.input_depth != 1
        || !i.greyhound.rows.is_power_of_two()
        || mul(i.greyhound.rows, i.greyhound.cols)? != packed
    {
        return Err(ProofError::InvalidProof);
    }
    // Initial ring, projected ring commitments, successor commitment, both
    // sumchecks (including endpoints), line polynomial, successor target.
    let fq = Fq::from(0u64).compressed_size();
    let ext = Fq4::ZERO.compressed_size();
    let ring = R64::ZERO.compressed_size();
    let scalar_words = mul(N, add(1, i.key.d.len())?)?;
    let sumcheck = add(mul(rounds, add(small, 4)?)?, add(3, rounds)?)?;
    let size = add(
        add(mul(scalar_words, fq)?, mul(i.greyhound.b.len(), ring)?)?,
        add(mul(sumcheck, ext)?, ring)?,
    )?;
    Ok((rounds, small, size))
}
fn encode_hachi<const N: usize>(i: &HachiIndex<N>, m: &HachiMessage<N>) -> ProofResult<Vec<u8>> {
    let (rounds, small, size) = hachi_layout(i)?;
    if m.fold.projected_commitment.len() != i.key.d.len()
        || m.fold.projected_commitment.iter().any(|p| p.len() != N)
        || m.fold.linear.commitment.len() != i.greyhound.b.len()
        || m.fold.linear.line.len() != rounds + 1
    {
        return Err(ProofError::InvalidProof);
    }
    let mut out = Vec::with_capacity(size);
    put(&mut out, &m.initial.ring_value)?;
    for p in &m.fold.projected_commitment {
        for x in p {
            put(&mut out, x)?;
        }
    }
    for c in &m.fold.linear.commitment {
        put(&mut out, c)?;
    }
    for (s, degree) in [
        (&m.fold.linear.smallness, small),
        (&m.fold.linear.relation, 2),
    ] {
        if s.rounds.len() != rounds || s.rounds.iter().any(|r| r.len() != degree + 1) {
            return Err(ProofError::InvalidProof);
        }
        for r in &s.rounds {
            for x in r {
                put(&mut out, x)?;
            }
        }
        put(&mut out, &s.endpoint)?;
    }
    for x in &m.fold.linear.line {
        put(&mut out, x)?;
    }
    put(&mut out, &m.adapter.ring_target)?;
    if out.len() != size {
        return Err(ProofError::InvalidProof);
    }
    Ok(out)
}
fn decode_hachi<const N: usize>(i: &HachiIndex<N>, bytes: &[u8]) -> ProofResult<HachiMessage<N>> {
    let (rounds, small, size) = hachi_layout(i)?;
    if bytes.len() != size {
        return Err(ProofError::InvalidProof);
    }
    let mut input = bytes;
    let initial = InitialMessage {
        ring_value: get(&mut input)?,
    };
    let projected_commitment = (0..i.key.d.len())
        .map(|_| (0..N).map(|_| get(&mut input)).collect())
        .collect::<ProofResult<_>>()?;
    let commitment = (0..i.greyhound.b.len())
        .map(|_| get(&mut input))
        .collect::<ProofResult<_>>()?;
    let mut sumchecks = Vec::new();
    for degree in [small, 2] {
        let messages = (0..rounds)
            .map(|_| (0..=degree).map(|_| get(&mut input)).collect())
            .collect::<ProofResult<_>>()?;
        sumchecks.push(Sumcheck {
            rounds: messages,
            endpoint: get(&mut input)?,
        });
    }
    let line = (0..=rounds)
        .map(|_| get(&mut input))
        .collect::<ProofResult<_>>()?;
    let adapter = AdapterProof {
        ring_target: get(&mut input)?,
    };
    if !input.is_empty() {
        return Err(ProofError::InvalidProof);
    }
    let relation = sumchecks.pop().ok_or(ProofError::InvalidProof)?;
    let smallness = sumchecks.pop().ok_or(ProofError::InvalidProof)?;
    Ok(HachiMessage {
        initial,
        fold: FoldProof {
            projected_commitment,
            linear: LinearProof {
                commitment,
                smallness,
                relation,
                line,
            },
        },
        adapter,
    })
}

pub struct HachiToGreyhound<const N: usize>;
impl<const N: usize> HachiToGreyhound<N> {
    pub fn message_size(i: &HachiIndex<N>) -> ProofResult<usize> {
        hachi_layout(i).map(|x| x.2)
    }
}
impl<const N: usize> Reduction<HachiOpening<N>, GreyhoundOpening<R64>, ()> for HachiToGreyhound<N> {
    fn iopattern(_: &(), i: &HachiIndex<N>, _: &HachiClaim<N>) -> IOPattern {
        frame(
            "lattirust/hachi/reduction-transport/v1",
            hachi_layout(i).map(|x| x.2),
        )
    }
    fn prove(
        _: &(),
        i: &HachiIndex<N>,
        x: &HachiClaim<N>,
        w: &Vec<Fq>,
        merlin: &mut Merlin,
    ) -> ProofResult<(Arc<GreyhoundKey<R64>>, GreyhoundClaim<R64>, Vec<Vec<R64>>)> {
        HachiOpening::is_well_defined_err(i, x, Some(w)).map_err(invalid)?;
        let (value, message, claim, table, context) = hachi_initial::prepare_hachi_reference(
            &i.key,
            &i.greyhound,
            &x.commitment,
            i.row_log,
            &x.point,
            w,
            &x.context,
            &|v, d| canonical_digits(v, i.greyhound.radix, d),
        )
        .map_err(invalid)?;
        if value != x.value {
            return Err(ProofError::InvalidProof);
        }
        let bytes = encode_hachi(i, &message)?;
        send(merlin, &bytes, hachi_layout(i)?.2)?;
        Ok((
            i.greyhound.clone(),
            GreyhoundClaim {
                claim,
                context: context.to_vec(),
            },
            table,
        ))
    }
    fn verify(
        _: &(),
        i: &HachiIndex<N>,
        x: &HachiClaim<N>,
        arthur: &mut Arthur,
    ) -> ProofResult<(Arc<GreyhoundKey<R64>>, GreyhoundClaim<R64>)> {
        HachiOpening::is_well_defined_err(i, x, None).map_err(invalid)?;
        let bytes = receive(arthur, hachi_layout(i)?.2)?;
        let message = decode_hachi(i, &bytes)?;
        let (claim, context) = hachi_initial::reduce_hachi(
            &i.key,
            &i.greyhound,
            &x.commitment,
            i.row_log,
            &x.point,
            x.value,
            &x.context,
            &message,
        )
        .map_err(invalid)?;
        Ok((
            i.greyhound.clone(),
            GreyhoundClaim {
                claim,
                context: context.to_vec(),
            },
        ))
    }
}

pub struct GreyhoundToShort<R>(PhantomData<R>);
impl<R: PolyRing> GreyhoundToShort<R> {
    pub fn message_size(i: &GreyhoundKey<R>, x: &GreyhoundClaim<R>) -> ProofResult<usize> {
        greyhound_size(i, x)
    }
}
fn greyhound_size<R: PolyRing>(i: &GreyhoundKey<R>, x: &GreyhoundClaim<R>) -> ProofResult<usize> {
    mul(
        mul(i.d.len(), x.claim.evaluations.len())?,
        R::ZERO.compressed_size(),
    )
}
impl<R: PolyRing> Reduction<GreyhoundOpening<R>, ShortAffine<R>, ()> for GreyhoundToShort<R>
where
    R::BaseRing: PrimeField,
{
    fn iopattern(_: &(), i: &Arc<GreyhoundKey<R>>, x: &GreyhoundClaim<R>) -> IOPattern {
        frame(
            "lattirust/greyhound/reduction-transport/v1",
            greyhound_size(i, x),
        )
    }
    fn prove(
        _: &(),
        i: &Arc<GreyhoundKey<R>>,
        x: &GreyhoundClaim<R>,
        w: &Vec<Vec<R>>,
        merlin: &mut Merlin,
    ) -> ProofResult<((), ShortRelation<R>, Vec<R>)> {
        GreyhoundOpening::is_well_defined_err(i, x, Some(w)).map_err(invalid)?;
        let prepared = greyhound::prepare_opening(
            i,
            &x.claim,
            w,
            &|v, d| canonical_digits(v, i.radix, d),
            &x.context,
        )
        .map_err(invalid)?;
        let mut bytes = Vec::new();
        for row in &prepared.message.auxiliary_commitments {
            for v in row {
                put(&mut bytes, v)?;
            }
        }
        send(merlin, &bytes, greyhound_size(i, x)?)?;
        Ok(((), prepared.relation, prepared.short_witness))
    }
    fn verify(
        _: &(),
        i: &Arc<GreyhoundKey<R>>,
        x: &GreyhoundClaim<R>,
        arthur: &mut Arthur,
    ) -> ProofResult<((), ShortRelation<R>)> {
        GreyhoundOpening::is_well_defined_err(i, x, None).map_err(invalid)?;
        let bytes = receive(arthur, greyhound_size(i, x)?)?;
        let mut input = bytes.as_slice();
        let auxiliary_commitments = (0..x.claim.evaluations.len())
            .map(|_| (0..i.d.len()).map(|_| get(&mut input)).collect())
            .collect::<ProofResult<_>>()?;
        let relation = greyhound::reduce(
            i,
            &x.claim,
            &OpeningMessage {
                auxiliary_commitments,
            },
            &x.context,
        )
        .map_err(invalid)?;
        Ok(((), relation))
    }
}

/// Deterministic, zero-message conversion. Padding remains in the norm bound.
pub struct ShortToPrincipal<R>(PhantomData<R>);
impl<R: PolyRing> Reduction<ShortAffine<R>, PrincipalRelation<R>, usize> for ShortToPrincipal<R>
where
    R::BaseRing: lattirust_arithmetic::ring::representatives::WithSignedRepresentative,
{
    fn iopattern(_: &usize, _: &(), _: &ShortRelation<R>) -> IOPattern {
        IOPattern::new("lattirust/short-affine/principal/v1")
    }
    fn prove(
        width: &usize,
        i: &(),
        x: &ShortRelation<R>,
        w: &Vec<R>,
        _: &mut Merlin,
    ) -> ProofResult<(Index<R>, Instance<R>, Witness<R>)> {
        ShortAffine::is_well_defined_err(i, x, Some(w)).map_err(invalid)?;
        let (index, instance) = labrador_relation::convert(x, *width).map_err(invalid)?;
        let rows = w
            .chunks(index.n)
            .map(|chunk| {
                let mut row = Vec::with_capacity(index.n);
                row.extend_from_slice(chunk);
                row.resize(index.n, R::ZERO);
                Vector::from(row)
            })
            .collect();
        let witness = Witness::new(rows);
        Ok((index, instance, witness))
    }
    fn verify(
        width: &usize,
        i: &(),
        x: &ShortRelation<R>,
        _: &mut Arthur,
    ) -> ProofResult<(Index<R>, Instance<R>)> {
        ShortAffine::is_well_defined_err(i, x, None).map_err(invalid)?;
        labrador_relation::convert(x, *width).map_err(invalid)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Fq2;
    fn fixture() -> HachiIndex<8> {
        let s = |v: u64| R64::from(Fq::from(v));
        HachiIndex {
            key: Arc::new(FoldParams {
                ring_degree: 8,
                base: 16,
                folded_digits: 2,
                challenge_weight: 1,
                a: vec![(1..=32).map(|i| vec![Fq::from(i as u64)]).collect()],
                b: vec![(9..=16).map(|i| vec![Fq::from(i as u64)]).collect()],
                d: vec![(17..=24).map(|i| vec![Fq::from(i as u64)]).collect()],
            }),
            greyhound: Arc::new(GreyhoundKey {
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
            }),
            coefficients: 8,
            row_log: 0,
        }
    }
    #[test]
    fn existing_framework_completeness_checks_all_three_handoffs() {
        let i = fixture();
        HachiToGreyhound::<8>::test_completeness(&(), &i).unwrap();
        GreyhoundToShort::<R64>::test_completeness(&(), &i.greyhound).unwrap();
        let (_, x, w) = GreyhoundOpening::generate_satisfied_instance(&i.greyhound);
        let io = GreyhoundToShort::iopattern(&(), &i.greyhound, &x);
        let (_, short, _) =
            GreyhoundToShort::prove(&(), &i.greyhound, &x, &w, &mut io.to_merlin()).unwrap();
        ShortToPrincipal::<R64>::test_completeness(&32, &short).unwrap();
    }
    #[test]
    fn nonzero_handoffs_agree_and_preserve_integer_norm_and_padding() {
        let i = fixture();
        let w: Vec<_> = (0..8).map(|i| Fq::from(i as u64)).collect();
        let a = Fq4::new(
            Fq2::new(4u64.into(), 2u64.into()),
            Fq2::new(8u64.into(), 3u64.into()),
        );
        let point = vec![a, a + Fq4::ONE, a];
        let (value, _) = hachi_initial::prepare::<8>(&w, &point).unwrap();
        let x = HachiClaim {
            commitment: hachi_initial::commit::<8>(&i.key, &w, 0).unwrap(),
            point,
            value,
            context: b"reduction-framework/nonzero".to_vec(),
        };
        assert!(HachiOpening::is_satisfied(&i, &x, &w));
        let io = HachiToGreyhound::iopattern(&(), &i, &x);
        let mut merlin = io.to_merlin();
        let (key, query, table) = HachiToGreyhound::prove(&(), &i, &x, &w, &mut merlin).unwrap();
        let (vkey, vquery) =
            HachiToGreyhound::verify(&(), &i, &x, &mut io.to_arthur(merlin.transcript())).unwrap();
        assert_eq!(key, vkey);
        assert_eq!(query, vquery);
        assert!(GreyhoundOpening::is_satisfied(&key, &query, &table));
        let gio = GreyhoundToShort::iopattern(&(), &key, &query);
        let mut gm = gio.to_merlin();
        let (_, short, sw) = GreyhoundToShort::prove(&(), &key, &query, &table, &mut gm).unwrap();
        let (_, vs) =
            GreyhoundToShort::verify(&(), &key, &query, &mut gio.to_arthur(gm.transcript()))
                .unwrap();
        assert_eq!(short, vs);
        assert!(ShortAffine::is_satisfied(&(), &short, &sw));
        let io = ShortToPrincipal::iopattern(&32, &(), &short);
        let mut empty = io.to_merlin();
        let (pi, px, pw) = ShortToPrincipal::prove(&32, &(), &short, &sw, &mut empty).unwrap();
        assert!(empty.transcript().is_empty());
        assert_eq!(
            (pi.clone(), px.clone()),
            ShortToPrincipal::verify(&32, &(), &short, &mut io.to_arthur(&[])).unwrap()
        );
        assert_eq!(pi.norm_bound_squared, short.norm_squared as f64);
        assert_eq!(pi.r * pi.n, 96);
        assert_eq!(sw.len(), 80);
        assert!(PrincipalRelation::is_satisfied(&pi, &px, &pw));
        let mut wrong = sw;
        wrong[0] += R64::from(Fq::from(100_000u64));
        assert!(!ShortAffine::is_satisfied(&(), &short, &wrong));
        let (_, _, bad) =
            ShortToPrincipal::prove(&32, &(), &short, &wrong, &mut io.to_merlin()).unwrap();
        assert!(!PrincipalRelation::is_satisfied(&pi, &px, &bad));
        let mut wrong = x.clone();
        wrong.value += Fq4::ONE;
        assert!(HachiToGreyhound::verify(
            &(),
            &i,
            &wrong,
            &mut io_for_hachi(&i, &wrong).to_arthur(merlin.transcript())
        )
        .is_err());
        let mut bad_index = i.clone();
        Arc::make_mut(&mut bad_index.greyhound).cols = 15;
        assert!(HachiToGreyhound::prove(
            &(),
            &bad_index,
            &x,
            &w,
            &mut HachiToGreyhound::iopattern(&(), &bad_index, &x).to_merlin()
        )
        .is_err());
    }
    fn io_for_hachi(i: &HachiIndex<8>, x: &HachiClaim<8>) -> IOPattern {
        HachiToGreyhound::iopattern(&(), i, x)
    }
    #[test]
    fn canonical_digit_helper_rejects_short_decompositions() {
        let value = R64::from(Fq::from(256u64));
        assert!(canonical_digits(value, 16, 2).is_err());
        let digits = canonical_digits(value, 16, 8).unwrap();
        assert_eq!(digits[2], R64::ONE);
        assert_eq!(digits[0], R64::ZERO);
        let decimal = canonical_digits(value, 10, 3).unwrap();
        assert_eq!(
            decimal,
            vec![
                R64::from(Fq::from(6u64)),
                R64::from(Fq::from(5u64)),
                R64::from(Fq::from(2u64))
            ]
        );
    }
}
