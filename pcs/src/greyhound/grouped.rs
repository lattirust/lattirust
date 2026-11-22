//! Greyhound Fig. 2: each group shares its evaluation weights, but consists
//! of independently folded tables. One aggregate commitment u and one v.
use super::*;

#[derive(Clone, Debug)]
pub struct GroupKey<R> {
    /// A column block per table, of length rows*rank(A)*depth.
    pub b: Vec<Vec<R>>,
    /// A column block per table, of length rows*depth.
    pub d: Vec<Vec<R>>,
}

#[derive(Clone, Debug)]
pub struct GroupedKey<R> {
    pub rows: usize,
    pub cols: usize,
    pub radix: u32,
    pub depth: usize,
    pub input_depth: usize,
    pub challenge_weight: usize,
    pub norm_squared: u64,
    pub a: Vec<Vec<R>>,
    pub groups: Vec<GroupKey<R>>,
}

#[derive(Clone, Debug)]
pub struct Group<R> {
    pub inner: Vec<R>,
    pub outer: Vec<R>,
    /// One claimed evaluation per table in this group.
    pub values: Vec<R>,
}

#[derive(Clone, Debug)]
pub struct GroupedClaim<R> {
    pub commitment: Vec<R>,
    pub groups: Vec<Group<R>>,
}

/// Actual Figure 2 aggregate relation. The order of the returned witness is
/// [all z by group, all t_hat by group/table, all w_hat by group/table].
pub fn reduce_grouped<R: PolyRing>(
    key: &GroupedKey<R>,
    expected: &GroupedClaim<R>,
    auxiliary_commitment: &[R],
    context: &[u8],
) -> Result<ShortRelation<R>, Error> {
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
        || k.groups.is_empty()
    {
        return Err(Error::Parameters);
    }
    if k.groups.len() != expected.groups.len()
        || expected.commitment.is_empty()
        || auxiliary_commitment.is_empty()
    {
        return Err(Error::Shape);
    }
    let z_len = checked_mul(k.cols, k.input_depth)?;
    let t_block = checked_mul(k.a.len(), k.depth)?;
    let t_len = checked_mul(k.rows, t_block)?;
    let w_len = checked_mul(k.rows, k.depth)?;
    let t_start = checked_mul(k.groups.len(), z_len)?;
    let tables = expected.groups.iter().try_fold(0usize, |acc, g| {
        acc.checked_add(g.values.len()).ok_or(Error::Parameters)
    })?;
    let w_start = t_start
        .checked_add(checked_mul(tables, t_len)?)
        .ok_or(Error::Parameters)?;
    let witness_len = w_start
        .checked_add(checked_mul(tables, w_len)?)
        .ok_or(Error::Parameters)?;
    if k.a.iter().any(|row| row.len() != z_len) {
        return Err(Error::Shape);
    }
    let mut transcript = Shake256::default();
    bytes(
        &mut transcript,
        b"lattirust/greyhound/grouped-short-relation/v2",
    )?;
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
    values(&mut transcript, &expected.commitment)?;
    values(&mut transcript, auxiliary_commitment)?;
    number(&mut transcript, k.groups.len())?;
    for (keys, claims) in k.groups.iter().zip(&expected.groups) {
        let t_group = checked_mul(claims.values.len(), t_len)?;
        let w_group = checked_mul(claims.values.len(), w_len)?;
        if claims.values.is_empty()
            || claims.inner.len() != k.cols
            || claims.outer.len() != k.rows
            || keys.b.len() != expected.commitment.len()
            || keys.d.len() != auxiliary_commitment.len()
            || keys.b.iter().any(|row| row.len() != t_group)
            || keys.d.iter().any(|row| row.len() != w_group)
        {
            return Err(Error::Shape);
        }
        matrix(&mut transcript, &keys.b)?;
        matrix(&mut transcript, &keys.d)?;
        values(&mut transcript, &claims.inner)?;
        values(&mut transcript, &claims.outer)?;
        values(&mut transcript, &claims.values)?;
    }
    let mut reader = transcript.finalize_xof();
    let mut transcript_digest = [0u8; 32];
    reader.read(&mut transcript_digest);
    let fold = (0..checked_mul(tables, k.rows)?)
        .map(|_| sparse::<R>(&mut reader, k.challenge_weight))
        .collect::<Result<Vec<_>, _>>()?;
    let g = radix_powers::<R>(k.radix, k.depth)?;
    let g_in = radix_powers::<R>(k.radix, k.input_depth)?;
    let make_row = |terms, value| Row {
        terms,
        value,
        projection: Projection::Ring,
    };
    let mut offset = 0;
    let mut b_rows: Vec<_> = expected
        .commitment
        .iter()
        .map(|&x| make_row(Vec::new(), x))
        .collect();
    let mut d_rows: Vec<_> = auxiliary_commitment
        .iter()
        .map(|&x| make_row(Vec::new(), x))
        .collect();
    let mut rows = Vec::new();
    for (group_idx, (keys, claims)) in k.groups.iter().zip(&expected.groups).enumerate() {
        for (row, b) in b_rows.iter_mut().zip(&keys.b) {
            row.terms.extend(
                b.iter()
                    .enumerate()
                    .map(|(i, &x)| (t_start + offset * t_len + i, x)),
            );
        }
        for (row, d) in d_rows.iter_mut().zip(&keys.d) {
            row.terms.extend(
                d.iter()
                    .enumerate()
                    .map(|(i, &x)| (w_start + offset * w_len + i, x)),
            );
        }
        let mut a_rows: Vec<_> =
            k.a.iter()
                .map(|a| {
                    make_row(
                        a.iter()
                            .enumerate()
                            .map(|(i, &x)| (group_idx * z_len + i, x))
                            .collect(),
                        R::ZERO,
                    )
                })
                .collect();
        let mut folded_evaluation = make_row(Vec::new(), R::ZERO);
        for (i, &x) in claims.inner.iter().enumerate() {
            for (j, &scale) in g_in.iter().enumerate() {
                folded_evaluation
                    .terms
                    .push((group_idx * z_len + i * k.input_depth + j, x * scale));
            }
        }
        for (table_idx, &value) in claims.values.iter().enumerate() {
            let table_idx = offset + table_idx;
            let mut evaluation = make_row(Vec::new(), value);
            for (i, &outer) in claims.outer.iter().enumerate() {
                let c = fold[table_idx * k.rows + i];
                for (j, &scale) in g.iter().enumerate() {
                    let w_index = w_start + table_idx * w_len + i * k.depth + j;
                    evaluation.terms.push((w_index, outer * scale));
                    folded_evaluation.terms.push((w_index, -(c * scale)));
                    for (a_idx, row) in a_rows.iter_mut().enumerate() {
                        let t_index =
                            t_start + table_idx * t_len + i * t_block + a_idx * k.depth + j;
                        row.terms.push((t_index, -(c * scale)));
                    }
                }
            }
            rows.push(evaluation);
        }
        rows.push(folded_evaluation);
        rows.extend(a_rows);
        offset += claims.values.len();
    }
    rows.extend(b_rows);
    rows.extend(d_rows);
    Ok(ShortRelation {
        witness_len,
        rows,
        norm_squared: k.norm_squared,
        fold,
        transcript_digest,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use lattirust_arithmetic::ring::{f_p::Fq, Pow2CyclotomicPolyRing};

    type F = Fq<4294967197>;
    type R = Pow2CyclotomicPolyRing<F, 8>;
    fn r(n: u32) -> R {
        R::from_scalar(F::from(n))
    }

    #[test]
    fn figure_two_binds_global_commitments_and_separate_group_folds() {
        let keys = GroupedKey {
            rows: 1,
            cols: 1,
            radix: 16,
            depth: 1,
            input_depth: 1,
            challenge_weight: 2,
            norm_squared: 1000,
            a: vec![vec![r(2)]],
            groups: vec![
                GroupKey {
                    b: vec![vec![r(2), r(3)]],
                    d: vec![vec![r(1), r(2)]],
                },
                GroupKey {
                    b: vec![vec![r(4)]],
                    d: vec![vec![r(3)]],
                },
            ],
        };
        // Three independent 1×1 tables, entries 2,3,4. t_hat=(4,6,8).
        // First group inner 2, outer 3: w_hat=(4,6); second inner 1, outer 5: w_hat=(4).
        let mut claim = GroupedClaim {
            commitment: vec![r(58)],
            groups: vec![
                Group {
                    inner: vec![r(2)],
                    outer: vec![r(3)],
                    values: vec![r(12), r(18)],
                },
                Group {
                    inner: vec![r(1)],
                    outer: vec![r(5)],
                    values: vec![r(20)],
                },
            ],
        };
        let v = [r(28)];
        let rel = reduce_grouped(&keys, &claim, &v, b"grouped fixture").unwrap();
        assert_eq!(rel.rows.len(), 9);
        assert_eq!(rel.fold.len(), 3);
        let mut witness = vec![rel.fold[0] * r(2) + rel.fold[1] * r(3), rel.fold[2] * r(4)];
        witness.extend([r(4), r(6), r(8), r(4), r(6), r(4)]);
        assert!(check_disclosed_witness(&rel, &witness));
        claim.groups[0].values[1] += r(1);
        let changed = reduce_grouped(&keys, &claim, &v, b"grouped fixture").unwrap();
        assert_ne!(rel.transcript_digest, changed.transcript_digest);
        assert!(!check_disclosed_witness(&changed, &witness));
        assert!(matches!(
            reduce_grouped(&keys, &claim, &[], b"grouped fixture"),
            Err(Error::Shape)
        ));
    }
}
