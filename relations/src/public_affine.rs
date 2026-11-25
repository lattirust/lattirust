//! Public coefficient generation for affine checks on a bit-packed tape.
//!
//! These helpers do not authenticate the tape, prove bitness or ranges, or
//! sample a batching challenge. The caller must provide those guarantees.
//! The output can feed a sparse constant-coefficient relation builder.
use std::fmt;

const BIT_LIMIT: u64 = 1 << 32;
const DIGITS: usize = 4;
const MAX_ACCUMULATED_WEIGHT: u64 = 1 + (DIGITS as u64) * 255;

/// Malformed public geometry or an invalid arithmetic bound.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AffineError {
    Modulus,
    Coefficient,
    Index,
    Layout,
    Overlap,
    UnsupportedBit,
    Overflow,
    WrapPossible,
}

impl fmt::Display for AffineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "public affine coefficients: {self:?}")
    }
}

impl std::error::Error for AffineError {}

/// A canonical coefficient on an absolute bit position, not a ring index.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BitCoeff {
    pub bit: u32,
    pub coeff: u32,
}

/// A little-endian binary limb. Each lane weights this limb by a scalar
/// in the source coefficient ring. A raw tape bit is a run of width one.
#[derive(Clone, Copy, Debug)]
pub struct BitRun<const LANES: usize> {
    pub start: u32,
    pub width: u8,
    pub coeff: [u32; LANES],
}

fn validate_p(p: u32) -> Result<(), AffineError> {
    if p < 3 || p % 2 == 0 {
        Err(AffineError::Modulus)
    } else {
        Ok(())
    }
}

fn validate_coeff<const LANES: usize>(p: u32, values: &[u32; LANES]) -> Result<(), AffineError> {
    if values.iter().any(|&x| x >= p) {
        Err(AffineError::Coefficient)
    } else {
        Ok(())
    }
}

fn expand(start: u32, width: u8, p: u32, c: u32) -> impl Iterator<Item = BitCoeff> {
    let mut coeff = u64::from(c);
    let modulus = u64::from(p);
    (0..u32::from(width)).map(move |j| {
        let out = BitCoeff {
            bit: start + j,
            coeff: coeff as u32,
        };
        coeff *= 2;
        if coeff >= modulus {
            coeff -= modulus;
        }
        out
    })
}

/// Validate one disjoint typed plan once, then emit any lane without a sort.
///
/// The caller may sort source runs by `start` first. If runs alias the same
/// bit, use a separate modular deduplication path; this type rejects overlap
/// rather than silently replacing one of the terms.
pub struct SortedBitRuns<'a, const LANES: usize> {
    p: u32,
    runs: &'a [BitRun<LANES>],
}

impl<'a, const LANES: usize> SortedBitRuns<'a, LANES> {
    pub fn new(p: u32, runs: &'a [BitRun<LANES>]) -> Result<Self, AffineError> {
        validate_p(p)?;
        let mut last_end = 0;
        for run in runs {
            let start = u64::from(run.start);
            let end = start + u64::from(run.width);
            if run.width == 0 || run.width > 32 || end > BIT_LIMIT {
                return Err(AffineError::Layout);
            }
            if start < last_end {
                return Err(AffineError::Overlap);
            }
            validate_coeff(p, &run.coeff)?;
            last_end = end;
        }
        Ok(Self { p, runs })
    }

    /// Ascending, nonzero terms, emitted without intermediate allocation.
    pub fn lane(&self, lane: usize) -> Result<impl Iterator<Item = BitCoeff> + 'a, AffineError> {
        if lane >= LANES {
            return Err(AffineError::Index);
        }
        let p = self.p;
        Ok(self
            .runs
            .iter()
            .filter(move |run| run.coeff[lane] != 0)
            .flat_map(move |run| expand(run.start, run.width, p, run.coeff[lane])))
    }
}

/// Geometry of repeated little-endian digit groups. For group g and block
/// j the first bit is `start + g*group_stride + j*block_stride`.
#[derive(Clone, Copy, Debug)]
pub struct StridedLayout {
    pub start: u32,
    pub group_stride: u32,
    pub block_stride: u32,
    pub blocks: u16,
    pub width: u8,
}

/// One coefficient vector per group, reused on that group's digit blocks.
/// This covers regular layouts without materializing individual bit runs.
pub struct StridedBitRuns<'a, const LANES: usize> {
    p: u32,
    layout: StridedLayout,
    weights: &'a [[u32; LANES]],
}

impl<'a, const LANES: usize> StridedBitRuns<'a, LANES> {
    pub fn new(
        p: u32,
        layout: StridedLayout,
        weights: &'a [[u32; LANES]],
    ) -> Result<Self, AffineError> {
        validate_p(p)?;
        let width = u64::from(layout.width);
        let block_stride = u64::from(layout.block_stride);
        let group_stride = u64::from(layout.group_stride);
        if width == 0 || width > 32 || layout.blocks == 0 || block_stride < width {
            return Err(AffineError::Layout);
        }
        let span = u64::from(layout.blocks - 1) * block_stride + width;
        if span > group_stride {
            return Err(AffineError::Overlap);
        }
        if !weights.is_empty() {
            let end = u64::try_from(weights.len() - 1)
                .map_err(|_| AffineError::Overflow)?
                .checked_mul(group_stride)
                .and_then(|x| x.checked_add(u64::from(layout.start)))
                .and_then(|x| x.checked_add(span))
                .ok_or(AffineError::Overflow)?;
            if end > BIT_LIMIT {
                return Err(AffineError::Layout);
            }
        }
        for v in weights {
            validate_coeff(p, v)?;
        }
        Ok(Self { p, layout, weights })
    }

    pub fn lane(&self, lane: usize) -> Result<impl Iterator<Item = BitCoeff> + 'a, AffineError> {
        if lane >= LANES {
            return Err(AffineError::Index);
        }
        let p = self.p;
        let l = self.layout;
        Ok(self
            .weights
            .iter()
            .enumerate()
            .filter(move |(_, w)| w[lane] != 0)
            .flat_map(move |(i, w)| {
                let group = u64::from(l.start) + i as u64 * u64::from(l.group_stride);
                (0..u64::from(l.blocks)).flat_map(move |j| {
                    let start = group + j * u64::from(l.block_stride);
                    expand(start as u32, l.width, p, w[lane])
                })
            }))
    }
}

/// Rank on the union of public bit support. Word i covers bits 64i..64i+64
/// from its least significant end. Exactly one 32-bit prefix per word is
/// allocated; holes in the support do not allocate accumulator columns.
pub struct SupportRank<'a> {
    words: &'a [u64],
    prefix: Vec<u32>,
    len: usize,
}

impl<'a> SupportRank<'a> {
    pub fn new(words: &'a [u64]) -> Result<Self, AffineError> {
        if u64::try_from(words.len()).map_or(true, |n| n > BIT_LIMIT / 64) {
            return Err(AffineError::Layout);
        }
        let mut prefix = Vec::with_capacity(words.len());
        let mut count = 0u32;
        for &word in words {
            prefix.push(count);
            count = count
                .checked_add(word.count_ones())
                .ok_or(AffineError::Overflow)?;
        }
        Ok(Self {
            words,
            prefix,
            len: usize::try_from(count).map_err(|_| AffineError::Overflow)?,
        })
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Zero-based address among present bits, or None for a hole.
    pub fn position(&self, bit: u32) -> Option<usize> {
        let i = (bit >> 6) as usize;
        let w = *self.words.get(i)?;
        let mask = 1u64 << (bit & 63);
        if w & mask == 0 {
            return None;
        }
        Some((self.prefix[i] + (w & (mask - 1)).count_ones()) as usize)
    }
}

/// A private integer dot sum is typically committed as this contiguous
/// unsigned binary span (signed/biased encodings use a separate constant).
#[derive(Clone, Copy, Debug)]
pub struct SumBits {
    pub start: u32,
    pub width: u32,
}

/// Four base-256 planes, and LANES outer scalar rows. Coefficients are
/// supplied as canonical native p32 representatives. This implements a
/// public linear map, not the proof or soundness of batching the planes.
pub struct Radix256<const LANES: usize> {
    p: u32,
    q: u64,
    weights: [[u64; LANES]; DIGITS],
}

impl<const LANES: usize> Radix256<LANES> {
    pub fn new(p: u32, q: u64, weights: [[u64; LANES]; DIGITS]) -> Result<Self, AffineError> {
        validate_p(p)?;
        if q < 2 || q > u64::MAX / MAX_ACCUMULATED_WEIGHT {
            return Err(AffineError::Modulus);
        }
        if weights.iter().flatten().any(|&a| a >= q) {
            return Err(AffineError::Coefficient);
        }
        Ok(Self { p, q, weights })
    }

    /// Build the four digit-plane contributions in a single pass.
    /// `out` is support-major; each entry holds all ambient lanes. A call
    /// returning Err may already have updated a prefix of the accumulator.
    pub fn add_terms(
        &self,
        support: &SupportRank<'_>,
        terms: impl IntoIterator<Item = BitCoeff>,
        out: &mut [[u64; LANES]],
    ) -> Result<(), AffineError> {
        if out.len() != support.len() {
            return Err(AffineError::Layout);
        }
        for row in terms {
            let i = support
                .position(row.bit)
                .ok_or(AffineError::UnsupportedBit)?;
            if row.coeff >= self.p {
                return Err(AffineError::Coefficient);
            }
            let d = row.coeff.to_le_bytes();
            for (lane, c) in out[i].iter_mut().enumerate() {
                if *c >= self.q {
                    return Err(AffineError::Coefficient);
                }
                let z = *c
                    + u64::from(d[0]) * self.weights[0][lane]
                    + u64::from(d[1]) * self.weights[1][lane]
                    + u64::from(d[2]) * self.weights[2][lane]
                    + u64::from(d[3]) * self.weights[3][lane];
                *c = z % self.q;
            }
        }
        Ok(())
    }

    /// Subtract each committed dot-sum expansion from its own digit
    /// plane. If the encoding is biased, the caller adds its known
    /// constant to the corresponding public equation independently.
    /// A call returning Err may have updated a prefix of the accumulator.
    pub fn subtract_sums(
        &self,
        support: &SupportRank<'_>,
        spans: [SumBits; DIGITS],
        out: &mut [[u64; LANES]],
    ) -> Result<(), AffineError> {
        if out.len() != support.len() {
            return Err(AffineError::Layout);
        }
        for (j, span) in spans.into_iter().enumerate() {
            if span.width == 0 || u64::from(span.start) + u64::from(span.width) > BIT_LIMIT {
                return Err(AffineError::Layout);
            }
            let mut z = self.weights[j].map(|c| if c == 0 { 0 } else { self.q - c });
            for b in 0..span.width {
                let i = support
                    .position(span.start + b)
                    .ok_or(AffineError::UnsupportedBit)?;
                for (c, d) in out[i].iter_mut().zip(z.iter_mut()) {
                    if *c >= self.q {
                        return Err(AffineError::Coefficient);
                    }
                    *c += *d;
                    if *c >= self.q {
                        *c -= self.q;
                    }
                    *d *= 2;
                    if *d >= self.q {
                        *d -= self.q;
                    }
                }
            }
        }
        Ok(())
    }
}

/// A contribution `coeff * x` to the integer lift of a carry equation.
/// The caller must prove that x lies in [min, max].
#[derive(Clone, Copy, Debug)]
pub struct BoundedTerm {
    pub coeff: i128,
    pub min: i128,
    pub max: i128,
}

/// Check that modular equality to zero forces integer equality to zero.
///
/// Given the stated per-variable intervals, this returns the conservative
/// residual interval only when it lies strictly inside (-q, q). It does
/// not establish the intervals, canonical encoding, or carry recurrence.
pub fn no_wrap_residual(
    q: u64,
    constant: i128,
    terms: impl IntoIterator<Item = BoundedTerm>,
) -> Result<(i128, i128), AffineError> {
    if q < 2 {
        return Err(AffineError::Modulus);
    }
    let (mut lo, mut hi) = (constant, constant);
    for t in terms {
        if t.min > t.max {
            return Err(AffineError::Layout);
        }
        let a = t.coeff.checked_mul(t.min).ok_or(AffineError::Overflow)?;
        let b = t.coeff.checked_mul(t.max).ok_or(AffineError::Overflow)?;
        lo = lo.checked_add(a.min(b)).ok_or(AffineError::Overflow)?;
        hi = hi.checked_add(a.max(b)).ok_or(AffineError::Overflow)?;
    }
    let q = i128::from(q);
    if lo <= -q || hi >= q {
        return Err(AffineError::WrapPossible);
    }
    Ok((lo, hi))
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: u32 = u32::MAX - 98;
    const Q: u64 = (1 << 40) - 195;

    fn slow_bits(start: u32, width: u8, c: u32) -> Vec<BitCoeff> {
        (0..u32::from(width))
            .filter(|_| c != 0)
            .map(|j| BitCoeff {
                bit: start + j,
                coeff: (u64::from(c) * (1u64 << j) % u64::from(P)) as u32,
            })
            .collect()
    }

    #[test]
    fn typed_runs_match_reference_at_high_addresses() {
        let input = [
            BitRun {
                start: 9,
                width: 1,
                coeff: [0, 42],
            },
            BitRun {
                start: 147,
                width: 32,
                coeff: [P - 1, P - 42],
            },
            BitRun {
                start: u32::MAX - 3,
                width: 4,
                coeff: [P - 2, 0],
            },
        ];
        let view = SortedBitRuns::new(P, &input).unwrap();
        for lane in 0..2 {
            let expect: Vec<_> = input
                .iter()
                .flat_map(|r| slow_bits(r.start, r.width, r.coeff[lane]))
                .collect();
            assert_eq!(view.lane(lane).unwrap().collect::<Vec<_>>(), expect);
        }
        assert_eq!(view.lane(2).err(), Some(AffineError::Index));
        assert_eq!(
            SortedBitRuns::new(P, &[input[1], input[1]]).err(),
            Some(AffineError::Overlap)
        );
        assert_eq!(
            SortedBitRuns::new(
                P,
                &[BitRun {
                    start: u32::MAX,
                    width: 2,
                    coeff: [3]
                }]
            )
            .err(),
            Some(AffineError::Layout)
        );
        assert_eq!(
            SortedBitRuns::new(
                P,
                &[BitRun {
                    start: 1,
                    width: 2,
                    coeff: [P]
                }]
            )
            .err(),
            Some(AffineError::Coefficient)
        );
    }

    #[test]
    fn strided_digits_are_canonical_without_sorting() {
        let fields = [[P - 9, 0], [11, P - 4], [0, 17]];
        let layout = StridedLayout {
            start: u32::MAX - 191,
            group_stride: 64,
            block_stride: 8,
            blocks: 8,
            width: 4,
        };
        let a = StridedBitRuns::new(P, layout, &fields).unwrap();
        for lane in 0..2 {
            let expect: Vec<_> = fields
                .iter()
                .enumerate()
                .flat_map(|(i, v)| {
                    (0..8).flat_map(move |j| {
                        slow_bits(layout.start + 64 * i as u32 + 8 * j, 4, v[lane])
                    })
                })
                .collect();
            assert_eq!(a.lane(lane).unwrap().collect::<Vec<_>>(), expect);
        }
        assert_eq!(
            StridedBitRuns::new(
                P,
                StridedLayout {
                    group_stride: 30,
                    ..layout
                },
                &fields
            )
            .err(),
            Some(AffineError::Overlap)
        );
    }

    #[test]
    fn fused_radix_matches_separate_planes() {
        let bits = [!0u64 << 5, 0x5000_0000_0000_002a, !0u64 >> 11];
        let support = SupportRank::new(&bits).unwrap();
        let ids: Vec<u32> = (0..192)
            .filter(|&i| support.position(i).is_some())
            .collect();
        assert_eq!(ids.len(), support.len());
        let spans = [
            SumBits {
                start: 130,
                width: 17,
            },
            SumBits {
                start: 147,
                width: 17,
            },
            SumBits {
                start: 164,
                width: 17,
            },
            SumBits {
                start: 33,
                width: 28,
            },
        ];
        let weights = [
            [Q - 1, 2, Q - 3, 0],
            [17, 3, 0, Q - 17],
            [Q - 310, Q - 1, 0, 1],
            [1, Q - 1, Q - 1, Q - 1],
        ];
        let helper = Radix256::new(P, Q, weights).unwrap();
        let mut out = vec![[0; 4]; ids.len()];
        let mut expect = out.clone();
        for batch in 1..=5u32 {
            let terms: Vec<_> = ids
                .iter()
                .step_by(batch as usize)
                .map(|&bit| BitCoeff {
                    bit,
                    coeff: P - bit - batch,
                })
                .collect();
            helper
                .add_terms(&support, terms.iter().copied(), &mut out)
                .unwrap();
            helper.subtract_sums(&support, spans, &mut out).unwrap();
            for (plane, span) in spans.iter().enumerate() {
                for row in &terms {
                    let i = ids.binary_search(&row.bit).unwrap();
                    let c = u64::from(row.coeff.to_le_bytes()[plane]);
                    for lane in 0..4 {
                        expect[i][lane] = (expect[i][lane] + c * weights[plane][lane]) % Q;
                    }
                }
                for k in 0..span.width {
                    let i = ids.binary_search(&(span.start + k)).unwrap();
                    for lane in 0..4 {
                        let x = ((u128::from(weights[plane][lane]) << k) % u128::from(Q)) as u64;
                        expect[i][lane] = (expect[i][lane] + Q - x) % Q;
                    }
                }
            }
            assert_eq!(out, expect);
        }
        assert_eq!(
            helper.add_terms(&support, [BitCoeff { bit: 0, coeff: 1 }], &mut out),
            Err(AffineError::UnsupportedBit)
        );
        assert_eq!(
            Radix256::new(P, u64::MAX, weights).err(),
            Some(AffineError::Modulus)
        );
    }

    #[test]
    fn carry_bound_is_strict_and_uses_integer_lifts() {
        assert_eq!(
            no_wrap_residual(
                Q,
                1,
                [
                    BoundedTerm {
                        coeff: 255,
                        min: 0,
                        max: 1
                    },
                    BoundedTerm {
                        coeff: -256,
                        min: -8,
                        max: 7
                    },
                ]
            ),
            Ok((-1791, 2304))
        );
        for c in [-i128::from(Q), i128::from(Q)] {
            assert_eq!(no_wrap_residual(Q, c, []), Err(AffineError::WrapPossible));
        }
        assert_eq!(
            no_wrap_residual(
                Q,
                0,
                [BoundedTerm {
                    coeff: i128::MAX,
                    min: 0,
                    max: 2
                }]
            ),
            Err(AffineError::Overflow)
        );
    }
}
