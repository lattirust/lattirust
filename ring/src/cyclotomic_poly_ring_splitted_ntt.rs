use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Valid};
use ark_std::fmt::{Debug, Display};
use ark_std::hash::Hash;
use ark_std::iter::{Product, Sum};
use ark_std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use ark_std::UniformRand;
use derive_more::{Add, AddAssign, From, Sub, SubAssign, Sum};
use num_traits::{One, Zero};

use crate::partial_ntt::PartialNTT;
use crate::traits::{FromRandomBytes, WithL2Norm, WithLinfNorm};
use crate::{OverField, PolyRing};
use crate::{Ring, Zq};
use lattirust_linear_algebra::SVector;

use super::poly_ring::WithRot;

/// For an m-th cyclotomic ring its minimal polynomail has the form
/// `Phi_m = Prod_{j=1}^phi(m)(x - r_j)`
/// Where `ord_p(r_j) = m` and `r_j = w^i` where `i \in {i: gcd(i,m) = 1}`
/// and `j` is its enumaration
/// Note that this requires that `p = 1 mod m`
///
/// Another option is to split the minimal polynomial as
/// `Phi_m = Prod_{j = 1}^phi(z) (x^{m/z} - r_j)`
/// Where `ord_p(r_j) = z` and `r_j = w^i` where `i \in {i: gcd(i,z) = 1}`
/// this requires that `p = q mod z` so `m` doesn't have to divide `p-1` but
/// z has to share the prime factorizations with diferent powers
/// For more see:
///     Short invertible elements in partially splitting cyclotomic rings
///
/// We have that:
/// R = DirProd R_j and R_j \in Z_p[X]/ (X^{m/z} - r_j)
/// This struct can be thought as the concatenation of R_j elements
/// where the operations over each R_j are element-wise
/// let `D` = m/z so we have phi(z) R_j components each with D components
/// also phi(z) = N/D
/// For more details look at [Short, Invertible Elements in Partially Splitting Cyclotomic Rings and Applications to Lattice-Based Zero-Knowledge Proofs](https://eprint.iacr.org/2017/523)
#[derive(PartialEq, Eq, Debug, Clone, Copy, Hash, Add, AddAssign, Sum, Sub, SubAssign, From)]
pub struct CyclotomicPolyRingSplittedNTT<
    const Q: u64,
    const ROU: u64,
    const N: usize,
    const D: usize,
    const Z: usize,
    const PHI_Z: usize,
>(SVector<Zq<Q>, N>);

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    #[allow(dead_code)]
    pub(crate) type Inner = SVector<Zq<Q>, N>;

    /// Constructs a polynomial from an array of coefficients in NTT form.
    /// This function is private since we can't enforce coeffs being in NTT form if called from outside.
    fn from_array(coeffs_ntt: [Zq<Q>; N]) -> Self {
        Self(Self::Inner::const_from_array(coeffs_ntt))
    }

    pub fn from_fn<F>(f: F) -> Self
    where
        F: FnMut(usize) -> Zq<Q>,
    {
        let mut coeffs = core::array::from_fn(f);
        Self::ntt(&mut coeffs, Zq::<Q>::from(ROU));
        Self::from_array(coeffs)
    }

    /// `coeffs` are the coefficients of the complete Ring not partially splitted
    pub fn from_coeffs(coeffs: [Zq<Q>; N]) -> Self {
        let mut coeffs = coeffs;
        Self::ntt(&mut coeffs, Zq::<Q>::from(ROU));
        Self::from_array(coeffs)
    }

    pub fn ntt_mul(&self, rhs: &Self, rou: Zq<Q>) -> Self {
        let mut temp = [Zq::<Q>::from(0); N];
        let components = coprimes_set::<Z, PHI_Z>();
        let num_chunks = N / D;
        for k in 0..num_chunks {
            for i in 0..D {
                for j in 0..D - i {
                    temp[k * D + i + j] += self.0[i + k * D] * rhs.0[j + k * D];
                }
            }
            let rj_power = components[k] as u64;
            let rj = Ring::pow(&rou, [rj_power]);
            for i in 1..D {
                for j in D - i..D {
                    temp[k * D + i + j - D] += self.0[i + k * D] * rhs.0[j + k * D] * rj;
                }
            }
        }
        Self::from_array(temp)
    }
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > PartialNTT<Q, N, D, Z, PHI_Z> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
}

const fn vec_from_element<const Q: u64, const N: usize>(elem: Zq<Q>) -> SVector<Zq<Q>, N> {
    SVector::<Zq<Q>, N>::const_from_array([elem; N])
}

const fn coprimes_set<const Z: usize, const PHI_Z: usize>() -> [usize; PHI_Z] {
    let mut i = 0;
    let mut j = 0;
    let mut set = [0; PHI_Z];
    while i < Z {
        if gcd::<Z>(i) {
            set[j] = i;
            j += 1;
        }
        i += 1;
    }
    set
}

const fn gcd<const Z: usize>(i: usize) -> bool {
    let mut a = i;
    let mut b = Z;
    while b != 0 {
        let temp = b;
        b = a % b;
        a = temp
    }
    a == 1
}

// ============================================================================

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > CanonicalSerialize for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn serialize_with_mode<W: ark_std::io::prelude::Write>(
        &self,
        writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        self.0.serialize_with_mode(writer, compress)
    }

    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        self.0.serialized_size(compress)
    }
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Valid for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        self.0.check()
    }
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > CanonicalDeserialize for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn deserialize_with_mode<R: ark_std::io::prelude::Read>(
        _reader: R,
        _compress: ark_serialize::Compress,
        _validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        todo!()
    }
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Ring for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    const ZERO: Self = Self(vec_from_element(<Zq<Q> as Ring>::ZERO));
    const ONE: Self = Self(vec_from_element(<Zq<Q> as Ring>::ONE));
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > FromRandomBytes<Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn byte_size() -> usize {
        N * Zq::<Q>::byte_size()
    }
    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        assert_eq!(bytes.len(), Self::byte_size());
        let coeffs = core::array::from_fn(|i| {
            Zq::<Q>::try_from_random_bytes(
                &bytes[i * Zq::<Q>::byte_size()..(i + 1) * Zq::<Q>::byte_size()],
            )
            .unwrap()
        });
        Some(Self::from_array(coeffs))
    }
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Default for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    #[inline(always)]
    fn default() -> Self {
        Self::zero()
    }
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Display for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
        ark_std::fmt::Display::fmt(&self.0, f)
    }
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Zero for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn zero() -> Self {
        Self::ZERO
    }
    fn is_zero(&self) -> bool {
        self.eq(&Self::ZERO)
    }
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > One for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn one() -> Self {
        Self::ONE
    }
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Mul<Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        self.ntt_mul(&rhs, Zq::<Q>::from(ROU))
    }
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Neg for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    type Output = Self;
    fn neg(self) -> Self::Output {
        self.0.neg().into()
    }
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > UniformRand for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn rand<R: rand::prelude::Rng + ?Sized>(rng: &mut R) -> Self {
        Self::from_fn(|_| Zq::<Q>::rand(rng))
    }
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > MulAssign<Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn mul_assign(&mut self, rhs: Self) {
        let res = self.ntt_mul(&rhs, Zq::<Q>::from(ROU));
        *self += res;
    }
}

impl<
        'a,
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Add<&'a Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    type Output = Self;
    fn add(self, rhs: &'a Self) -> Self::Output {
        self.0.add(rhs.0).into()
    }
}

impl<
        'a,
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Sub<&'a Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    type Output = Self;
    fn sub(self, rhs: &'a Self) -> Self::Output {
        self.0.sub(rhs.0).into()
    }
}

impl<
        'a,
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Mul<&'a Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    type Output = Self;
    fn mul(self, rhs: &'a Self) -> Self::Output {
        self.ntt_mul(rhs, Zq::<Q>::from(ROU))
    }
}

impl<
        'a,
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > AddAssign<&'a Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn add_assign(&mut self, rhs: &'a Self) {
        self.0.add_assign(rhs.0)
    }
}

impl<
        'a,
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > SubAssign<&'a Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn sub_assign(&mut self, rhs: &'a Self) {
        self.0.sub_assign(rhs.0)
    }
}

impl<
        'a,
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > MulAssign<&'a Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn mul_assign(&mut self, rhs: &'a Self) {
        let res = self.ntt_mul(rhs, Zq::<Q>::from(ROU));
        self.0.add_assign(res.0)
    }
}

impl<
        'a,
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Add<&'a mut Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    type Output = Self;
    fn add(self, rhs: &'a mut Self) -> Self::Output {
        self.0.add(rhs.0).into()
    }
}

impl<
        'a,
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Sub<&'a mut Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    type Output = Self;
    fn sub(self, rhs: &'a mut Self) -> Self::Output {
        self.0.sub(rhs.0).into()
    }
}

impl<
        'a,
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Mul<&'a mut Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    type Output = Self;
    fn mul(self, rhs: &'a mut Self) -> Self::Output {
        self.ntt_mul(rhs, Zq::<Q>::from(ROU))
    }
}

impl<
        'a,
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > AddAssign<&'a mut Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn add_assign(&mut self, rhs: &'a mut Self) {
        self.0.add_assign(rhs.0)
    }
}

impl<
        'a,
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > SubAssign<&'a mut Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn sub_assign(&mut self, rhs: &'a mut Self) {
        self.0.sub_assign(rhs.0)
    }
}

impl<
        'a,
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > MulAssign<&'a mut Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn mul_assign(&mut self, rhs: &'a mut Self) {
        let res = self.ntt_mul(rhs, Zq::<Q>::from(ROU));
        self.0.add_assign(res.0)
    }
}

// From<CyclotomicPolyRing>

// Into

impl<
        'a,
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Mul<Zq<Q>> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    type Output = Self;
    fn mul(self, rhs: Zq<Q>) -> Self::Output {
        let constant = Self::from_scalar(rhs);
        self.ntt_mul(&constant, Zq::<Q>::from(ROU))
    }
}

macro_rules! impl_from_primitive_type {
    ($primitive_type: ty) => {
        impl<
                const Q: u64,
                const ROU: u64,
                const N: usize,
                const D: usize,
                const Z: usize,
                const PHI_Z: usize,
            > From<$primitive_type> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
        {
            fn from(value: $primitive_type) -> Self {
                Self::from_scalar(Zq::<Q>::from(value))
            }
        }
    };
}

impl_from_primitive_type!(u128);
impl_from_primitive_type!(u64);
impl_from_primitive_type!(u32);
impl_from_primitive_type!(u16);
impl_from_primitive_type!(u8);
impl_from_primitive_type!(bool);

impl<
        'a,
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Sum<&'a Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Product<Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<
        'a,
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > Product<&'a Self> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > PolyRing for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    type BaseRing = Zq<Q>;
    fn coeffs(&self) -> Vec<Self::BaseRing> {
        // Do INTT
        let rou = Zq::<Q>::from(ROU);
        let zetas = coprimes_set::<Z, PHI_Z>()
            .iter()
            .map(|&i| Ring::pow(&rou, [i as u64]))
            .collect::<Vec<_>>();
        let minimal_poly = minimal_polynomial::<Q, N, D>(&zetas);

        let mut ms = zetas
            .iter()
            .map(|&zeta| polynomial_division::<Q, N, D>(&minimal_poly, zeta))
            .collect::<Vec<_>>();
        for (i, m) in ms.iter_mut().enumerate() {
            let mut zetas_copy = zetas.clone();
            let zeta_i = zetas_copy.remove(i);
            let mut y = <Zq<Q> as Ring>::ONE;
            for zeta_j in zetas_copy {
                y *= zeta_i - zeta_j;
            }

            y = y.inverse().unwrap();

            for i in 0..m.len() {
                m[i] *= y;
            }
        }
        let components = ms
            .iter()
            .enumerate()
            .map(|(i, m)| {
                let mut slice = vec![<Zq::<Q> as Ring>::ZERO; D];
                for j in 0..D {
                    slice[j] = self.0[i * D + j]
                }
                polynomial_multiplication::<Q>(m, &slice)
            })
            .collect::<Vec<_>>();
        let sum = sum_polys(components.as_slice());
        polynomial_division_remainder(&sum, &minimal_poly)
    }

    fn dimension() -> usize {
        N
    }

    fn from_scalar(scalar: Self::BaseRing) -> Self {
        let mut array = [Zq::<Q>::from(0); N];
        for i in 0..Z {
            array[i * D] = scalar;
        }
        Self::from_array(array)
    }
}

fn sum_polys<const Q: u64>(polys: &[Vec<Zq<Q>>]) -> Vec<Zq<Q>> {
    let max_degree = polys.iter().map(|p| p.len()).max().unwrap_or(0);

    let mut result = vec![<Zq::<Q> as Ring>::ZERO; max_degree];

    for poly in polys {
        for (i, &coeff) in poly.iter().enumerate() {
            result[i] += coeff;
        }
    }

    result
}

fn polynomial_multiplication<const Q: u64>(poly1: &Vec<Zq<Q>>, poly2: &Vec<Zq<Q>>) -> Vec<Zq<Q>> {
    let n = poly1.len();
    let m = poly2.len();
    let result_size = n + m - 1;
    let mut result = vec![<Zq::<Q> as Ring>::ZERO; result_size];

    for i in 0..n {
        for j in 0..m {
            result[i + j] += poly1[i] * poly2[j];
        }
    }

    result
}

fn minimal_polynomial<const Q: u64, const N: usize, const D: usize>(
    zetas: &Vec<Zq<Q>>,
) -> Vec<Zq<Q>> {
    let mut minimal_poly = vec![<Zq::<Q> as Ring>::ZERO; N + 1];
    minimal_poly[0] = <Zq<Q> as Ring>::ONE;
    for (_i, &zeta) in zetas.iter().enumerate() {
        let temp = minimal_poly.clone();
        minimal_poly.rotate_right(D);
        for j in 0..minimal_poly.len() {
            minimal_poly[j] -= zeta * temp[j]; //This can be optimize
        }
    }
    minimal_poly
}

fn polynomial_division<const Q: u64, const N: usize, const D: usize>(
    minimal_polynomial: &Vec<Zq<Q>>,
    zeta: Zq<Q>,
) -> Vec<Zq<Q>> {
    assert_eq!(N + 1, minimal_polynomial.len());
    let mut quotient = vec![<Zq::<Q> as Ring>::ZERO; N - D + 1];
    for i in N - 2 * D + 1..N - D + 1 {
        quotient[i] = minimal_polynomial[i + D];
    }

    for i in (0..N - 2 * D + 1).rev() {
        quotient[i] = minimal_polynomial[i + D] + (zeta * quotient[i + D]);
    }

    println!();
    quotient
}

fn polynomial_division_remainder<const Q: u64>(
    dividend: &[Zq<Q>],
    divisor: &[Zq<Q>], //Minimal poly
) -> Vec<Zq<Q>> {
    if divisor.is_empty() || divisor.iter().all(|&x| x == <Zq<Q> as Ring>::ZERO) {
        panic!("Division by zero polynomial")
    }

    let mut remainder = dividend.to_vec();
    let divisor_degree = divisor.len() - 1;
    let divided_degree = dividend.len() - 1;

    if divided_degree < divisor_degree {
        return remainder;
    }

    for i in (0..=divided_degree - divisor_degree).rev() {
        let factor = remainder[i + divisor_degree] / divisor[divisor_degree];
        for j in 0..=divisor_degree {
            remainder[i + j] -= factor * divisor[j];
        }
    }

    while remainder.len() > 1 && remainder.last() == Some(&<Zq<Q> as Ring>::ZERO) {
        remainder.pop();
    }
    remainder
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > From<Vec<Zq<Q>>> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn from(value: Vec<Zq<Q>>) -> Self {
        let mut array = TryInto::<[Zq<Q>; N]>::try_into(value).unwrap();
        Self::ntt(&mut array, Zq::<Q>::from(ROU));
        Self::from_array(array)
    }
}
impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > From<Zq<Q>> for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn from(value: Zq<Q>) -> Self {
        Self::from_scalar(value)
    }
}
// impl<
//         const Q: u64,
//         const ROU: u64,
//         const N: usize,
//         const D: usize,
//         const Z: usize,
//         const PHI_Z: usize,
//     > WithConjugationAutomorphism for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
// {
//     fn sigma(&self) -> Self {
//         // What is sigma?
//         todo!()
//     }
// }
impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > WithL2Norm for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn l2_norm_squared(&self) -> u128 {
        self.coeffs().l2_norm_squared()
    }
}

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > WithLinfNorm for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn linf_norm(&self) -> u128 {
        self.coeffs().linf_norm()
    }
}
impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > WithRot for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
    fn multiply_by_xi(&self, i: usize) -> Vec<Self::BaseRing> {
        let mut xi = [<Zq<Q> as Ring>::ZERO; N];
        if i < N {
            xi[i] = <Zq<Q> as Ring>::ONE;
        } else {
            unimplemented!("No support for multiplying with a polynomial bigger than N");
        }
        let xi_poly = Self::from_fn(|i| xi[i]);
        let result = *self * xi_poly;
        result.0.iter().copied().collect::<Vec<_>>()
    }
}
// ============================================================================

impl<
        const Q: u64,
        const ROU: u64,
        const N: usize,
        const D: usize,
        const Z: usize,
        const PHI_Z: usize,
    > OverField for CyclotomicPolyRingSplittedNTT<Q, ROU, N, D, Z, PHI_Z>
{
}

#[cfg(test)]
mod tests {
    use ark_std::UniformRand;
    use rand::rngs::StdRng;

    use crate::{
        cyclotomic_poly_ring_splitted_ntt::{
            minimal_polynomial, polynomial_division_remainder, polynomial_multiplication,
        },
        PolyRing, Ring, Zq,
    };

    use super::CyclotomicPolyRingSplittedNTT;

    const Q: u64 = 15 * (1 << 27) + 1;
    const ROU: u64 = 284861408; // rou of the multiplicative subgroup of F_Q is
                                // 76160998 then ROU = 76160998^{ (Q-1)/Z }
    const N: usize = 1 << 4;
    const Z: usize = 1 << 2;
    const D: usize = 1 << (5 - 2);
    const PHI_Z: usize = 1 << 1;

    fn test_rng_helper() -> StdRng {
        // Should be good enough for testing
        use rand::SeedableRng;
        // arbitrary seed
        let seed = [
            1, 0, 0, 0, 23, 0, 0, 0, 200, 1, 0, 0, 210, 30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0,
        ];
        rand::rngs::StdRng::from_seed(seed)
    }

    #[test]
    fn ntt_test() {
        assert!(Q % Z as u64 == 1);
        let mut rng = test_rng_helper();
        let initial_coeffs = (0..N).map(|_| Zq::<Q>::rand(&mut rng)).collect::<Vec<_>>();
        // let initial_coeffs = (0..N).map(|i| Zq::<Q>::from(i as u64)).collect::<Vec<_>>();

        let poly =
            CyclotomicPolyRingSplittedNTT::<Q, ROU, N, D, Z, PHI_Z>::from_fn(|i| initial_coeffs[i]);

        let intt_coeffs = poly.coeffs();
        assert_eq!(initial_coeffs, intt_coeffs, "INNT incorrect");
    }
    #[test]
    fn tailored_mul_test() {
        assert!(Q % Z as u64 == 1);
        // let resize_power = (Q - 1) / Z as u64;
        // let rou = Zq::<Q>::from(76160998 as u64).pow([resize_power]);
        let rou = Zq::<Q>::from(ROU);
        assert_eq!(
            Ring::pow(&rou, [Z as u64]),
            <Zq::<Q> as Ring>::ONE,
            "ROU^Z != 1"
        );
        let rou3 = Ring::pow(&rou, [3]);
        let poly1 = CyclotomicPolyRingSplittedNTT::<Q, ROU, N, D, Z, PHI_Z>::from_fn(|i| {
            Zq::<Q>::from(i as u64)
        });
        let mut red_poly = [Zq::<Q>::from(0); N];
        for i in 0..N / 2 {
            red_poly[i] = Zq::<Q>::from(i as u64) + Zq::<Q>::from(8 + i as u64) * rou;
        }
        for i in 0..N / 2 {
            red_poly[N / 2 + i] = Zq::<Q>::from(i as u64) + Zq::<Q>::from(8 + i as u64) * rou3;
        }
        let poly_reduced =
            CyclotomicPolyRingSplittedNTT::<Q, ROU, N, D, Z, PHI_Z>::from_array(red_poly);
        assert_eq!(poly_reduced, poly1, "Poly NTT incorrect");
        let poly2 = CyclotomicPolyRingSplittedNTT::<Q, ROU, N, D, Z, PHI_Z>::from_fn(|i| {
            Zq::<Q>::from((N - 1 - i) as u64)
        });
        let final_poly = poly1.ntt_mul(&poly2, rou);

        let mut poly1_split1 = [Zq::<Q>::from(0); D];
        let mut poly2_split1 = [Zq::<Q>::from(0); D];
        let mut poly1_split2 = [Zq::<Q>::from(0); D];
        let mut poly2_split2 = [Zq::<Q>::from(0); D];
        for i in 0..D {
            poly1_split1[i] = poly1.0[i];
            poly2_split1[i] = poly2.0[i];

            poly1_split2[i] = poly1.0[i + D];
            poly2_split2[i] = poly2.0[i + D];
        }

        let naive_mul_poly_split1 = poly_mul(&poly1_split1, &poly2_split1);
        let naive_mul_poly_split2 = poly_mul(&poly1_split2, &poly2_split2);
        let split1 = reduce(naive_mul_poly_split1, rou);
        let split2 = reduce(naive_mul_poly_split2, rou3);
        let mut naive_coeffs = [Zq::<Q>::from(0); N];
        for i in 0..D {
            naive_coeffs[i] = split1[i];
            naive_coeffs[i + D] = split2[i];
        }

        let naive_poly =
            CyclotomicPolyRingSplittedNTT::<Q, ROU, N, D, Z, PHI_Z>::from_array(naive_coeffs);
        assert_eq!(
            naive_poly, final_poly,
            "Splitted NTT multiplication incorrect"
        );
    }

    #[test]
    fn mul_test() {
        assert!(Q % Z as u64 == 1);
        // let resize_power = (Q - 1) / Z as u64;
        // let rou = Zq::<Q>::from(76160998 as u64).pow([resize_power]);
        let rou = Zq::<Q>::from(ROU);
        assert_eq!(
            Ring::pow(&rou, [Z as u64]),
            <Zq::<Q> as Ring>::ONE,
            "ROU^Z != 1"
        );
        let mut rng = test_rng_helper();
        let coeffs1 = (0..N).map(|_| Zq::<Q>::rand(&mut rng)).collect::<Vec<_>>();
        let poly1 =
            CyclotomicPolyRingSplittedNTT::<Q, ROU, N, D, Z, PHI_Z>::from_fn(|i| coeffs1[i]);
        let coeffs2 = (0..N).map(|_| Zq::<Q>::rand(&mut rng)).collect::<Vec<_>>();
        let poly2 =
            CyclotomicPolyRingSplittedNTT::<Q, ROU, N, D, Z, PHI_Z>::from_fn(|i| coeffs2[i]);
        let final_poly = poly1 * poly2;

        let true_coeffs = polynomial_multiplication(&coeffs1, &coeffs2);
        let zetas = vec![rou, Ring::pow(&rou, [3])];
        let minimal_poly = minimal_polynomial::<Q, N, D>(&zetas);
        let true_coeffs = polynomial_division_remainder(&true_coeffs, &minimal_poly);
        let final_coeffs = final_poly.coeffs();
        assert_eq!(true_coeffs, final_coeffs);
    }

    fn poly_mul(a: &[Zq<Q>; D], b: &[Zq<Q>; D]) -> [Zq<Q>; 2 * D] {
        let mut temp = [Zq::<Q>::from(0); 2 * D];
        for i in 0..D {
            for j in 0..D {
                temp[i + j] += a[i] * b[j];
            }
        }
        temp
    }

    fn reduce(a: [Zq<Q>; 2 * D], rou: Zq<Q>) -> [Zq<Q>; D] {
        let mut temp = [Zq::<Q>::from(0); D];
        for i in 0..D {
            temp[i] += a[i];
        }
        for i in 0..D {
            temp[i] += a[i + D] * rou;
        }
        temp
    }
}
