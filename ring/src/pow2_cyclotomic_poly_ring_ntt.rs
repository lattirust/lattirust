use ark_std::fmt::{Debug, Display, Formatter};
use ark_std::hash::Hash;
use ark_std::io::{Read, Write};
use ark_std::iter::{Product, Sum};
use ark_std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use ark_std::rand::Rng;
use ark_std::UniformRand;
use derive_more::{From, Into};
use lattirust_linear_algebra::SVector;
use num_traits::{One, Zero};

use crate::pow2_cyclotomic_poly_ring::CyclotomicPolyRingGeneral;
use crate::ring_config::{Pow2Rp64Config, RpConfig};
use crate::traits::{FromRandomBytes, WithL2Norm, WithLinfNorm};
use crate::PolyRing;
use crate::Ring;

use super::poly_ring::WithRot;

use ark_ff::{Field, Fp, Fp64, FpConfig};

/// A cyclotomic ring of the form Fp<C,N>[X]/<Phi_D(X)>
/// which fully splits in its CRT representation.
#[derive(From, Into)]
pub struct CyclotomicPolyRingNTTGeneral<C: RpConfig<N>, const N: usize, const PHI_D: usize>(
    SVector<Fp<C::FpConfig, N>, PHI_D>,
);

pub type Fp64Pow2<const Q: u64, const PHI_D: usize> =
    Fp64<<Pow2Rp64Config<Q, PHI_D> as RpConfig<1>>::FpConfig>; // This looks ugly change this

/// A cyclotomic ring with a cyclotomic polynomial degree of a power of two
/// and a modulus less than 2^64
pub type Pow2CyclotomicPolyRingNTT<const Q: u64, const PHI_D: usize> =
    CyclotomicPolyRingNTTGeneral<Pow2Rp64Config<Q, PHI_D>, 1, PHI_D>;

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> CyclotomicPolyRingNTTGeneral<C, N, PHI_D> {
    /// Constructs a polynomial from a function specifying coefficients in non-NTT form.
    fn from_coeffs_vec(coeffs: Vec<Fp<C::FpConfig, N>>) -> Self {
        let mut coeffs = coeffs;
        if coeffs.len() > PHI_D {
            C::reduce_in_place(&mut coeffs);
        }
        assert_eq!(
            coeffs.len(),
            PHI_D,
            "Reduce in place not reducing to correct length"
        );
        C::crt_in_place(&mut coeffs);
        Self::from_array(coeffs.try_into().map_err(|v: Vec<_>| v).unwrap())
    }
    fn from_fn<F>(f: F) -> Self
    where
        F: FnMut(usize) -> Fp<C::FpConfig, N>,
    {
        let mut coeffs = Vec::from(core::array::from_fn::<_, PHI_D, _>(f));
        assert_eq!(
            coeffs.len(),
            PHI_D,
            "Incorrect length of the coefficient vector"
        );
        C::crt_in_place(&mut coeffs);
        Self(
            coeffs
                .try_into()
                .map_err(|v: Vec<Fp<C::FpConfig, N>>| v)
                .unwrap(),
        )
    }

    fn from_array(ntt_coeffs: [Fp<C::FpConfig, N>; PHI_D]) -> Self {
        Self(SVector::const_from_array(ntt_coeffs))
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> PartialEq
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Eq
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Clone
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Copy
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Debug
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> ark_std::fmt::Result {
        Debug::fmt(&self.0, f)
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Display
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> ark_std::fmt::Result {
        write!(f, "CyclotomicPolyRingNTTGeneral(")?;
        let mut iter = self.0.iter();
        if let Some(first) = iter.next() {
            write!(f, "{}", first)?;
            for field_element in iter {
                write!(f, ", {}", field_element)?;
            }
        }
        write!(f, ")")
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Hash
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn hash<H: ark_std::hash::Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

const fn vec_from_element<FP: FpConfig<N>, const N: usize, const PHI_D: usize>(
    elem: Fp<FP, N>,
) -> SVector<Fp<FP, N>, PHI_D> {
    SVector::const_from_array([elem; PHI_D])
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> CanonicalSerialize
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn serialize_with_mode<W: Write>(
        &self,
        writer: W,
        compress: Compress,
    ) -> Result<(), SerializationError> {
        self.0.serialize_with_mode(writer, compress)
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        self.0.serialized_size(compress)
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Valid
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn check(&self) -> Result<(), SerializationError> {
        self.0.check()
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> CanonicalDeserialize
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn deserialize_with_mode<R: Read>(
        reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        SVector::<Fp<C::FpConfig, N>, PHI_D>::deserialize_with_mode(reader, compress, validate)
            .map(Self)
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Ring
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    const ZERO: Self = Self(vec_from_element(<Fp<C::FpConfig, N> as Field>::ZERO));
    const ONE: Self = Self(vec_from_element(<Fp<C::FpConfig, N> as Field>::ONE));
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> FromRandomBytes<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn byte_size() -> usize {
        PHI_D * Fp::<C::FpConfig, N>::byte_size()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        assert_eq!(bytes.len(), Self::byte_size());
        let coeffs = core::array::from_fn(|i| {
            Fp::<C::FpConfig, N>::try_from_random_bytes(
                &bytes[i * Fp::<C::FpConfig, N>::byte_size()
                    ..(i + 1) * Fp::<C::FpConfig, N>::byte_size()],
            )
            .unwrap()
        });
        Some(Self::from_array(coeffs))
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Default
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    #[inline(always)]
    fn default() -> Self {
        Self::zero()
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Zero
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    #[inline(always)]
    fn zero() -> Self {
        Self::ZERO
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.eq(&Self::ZERO)
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> One
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    #[inline(always)]
    fn one() -> Self {
        Self::ONE
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Mul<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0.component_mul(&rhs.0))
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Neg
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(self.0.neg())
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> UniformRand
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self::from_fn(|_| Fp::<C::FpConfig, N>::rand(rng))
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> MulAssign<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn mul_assign(&mut self, rhs: Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> Add<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> Sub<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> Add<&'a mut Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn add(self, rhs: &'a mut Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> Sub<&'a mut Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn sub(self, rhs: &'a mut Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> Mul<&'a mut Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn mul(self, rhs: &'a mut Self) -> Self::Output {
        Self(self.0.component_mul(&rhs.0))
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> AddAssign<&'a mut Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn add_assign(&mut self, rhs: &'a mut Self) {
        self.0 += rhs.0;
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> SubAssign<&'a mut Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn sub_assign(&mut self, rhs: &'a mut Self) {
        self.0 -= rhs.0;
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> MulAssign<&'a mut Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn mul_assign(&mut self, rhs: &'a mut Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize>
    From<CyclotomicPolyRingGeneral<C, N, PHI_D>> for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn from(value: CyclotomicPolyRingGeneral<C, N, PHI_D>) -> Self {
        let mut coeffs: Vec<Fp<C::FpConfig, N>> = value.coeffs();
        C::crt_in_place(&mut coeffs);
        Self(
            coeffs
                .try_into()
                .map_err(|v: Vec<Fp<C::FpConfig, N>>| v)
                .unwrap(),
        )
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Mul<Fp<C::FpConfig, N>>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn mul(self, rhs: Fp<C::FpConfig, N>) -> Self::Output {
        self.mul(Self::from_scalar(rhs))
    }
}

macro_rules! impl_from_primitive_type {
    ($primitive_type: ty) => {
        impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> From<$primitive_type>
            for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
        {
            fn from(value: $primitive_type) -> Self {
                Self::from_scalar(Fp::<C::FpConfig, N>::from(value))
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

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> Mul<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        Self(self.0.component_mul(&rhs.0))
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> AddAssign<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn add_assign(&mut self, rhs: &'a Self) {
        self.0 += rhs.0;
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> SubAssign<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn sub_assign(&mut self, rhs: &'a Self) {
        self.0 -= rhs.0;
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> MulAssign<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn mul_assign(&mut self, rhs: &'a Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Add<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Sub<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> AddAssign<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> SubAssign<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Sum<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<'a, C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Sum<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Product<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<'a, C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Product<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> PolyRing
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    type BaseRing = Fp<C::FpConfig, N>;
    fn coeffs(&self) -> Vec<Fp<C::FpConfig, N>> {
        let mut evaluations = self.0.as_slice().to_vec();
        C::icrt_in_place(&mut evaluations);
        if evaluations.is_empty() {
            // When coeffs is zero it returns an empty vec
            evaluations = vec![<Fp<C::FpConfig, N> as Field>::ZERO];
        }
        evaluations
    }
    fn dimension() -> usize {
        PHI_D
    }

    fn from_scalar(v: Self::BaseRing) -> Self {
        // NTT([v, 0, ..., 0]) = ([v, ..., v])
        Self::from_array([v; PHI_D])
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> From<Vec<Fp<C::FpConfig, N>>>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn from(mut value: Vec<Fp<C::FpConfig, N>>) -> Self {
        value.resize_with(PHI_D, Fp::zero);
        C::crt_in_place(&mut value);
        Self(
            value
                .try_into()
                .map_err(|v: Vec<Fp<C::FpConfig, N>>| v)
                .unwrap(),
        )
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> From<Fp<C::FpConfig, N>>
    for CyclotomicPolyRingNTTGeneral<C, N, PHI_D>
{
    fn from(value: Fp<C::FpConfig, N>) -> Self {
        Self::from_scalar(value)
    }
}

// impl<const Q: u64, const N: usize> WithConjugationAutomorphism for Pow2CyclotomicPolyRingNTT<Q, N> {
//     fn sigma(&self) -> Self {
//         // TODO: can we implement the automorphism directly in NTT form?
//         Into::<Pow2CyclotomicPolyRing<Zq<Q>, N>>::into(*self)
//             .sigma()
//             .into()
//     }
// }

impl<const Q: u64, const PHI_D: usize> WithL2Norm for Pow2CyclotomicPolyRingNTT<Q, PHI_D> {
    fn l2_norm_squared(&self) -> u128 {
        self.coeffs().l2_norm_squared()
    }
}

impl<const Q: u64, const PHI_D: usize> WithLinfNorm for Pow2CyclotomicPolyRingNTT<Q, PHI_D> {
    fn linf_norm(&self) -> u128 {
        self.coeffs().linf_norm()
    }
}

impl<const Q: u64, const PHI_D: usize> WithRot for Pow2CyclotomicPolyRingNTT<Q, PHI_D> {
    fn multiply_by_xi(&self, i: usize) -> Vec<Self::BaseRing> {
        let mut xi = if i < PHI_D {
            vec![<Fp64Pow2<Q, PHI_D> as Field>::ZERO; PHI_D]
        } else {
            vec![<Fp64Pow2<Q, PHI_D> as Field>::ZERO; i + 1] // Shouldn't get to this point for the
                                                             // purpose of RotSum
        };
        xi[i] = <Fp64Pow2<Q, PHI_D> as Field>::ONE;

        let xi_poly = Pow2CyclotomicPolyRingNTT::from_coeffs_vec(xi);
        let result = (*self * xi_poly).0;
        result.iter().copied().collect::<Vec<_>>()
    }
}

#[cfg(test)]
mod tests {
    use crate::{z_q::FqConfig, PolyRing, Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT};
    use ark_ff::{Fp, MontBackend};
    use ark_std::UniformRand;
    use lattirust_linear_algebra::SVector;
    use rand::thread_rng;
    const FERMAT_Q: u64 = (1 << 16) + 1;
    type FermatFqConfig = FqConfig<FERMAT_Q>;

    const FERMAT_NS: [usize; 14] = [
        1 << 1,
        1 << 2,
        1 << 3,
        1 << 4,
        1 << 5,
        1 << 6,
        1 << 7,
        1 << 8,
        1 << 9,
        1 << 10,
        1 << 11,
        1 << 12,
        1 << 13,
        1 << 14,
    ];
    macro_rules! generate_ntt_form_match {
        ($n:expr, $initial_coeffs:expr, $($size:expr),+) => {
            match $n {
                $(
                    $size => {
                        let ntt_form: Pow2CyclotomicPolyRingNTT<FERMAT_Q, $size> =
                                Pow2CyclotomicPolyRingNTT::from_coeffs_vec($initial_coeffs.clone());
                        let intt_coeffs = ntt_form.coeffs();
                        assert_eq!($initial_coeffs, intt_coeffs);
                    },
                )+
                _ => unreachable!("Unsupported N value: {}", $n),
            }
        };
    }

    #[test]
    fn test_ntt_pow2() {
        let mut rng = thread_rng();
        for N in FERMAT_NS {
            let initial_coeffs = (0..N)
                .map(|_| Fp::<MontBackend<FermatFqConfig, 1>, 1>::rand(&mut rng))
                .collect::<Vec<_>>();
            generate_ntt_form_match!(
                N,
                initial_coeffs,
                2,
                4,
                8,
                16,
                32,
                64,
                128,
                256,
                512,
                1024,
                2048,
                4096,
                8192,
                16384
            );
        }
    }

    fn test_mul_ntt_pow2<const SIZE: usize>() {
        let mut rng = thread_rng();

        let coeff_1 = Pow2CyclotomicPolyRing::<FERMAT_Q, SIZE>::rand(&mut rng);
        let coeff_2 = Pow2CyclotomicPolyRing::<FERMAT_Q, SIZE>::rand(&mut rng);

        let ntt_form_1 = Pow2CyclotomicPolyRingNTT::from(coeff_1);
        let ntt_form_2 = Pow2CyclotomicPolyRingNTT::from(coeff_2);

        let ntt_mul = ntt_form_1 * ntt_form_2;
        let coeffs_mul = coeff_1 * coeff_2;
        // ntt_mul.coeffs() performs INTT while coeffs_mul.coeffs() just returns the coefficients
        assert_eq!(ntt_mul.coeffs(), coeffs_mul.coeffs());
    }

    #[test]
    fn test_mul_ntt_pow2_hardcoded() {
        let coeffs_1_vec_vec = [1, 2, 3, 4, 5, 6, 7, 8];
        let coeffs_1_vec = coeffs_1_vec_vec
            .into_iter()
            .map(Fp::<MontBackend<FermatFqConfig, 1>, 1>::from)
            .collect::<Vec<_>>();
        let mut coeffs_1 = Pow2CyclotomicPolyRing::<FERMAT_Q, 8>::from_coeffs_vec(coeffs_1_vec);
        let coeffs_2_vec_vec = [8, 7, 6, 5, 4, 3, 2, 1];
        let coeffs_2_vec = coeffs_2_vec_vec
            .into_iter()
            .map(Fp::<MontBackend<FermatFqConfig, 1>, 1>::from)
            .collect::<Vec<_>>();
        let coeffs_2 = Pow2CyclotomicPolyRing::<FERMAT_Q, 8>::from_coeffs_vec(coeffs_2_vec);

        let mut ntt_form_1 = Pow2CyclotomicPolyRingNTT::from(coeffs_1);
        ntt_form_1 *= Pow2CyclotomicPolyRingNTT::from(coeffs_2);

        coeffs_1 *= coeffs_2;
        assert_eq!(ntt_form_1.coeffs(), coeffs_1.coeffs());
    }
    // TODO: test mutable mul, i.e., mul_assign

    #[test]
    fn test_mul_ntt_pow2_multiple_sizes() {
        for size in [2, 4, 8, 16, 32] {
            match size {
                2 => test_mul_ntt_pow2::<2>(),
                4 => test_mul_ntt_pow2::<4>(),
                8 => test_mul_ntt_pow2::<8>(),
                16 => test_mul_ntt_pow2::<16>(),
                32 => test_mul_ntt_pow2::<32>(),
                _ => unreachable!("Unsupported size"),
            }
        }
    }

    fn test_ntt_pow2_hardcoded<const SIZE: usize>(
        coeffs: [u64; SIZE],
        expected_values: [u64; SIZE],
    ) {
        let coeffs = coeffs
            .into_iter()
            .map(Fp::<MontBackend<FermatFqConfig, 1>, 1>::from)
            .collect::<Vec<_>>();
        let expected_ntt = expected_values
            // Numbers obtained from Python library
            .into_iter()
            .map(Fp::<MontBackend<FermatFqConfig, 1>, 1>::from)
            .collect::<Vec<_>>();

        let ntt_form: Pow2CyclotomicPolyRingNTT<FERMAT_Q, SIZE> =
            Pow2CyclotomicPolyRingNTT::from_coeffs_vec(coeffs);

        let expected = SVector::const_from_array(expected_ntt.try_into().unwrap());
        assert_eq!(ntt_form.0, expected);
    }

    #[test]
    fn test_ntt_pow2_example_1() {
        test_ntt_pow2_hardcoded::<8>(
            [1, 2, 3, 4, 5, 6, 7, 8],
            [52195, 23595, 14578, 40635, 35584, 36407, 23601, 35561],
        );
    }

    #[test]
    fn test_ntt_pow2_example_2() {
        test_ntt_pow2_hardcoded::<16>(
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
            [
                1620, 16633, 40, 14048, 49787, 36415, 21641, 52183, 33435, 38424, 36421, 55728,
                54235, 65523, 46545, 1634,
            ],
        );
    }

    #[test]
    fn test_ntt_pow2_example_3() {
        test_ntt_pow2_hardcoded::<32>(
            [
                1182, 3320, 5933, 6237, 10981, 11828, 12004, 15261, 15742, 18544, 20536, 21395,
                22087, 22505, 22654, 23275, 24128, 33816, 40316, 41665, 41921, 41956, 42300, 46980,
                47250, 49871, 53569, 54471, 63813, 64688, 65223, 65455,
            ],
            [
                10587, 59994, 64892, 33848, 20090, 62609, 20516, 29481, 19130, 54823, 54845, 11706,
                56764, 10266, 60708, 1414, 43951, 59219, 20434, 2363, 55942, 26285, 6855, 31089,
                60856, 6403, 50979, 20512, 15763, 16492, 47941, 49659,
            ],
        );
    }
}
