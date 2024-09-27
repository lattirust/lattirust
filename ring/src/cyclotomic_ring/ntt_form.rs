use ark_ff::{Field, Fp};
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use ark_std::{
    fmt::{Debug, Display, Formatter},
    hash::Hash,
    io::{Read, Write},
    iter::{Product, Sum},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    rand::Rng,
    One, UniformRand, Zero,
};
use derive_more::{From, Into};

use super::{ring_config::CyclotomicConfig, CyclotomicPolyRingGeneral};
use crate::{traits::FromRandomBytes, PolyRing, Ring};
use lattirust_linear_algebra::SVector;

/// A cyclotomic ring Fp[X]/(Phi_m(X)) in the CRT-form.
/// * `C` is the configuration of the cyclotomic ring.
/// * `N` is the byte size of the underlying prime field.
/// * `D` is the number of factors in the CRT-representation of the ring.
#[derive(From, Into)]
pub struct CyclotomicPolyRingNTTGeneral<C: CyclotomicConfig<N>, const N: usize, const D: usize>(
    pub(crate) SVector<C::BaseCRTField, D>,
);

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> CyclotomicPolyRingNTTGeneral<C, N, D> {
    /// Constructs a polynomial from a function specifying coefficients in non-NTT form.
    fn from_coeffs_vec(coeffs: Vec<Fp<C::BaseFieldConfig, N>>) -> Self {
        let mut coeffs = coeffs;

        if coeffs.len() > D * C::CRT_FIELD_EXTENSION_DEGREE {
            C::reduce_in_place(&mut coeffs);
        } else {
            coeffs.resize(D * C::CRT_FIELD_EXTENSION_DEGREE, Fp::zero());
        }

        let evaluations = C::crt(coeffs);
        Self::from_array(
            evaluations
                .try_into()
                .expect("the output of CRT has incorrect length"),
        )
    }

    fn from_fn<F>(f: F) -> Self
    where
        F: FnMut(usize) -> Fp<C::BaseFieldConfig, N>,
    {
        let coeffs = Vec::from(core::array::from_fn::<_, D, _>(f));

        let evaluations = C::crt(coeffs);

        Self(
            evaluations
                .try_into()
                .expect("the output of CRT has incorrect length"),
        )
    }

    fn from_array(ntt_coeffs: [C::BaseCRTField; D]) -> Self {
        Self(SVector::const_from_array(ntt_coeffs))
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> PartialEq
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Eq
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Clone
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Copy
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Debug
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> ark_std::fmt::Result {
        Debug::fmt(&self.0, f)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Display
    for CyclotomicPolyRingNTTGeneral<C, N, D>
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

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Hash
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn hash<H: ark_std::hash::Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

const fn vec_from_element<F: Field, const D: usize>(elem: F) -> SVector<F, D> {
    SVector::const_from_array([elem; D])
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> CanonicalSerialize
    for CyclotomicPolyRingNTTGeneral<C, N, D>
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

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Valid
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn check(&self) -> Result<(), SerializationError> {
        self.0.check()
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> CanonicalDeserialize
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn deserialize_with_mode<R: Read>(
        reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        SVector::<C::BaseCRTField, D>::deserialize_with_mode(reader, compress, validate).map(Self)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Ring
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    const ZERO: Self = Self(vec_from_element(<C::BaseCRTField as Field>::ZERO));
    const ONE: Self = Self(vec_from_element(<C::BaseCRTField as Field>::ONE));
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> FromRandomBytes<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn byte_size() -> usize {
        D * C::BaseCRTField::byte_size()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        assert_eq!(bytes.len(), Self::byte_size());

        let coeffs = core::array::from_fn(|i| {
            C::BaseCRTField::try_from_random_bytes(
                &bytes[i * C::BaseCRTField::byte_size()..(i + 1) * C::BaseCRTField::byte_size()],
            )
            .unwrap()
        });
        Some(Self::from_array(coeffs))
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Default
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    #[inline(always)]
    fn default() -> Self {
        Self::zero()
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Zero
    for CyclotomicPolyRingNTTGeneral<C, N, D>
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

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> One
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    #[inline(always)]
    fn one() -> Self {
        Self::ONE
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Mul<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0.component_mul(&rhs.0))
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Neg
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(self.0.neg())
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> UniformRand
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self::from_fn(|_| Fp::<C::BaseFieldConfig, N>::rand(rng))
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> MulAssign<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn mul_assign(&mut self, rhs: Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Add<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Sub<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Add<&'a mut Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn add(self, rhs: &'a mut Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Sub<&'a mut Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn sub(self, rhs: &'a mut Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Mul<&'a mut Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn mul(self, rhs: &'a mut Self) -> Self::Output {
        Self(self.0.component_mul(&rhs.0))
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> AddAssign<&'a mut Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn add_assign(&mut self, rhs: &'a mut Self) {
        self.0 += rhs.0;
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> SubAssign<&'a mut Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn sub_assign(&mut self, rhs: &'a mut Self) {
        self.0 -= rhs.0;
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> MulAssign<&'a mut Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn mul_assign(&mut self, rhs: &'a mut Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize>
    From<CyclotomicPolyRingGeneral<C, N, { D * C::CRT_FIELD_EXTENSION_DEGREE }>>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn from(value: CyclotomicPolyRingGeneral<C, N, { D * C::CRT_FIELD_EXTENSION_DEGREE }>) -> Self {
        Self::from_coeffs_vec(value.coeffs())
    }
}

macro_rules! impl_from_primitive_type {
    ($primitive_type: ty) => {
        impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> From<$primitive_type>
            for CyclotomicPolyRingNTTGeneral<C, N, D>
        {
            fn from(value: $primitive_type) -> Self {
                Self::from_scalar(C::BaseCRTField::from_base_prime_field(Fp::<
                    C::BaseFieldConfig,
                    N,
                >::from(
                    value
                )))
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

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Mul<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        Self(self.0.component_mul(&rhs.0))
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> AddAssign<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn add_assign(&mut self, rhs: &'a Self) {
        self.0 += rhs.0;
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> SubAssign<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn sub_assign(&mut self, rhs: &'a Self) {
        self.0 -= rhs.0;
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> MulAssign<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn mul_assign(&mut self, rhs: &'a Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Add<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Sub<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> AddAssign<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> SubAssign<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Sum<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Sum<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Product<Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Product<&'a Self>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> PolyRing
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type BaseRing = C::BaseCRTField;

    fn coeffs(&self) -> Vec<C::BaseCRTField> {
        self.0.as_slice().into()
    }
    fn dimension() -> usize {
        D
    }

    fn from_scalar(v: Self::BaseRing) -> Self {
        // NTT([v, 0, ..., 0]) = ([v, ..., v])
        Self::from_array([v; D])
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> From<Vec<C::BaseCRTField>>
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn from(value: Vec<C::BaseCRTField>) -> Self {
        Self(value.try_into().expect("Should be of correct length"))
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::{Fp, MontBackend};
    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::{
        cyclotomic_ring::models::pow2_debug::{Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT},
        zn::z_q::FqConfig,
        PolyRing,
    };
    use lattirust_linear_algebra::SVector;

    const FERMAT_Q: u64 = (1 << 16) + 1;
    type FermatFqConfig = FqConfig<FERMAT_Q>;

    const FERMAT_NS: [usize; 13] = [
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
        // This now doesn't work because we initialise
        // Pow2CyclotomicPolyRing and it run from_coeffs_vec
        // which in turn casts a vector to [F; N], an array on stack,
        // causing an overflow.
        //
        // 1 << 14,
    ];
    macro_rules! generate_ntt_form_match {
        ($n:expr, $initial_coeffs:expr, $($size:expr),+) => {
            match $n {
                $(
                    $size => {
                        let ntt_form: Pow2CyclotomicPolyRingNTT<FERMAT_Q, $size> =
                                Pow2CyclotomicPolyRingNTT::<FERMAT_Q, $size>::from(Pow2CyclotomicPolyRing::from($initial_coeffs.clone()));
                        let intt_coeffs = Pow2CyclotomicPolyRing::<FERMAT_Q, $size>::from(ntt_form).coeffs();//
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

    fn test_mul_ntt_pow2<const SIZE: usize>()
    where
        Pow2CyclotomicPolyRingNTT<FERMAT_Q, SIZE>: From<Pow2CyclotomicPolyRing<FERMAT_Q, SIZE>>,
        Pow2CyclotomicPolyRing<FERMAT_Q, SIZE>: From<Pow2CyclotomicPolyRingNTT<FERMAT_Q, SIZE>>,
    {
        let mut rng = thread_rng();

        let coeff_1 = Pow2CyclotomicPolyRing::rand(&mut rng);
        let coeff_2 = Pow2CyclotomicPolyRing::<FERMAT_Q, { SIZE }>::rand(&mut rng);

        let ntt_form_1 = Pow2CyclotomicPolyRingNTT::from(coeff_1);
        let ntt_form_2 = Pow2CyclotomicPolyRingNTT::from(coeff_2);

        let ntt_mul = ntt_form_1 * ntt_form_2;
        let coeffs_mul = coeff_1 * coeff_2;
        // ntt_mul.coeffs() performs INTT while coeffs_mul.coeffs() just returns the coefficients
        assert_eq!(
            Pow2CyclotomicPolyRing::from(ntt_mul).coeffs(),
            coeffs_mul.coeffs()
        );
    }

    #[test]
    fn test_mul_ntt_pow2_hardcoded() {
        let coeffs_1_vec_vec = [1, 2, 3, 4, 5, 6, 7, 8];
        let coeffs_1_vec = coeffs_1_vec_vec
            .into_iter()
            .map(Fp::<MontBackend<FermatFqConfig, 1>, 1>::from)
            .collect::<Vec<_>>();
        let mut coeffs_1 = Pow2CyclotomicPolyRing::<FERMAT_Q, 8>::from(coeffs_1_vec);
        let coeffs_2_vec_vec = [8, 7, 6, 5, 4, 3, 2, 1];
        let coeffs_2_vec = coeffs_2_vec_vec
            .into_iter()
            .map(Fp::<MontBackend<FermatFqConfig, 1>, 1>::from)
            .collect::<Vec<_>>();
        let coeffs_2 = Pow2CyclotomicPolyRing::<FERMAT_Q, 8>::from(coeffs_2_vec);

        let mut ntt_form_1 = Pow2CyclotomicPolyRingNTT::from(coeffs_1);
        ntt_form_1 *= Pow2CyclotomicPolyRingNTT::from(coeffs_2);

        coeffs_1 *= coeffs_2;
        assert_eq!(
            Pow2CyclotomicPolyRing::<FERMAT_Q, 8>::from(ntt_form_1).coeffs(),
            coeffs_1.coeffs()
        );
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
