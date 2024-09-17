use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
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
use num_traits::{One, Zero};

use crate::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
use crate::traits::FromRandomBytes;
use crate::z_q::FqConfig;
use crate::Ring;
use crate::{OverField, PolyRing};
use lattirust_linear_algebra::SVector;

use super::poly_ring::WithRot;

use ark_ff::{Fp, FpConfig, MontBackend};

/// A cyclotomic ring of the form Fp<C,N>/<X^D + 1>.
#[derive(From, Into)]
pub struct Pow2CyclotomicPolyRingNTTGeneral<C: FpConfig<N>, const N: usize, const D: usize>(
    SVector<Fp<C, N>, D>,
);

pub type Pow2CyclotomicPolyRingNTT<const Q: u64, const N: usize> =
    Pow2CyclotomicPolyRingNTTGeneral<MontBackend<FqConfig<Q>, 1>, 1, N>;

impl<C: FpConfig<N>, const N: usize, const D: usize> Pow2CyclotomicPolyRingNTTGeneral<C, N, D> {
    fn from_coefficients_vec(mut coeffs: Vec<Fp<C, N>>) -> Self {
        let eval_domain = Radix2EvaluationDomain::<Fp<C, N>>::new(2 * D);

        // We resize the coefficient vector with D zeros
        // to have a polynomial of "degree" 2D.
        coeffs.extend((0..D).map(|_| Fp::zero()));
        eval_domain.unwrap().fft_in_place(&mut coeffs);

        // Once we've done the NTT we remove the evaluations at even indices.
        // Those evaluations corresspond to the evaluations at non-primitive roots of unity.
        let coeffs: Vec<Fp<C, N>> = coeffs
            .into_iter()
            .enumerate()
            .filter_map(|(i, coeff)| if i % 2 != 0 { Some(coeff) } else { None })
            .collect();

        Self::from_array(coeffs.try_into().unwrap())
    }

    fn coefficients_vec(&self) -> Vec<Fp<C, N>> {
        let eval_domain = Radix2EvaluationDomain::<Fp<C, N>>::new(2 * D);

        // Enlarge the evaluation vector to be of length 2D.
        // Add dummy zero evaluations at the even indices.
        let mut evals: Vec<Fp<C, N>> = ((0..(2 * D)).map(|x| {
            if x % 2 != 0 {
                self.0[x / 2]
            } else {
                Fp::zero()
            }
        }))
        .collect();

        eval_domain.unwrap().ifft_in_place(&mut evals);

        // Reduce the resulting polynomial mod X^D + 1.
        let (left, right) = evals.split_at_mut(D);
        for (coeff_left, coeff_right) in left.iter_mut().zip(right) {
            *coeff_left -= coeff_right;
        }

        // Truncate the resulting vector to be of size D.
        evals.resize(D, Fp::zero());

        evals
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> PartialEq
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> Eq
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
}

impl<C: FpConfig<N>, const N: usize, const D: usize> Clone
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn clone(&self) -> Self {
        Self(self.0)
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> Copy
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
}

impl<C: FpConfig<N>, const N: usize, const D: usize> Debug
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> ark_std::fmt::Result {
        Debug::fmt(&self.0, f)
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> Display
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> ark_std::fmt::Result {
        Display::fmt(&self.0, f)
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> Hash
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn hash<H: ark_std::hash::Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> Pow2CyclotomicPolyRingNTTGeneral<C, N, D> {
    // const EVAL_DOMAIN: Option<Radix2EvaluationDomain<Fp<C, N>>> =
    //     Radix2EvaluationDomain::new(D);

    #[allow(dead_code)]
    pub(crate) type Inner = SVector<Fp<C, N>, D>;
    /// Constructs a polynomial from an array of coefficients in NTT form.
    /// This function is private since we can't enforce coeffs being in NTT form if called from outside.
    fn from_array(coeffs_ntt: [Fp<C, N>; D]) -> Self {
        Self(Self::Inner::const_from_array(coeffs_ntt))
    }

    /// Constructs a polynomial from a function specifying coefficients in non-NTT form.
    pub fn from_fn<F>(f: F) -> Self
    where
        F: FnMut(usize) -> Fp<C, N>,
    {
        let coeffs = Vec::from(core::array::from_fn::<_, D, _>(f));

        Self::from_coefficients_vec(coeffs)
    }
}

const fn vec_from_element<C: FpConfig<N>, const N: usize, const D: usize>(
    elem: Fp<C, N>,
) -> SVector<Fp<C, N>, D> {
    SVector::const_from_array([elem; D])
}

impl<C: FpConfig<N>, const N: usize, const D: usize> CanonicalSerialize
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
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

impl<C: FpConfig<N>, const N: usize, const D: usize> Valid
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn check(&self) -> Result<(), SerializationError> {
        self.0.check()
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> CanonicalDeserialize
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn deserialize_with_mode<R: Read>(
        reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        Self::Inner::deserialize_with_mode(reader, compress, validate).map(Self)
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> Ring
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    const ZERO: Self = Self(vec_from_element(<Fp<C, N> as Ring>::ZERO));
    const ONE: Self = Self(vec_from_element(<Fp<C, N> as Ring>::ONE));
}

impl<C: FpConfig<N>, const N: usize, const D: usize> FromRandomBytes<Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn byte_size() -> usize {
        D * Fp::<C, N>::byte_size()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        assert_eq!(bytes.len(), Self::byte_size());
        let coeffs = core::array::from_fn(|i| {
            Fp::<C, N>::try_from_random_bytes(
                &bytes[i * Fp::<C, N>::byte_size()..(i + 1) * Fp::<C, N>::byte_size()],
            )
            .unwrap()
        });
        Some(Self::from_array(coeffs))
    }
}

// impl<const Q: u64, const N: usize> Serialize for Pow2CyclotomicPolyRingNTT<Q, N> {
//     fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
//     where
//         S: Serializer,
//     {
//         ark_se(&self.coeffs(), serializer)
//     }
// }
//
// impl<'a, const Q: u64, const N: usize> Deserialize<'a> for Pow2CyclotomicPolyRingNTT<Q, N> {
//     fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
//     where
//         D: Deserializer<'a>,
//     {
//         ark_de(deserializer).map(|v: Vec<Fq<Q>>| Self::from(v))
//     }
// }

impl<C: FpConfig<N>, const N: usize, const D: usize> Default
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    #[inline(always)]
    fn default() -> Self {
        Self::zero()
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> Zero
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
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

impl<C: FpConfig<N>, const N: usize, const D: usize> One
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    #[inline(always)]
    fn one() -> Self {
        Self::ONE
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> Mul<Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.0.component_mul(&rhs.0).into()
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> Neg
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.0.neg().into()
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> UniformRand
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self::from_fn(|_| Fp::<C, N>::rand(rng))
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> MulAssign<Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn mul_assign(&mut self, rhs: Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> Add<&'a Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        self.0.add(rhs.0).into()
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> Sub<&'a Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output {
        self.0.sub(rhs.0).into()
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> Add<&'a mut Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn add(self, rhs: &'a mut Self) -> Self::Output {
        self.0.add(rhs.0).into()
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> Sub<&'a mut Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn sub(self, rhs: &'a mut Self) -> Self::Output {
        self.0.sub(rhs.0).into()
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> Mul<&'a mut Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn mul(self, rhs: &'a mut Self) -> Self::Output {
        self.0.component_mul(&rhs.0).into()
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> AddAssign<&'a mut Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn add_assign(&mut self, rhs: &'a mut Self) {
        self.0.add_assign(rhs.0)
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> SubAssign<&'a mut Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn sub_assign(&mut self, rhs: &'a mut Self) {
        self.0.sub_assign(rhs.0)
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> MulAssign<&'a mut Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn mul_assign(&mut self, rhs: &'a mut Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> From<Pow2CyclotomicPolyRing<Fp<C, N>, D>>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn from(value: Pow2CyclotomicPolyRing<Fp<C, N>, D>) -> Self {
        let coeffs: Vec<Fp<C, N>> = value.coeffs();

        Self::from_coefficients_vec(coeffs)
    }
}

// impl<const Q: u64, const N: usize> Into<Pow2CyclotomicPolyRing<Zq<Q>, N>>
//     for Pow2CyclotomicPolyRingNTT<Q, N>
// {
//     fn into(self) -> Pow2CyclotomicPolyRing<Zq<Q>, N> {
//         Pow2CyclotomicPolyRing::<Zq<Q>, N>::from(self.coeffs())
//     }
// }

impl<C: FpConfig<N>, const N: usize, const D: usize> Mul<Fp<C, N>>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn mul(self, rhs: Fp<C, N>) -> Self::Output {
        self.mul(Self::from_scalar(rhs))
    }
}

macro_rules! impl_from_primitive_type {
    ($primitive_type: ty) => {
        impl<C: FpConfig<N>, const N: usize, const D: usize> From<$primitive_type>
            for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
        {
            fn from(value: $primitive_type) -> Self {
                Self::from_scalar(Fp::<C, N>::from(value))
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

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> Mul<&'a Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        self.0.component_mul(&rhs.0).into()
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> AddAssign<&'a Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn add_assign(&mut self, rhs: &'a Self) {
        self.0.add_assign(rhs.0)
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> SubAssign<&'a Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn sub_assign(&mut self, rhs: &'a Self) {
        self.0.sub_assign(rhs.0)
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> MulAssign<&'a Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn mul_assign(&mut self, rhs: &'a Self) {
        self.0.component_mul_assign(&rhs.0);
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> Add<Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self.0.add(rhs.0).into()
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> Sub<Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self.0.sub(rhs.0).into()
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> AddAssign<Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn add_assign(&mut self, rhs: Self) {
        self.0.add_assign(rhs.0)
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> SubAssign<Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn sub_assign(&mut self, rhs: Self) {
        self.0.sub_assign(rhs.0)
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> Sum<Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> Sum<&'a Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> Product<Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<'a, C: FpConfig<N>, const N: usize, const D: usize> Product<&'a Self>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> PolyRing
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    type BaseRing = Fp<C, N>;
    fn coeffs(&self) -> Vec<Fp<C, N>> {
        self.coefficients_vec()
    }
    fn dimension() -> usize {
        D
    }

    fn from_scalar(v: Self::BaseRing) -> Self {
        // NTT([v, 0, ..., 0]) = ([v, ..., v])
        Self::from_array([v; D])
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> From<Vec<Fp<C, N>>>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn from(mut value: Vec<Fp<C, N>>) -> Self {
        value.resize_with(D, Fp::zero);

        Self::from_coefficients_vec(value)
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> From<Fp<C, N>>
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn from(value: Fp<C, N>) -> Self {
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

// impl<C: FpConfig<N>, const N: usize, const D: usize> WithL2Norm for Pow2CyclotomicPolyRingNTT<C, N, D> {
//     fn l2_norm_squared(&self) -> u128 {
//         self.coeffs().l2_norm_squared()
//     }
// }

// impl<C: FpConfig<N>, const N: usize, const D: usize> WithLinfNorm for Pow2CyclotomicPolyRingNTT<C, N, D> {
//     fn linf_norm(&self) -> u128 {
//         self.coeffs().linf_norm()
//     }
// }

impl<C: FpConfig<N>, const N: usize, const D: usize> WithRot
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn multiply_by_xi(&self, i: usize) -> Vec<Self::BaseRing> {
        let mut xi = vec![Self::BaseRing::ZERO; D];

        if i < D {
            xi[i] = Self::BaseRing::ONE;
        } else {
            unimplemented!("No support for multiplying with a polynomial bigger than N");
        }

        let xi_poly = Self::from_coefficients_vec(xi);

        let result = (*self * xi_poly).0;
        result.iter().copied().collect::<Vec<_>>()
    }
}

impl<C: FpConfig<N>, const N: usize, const D: usize> OverField
    for Pow2CyclotomicPolyRingNTTGeneral<C, N, D>
{
}

#[cfg(test)]
mod tests {
    use ark_ff::{Fp, MontBackend};

    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::{
        z_q::FqConfig, PolyRing, Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT,
        Pow2CyclotomicPolyRingNTTGeneral,
    };

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
                        let ntt_form = Pow2CyclotomicPolyRingNTT::<FERMAT_Q, $size>::from_fn(|i| $initial_coeffs[i]);
                        let intt_coeffs = ntt_form.coeffs();
                        assert_eq!($initial_coeffs, intt_coeffs);
                    },
                )+
                _ => panic!("Unsupported N value: {}", $n),
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

    #[test]
    fn test_mul_ntt_pow2() {
        // TODO: Add a couple more different dimensions.
        let mut rng = thread_rng();

        let coeff_1 =
            Pow2CyclotomicPolyRing::<Fp<MontBackend<FermatFqConfig, 1>, 1>, 16>::rand(&mut rng);
        let coeff_2 =
            Pow2CyclotomicPolyRing::<Fp<MontBackend<FermatFqConfig, 1>, 1>, 16>::rand(&mut rng);

        let ntt_form_1 = Pow2CyclotomicPolyRingNTTGeneral::from(coeff_1);
        let ntt_form_2 = Pow2CyclotomicPolyRingNTTGeneral::from(coeff_2);

        let ntt_mul = ntt_form_1 * ntt_form_2;
        let coeffs_mul = coeff_1 * coeff_2;
        // ntt_mul.coeffs() performs INTT while coeffs_mul_coeffs() just returns the coefficients
        assert_eq!(ntt_mul.coeffs(), coeffs_mul.coeffs());
    }
}
