use std::cmp::Ordering;
use std::convert::Into;
use std::fmt::{Debug, Display, Formatter};
use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{Num, One, Signed, ToPrimitive, Zero};
use rounded_div::RoundedDiv;

use crate::traits::Modulus;

/// A trait for types that can be represented as signed representatives, typically used for Z_q,
/// where representatives are in [-q/2, q/2-1].
///
/// If `Self::SignedRepresentative: Signed + Into<num_bigint::BigInt>`, then impls for
/// [WithL2Norm] and [WithLinfNorm] are automatically derived for `Self`.
pub trait WithSignedRepresentative: Sized + Clone {
    type SignedRepresentative: Signed
        + From<Self>
        + Into<Self>
        + TryFrom<i128, Error: Debug>
        + Clone
        + Debug
        + Send
        + Sync;

    fn as_signed_representative(&self) -> Self::SignedRepresentative {
        (*self).clone().into()
    }
}

/// For an odd prime [P], the signed representative of an element in Z_P is an integer in the range `[-floor(P/2), floor(P/2)]`.
#[derive(Default)]
pub struct SignedRepresentative<M: Modulus>(pub BigInt, std::marker::PhantomData<M>);

impl<M: Modulus> SignedRepresentative<M> {
    pub fn new(value: BigInt) -> Self {
        Self(value, std::marker::PhantomData)
    }

    #[inline]
    pub fn modulus() -> BigInt {
        M::modulus().into()
    }

    #[inline]
    pub fn min_inclusive() -> BigInt {
        -Self::max_inclusive()
    }

    #[inline]
    pub fn max_inclusive() -> BigInt {
        Self::modulus() / 2
    }
}

impl<M: Modulus> From<SignedRepresentative<M>> for BigInt {
    fn from(value: SignedRepresentative<M>) -> Self {
        value.0
    }
}

impl<M: Modulus> Debug for SignedRepresentative<M> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "FqSignedRepresentative({:?})", self.0)
    }
}

impl<M: Modulus> Display for SignedRepresentative<M> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<M: Modulus> Clone for SignedRepresentative<M> {
    fn clone(&self) -> Self {
        Self::new(self.0.clone())
    }
}

/// Maps [0, floor(O::modulus()/2)] to [0, floor(M::modulus()/2)] and [floor(O::modulus() / 2), O::modulus() - 1] to [-floor(M::modulus()/2), -1]
impl<M: Modulus, O: Modulus> From<O> for SignedRepresentative<M>
where
    BigUint: From<O>,
{
    fn from(value: O) -> Self {
        if M::modulus() < O::modulus() {
            panic!("Cannot convert {} to signed representative, as the modulus {} is larger than the signed representative modulus {}. ", std::any::type_name::<Self>(), O::modulus(), M::modulus());
        }
        let mut bigint = BigUint::from(value).into();
        let modulus: BigInt = O::modulus().into();
        let mod_half = modulus.clone().div(2);
        if bigint > mod_half {
            bigint -= modulus;
        }
        Self::new(bigint)
    }
}

impl<M: Modulus> Neg for SignedRepresentative<M> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.0)
    }
}

impl<M: Modulus> Add for SignedRepresentative<M> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut res = self.0.add(rhs.0);
        if res > Self::max_inclusive() {
            res -= Self::modulus();
        } else if res < Self::min_inclusive() {
            res += Self::modulus();
        }
        Self::new(res)
    }
}

impl<M: Modulus> Sub for SignedRepresentative<M> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut res = self.0.sub(rhs.0);
        if res > Self::max_inclusive() {
            res -= Self::modulus();
        } else if res < Self::min_inclusive() {
            res += Self::modulus();
        }
        Self::new(res)
    }
}

impl<M: Modulus> Mul for SignedRepresentative<M> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut res = self.0.mul(rhs.0);
        res %= &Self::modulus();

        if res.is_negative() {
            res += &Self::modulus();
        }

        if res > Self::max_inclusive() {
            res -= &Self::modulus();
        }
        Self::new(res)
    }
}

impl<M: Modulus> Div for SignedRepresentative<M> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self::new(self.0.div(rhs.0))
    }
}

impl<M: Modulus> Rem for SignedRepresentative<M> {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        Self::new(self.0.rem(rhs.0))
    }
}

impl<M: Modulus> RoundedDiv for SignedRepresentative<M> {
    type Output = Self;

    fn rounded_div(self, rhs: Self) -> Self::Output {
        let lhs = self.0;
        let rhs = rhs.0;

        if rhs.is_zero() {
            panic!("division by zero");
        }

        let (q, r) = lhs.div_rem(&rhs);

        let abs_rhs = rhs.abs();
        let double_r: BigInt = r.mul(2);

        let res = if double_r.abs() >= abs_rhs {
            if lhs.sign() == rhs.sign() {
                q + BigInt::one()
            } else {
                q - BigInt::one()
            }
        } else {
            q
        };
        Self::new(res)
    }
}

impl<M: Modulus> One for SignedRepresentative<M> {
    fn one() -> Self {
        Self::new(BigInt::one())
    }
}

impl<M: Modulus> Zero for SignedRepresentative<M> {
    fn zero() -> Self {
        Self::new(BigInt::zero())
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl<M: Modulus> PartialEq for SignedRepresentative<M> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<M: Modulus> PartialOrd for SignedRepresentative<M> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl<M: Modulus> Num for SignedRepresentative<M> {
    type FromStrRadixErr = ();

    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        let bigint = BigInt::from_str_radix(str, radix).map_err(|_| ())?;
        if bigint < Self::min_inclusive() || bigint > Self::max_inclusive() {
            return Err(());
        }
        Ok(Self::new(bigint))
    }
}

impl<M: Modulus> Signed for SignedRepresentative<M> {
    fn abs(&self) -> Self {
        Self::new(self.0.abs())
    }

    fn abs_sub(&self, other: &Self) -> Self {
        Self::new(self.0.abs_sub(&other.0))
    }

    fn signum(&self) -> Self {
        Self::new(self.0.signum())
    }

    fn is_positive(&self) -> bool {
        self.0.is_positive()
    }

    fn is_negative(&self) -> bool {
        self.0.is_negative()
    }
}

impl<M: Modulus> TryFrom<i128> for SignedRepresentative<M> {
    type Error = ();

    fn try_from(value: i128) -> Result<Self, Self::Error> {
        if value < Self::min_inclusive().to_i128().unwrap()
            || value > Self::max_inclusive().to_i128().unwrap()
        {
            return Err(());
        }
        Ok(Self::new(BigInt::from(value)))
    }
}

#[macro_export]
macro_rules! test_signed_representative {
    ($t:ty, $n:expr) => {

        fn test_vec() -> Vec<$t> {
            use ark_std::UniformRand;
            use num_traits::{One, Zero};

            let p = <$t>::modulus();
            let rng = &mut ark_std::test_rng();

            let mut vals: Vec<$t> = (0..$n).map(|_| <$t>::rand(rng)).collect();
            vals.extend(
                vec![
                    -<$t>::one(),
                    <$t>::zero(),
                    <$t>::one(),
                    <$t>::from(p.clone() / 2u64),
                    -<$t>::from(p.clone() / 2u64),
                    <$t>::from((p.clone() + 1u64) / 2u64),
                    -<$t>::from((p.clone() + 1u64) / 2u64),
                ]
                .into_iter(),
            );
            vals
        }

        #[test]
        fn test_signed_representative_from_into() {
            for value in test_vec() {
                let signed_repr: <$t as WithSignedRepresentative>::SignedRepresentative =
                    value.into();
                let value_from_repr: $t = signed_repr.clone().into();
                assert_eq!(value, value_from_repr);
            }
        }

        #[test]
        fn test_signed_representative_add() {
            use ark_std::UniformRand;
            let rng = &mut ark_std::test_rng();
            for a in test_vec() {
                let b = <$t>::rand(rng);
                let a_signed: <$t as WithSignedRepresentative>::SignedRepresentative = a.into();
                let b_signed: <$t as WithSignedRepresentative>::SignedRepresentative = b.into();

                let a_plus_b: $t = a + b;
                let a_plus_b_signed: <$t as WithSignedRepresentative>::SignedRepresentative =
                    a_signed + b_signed;
                let a_plus_b_from_repr: $t = a_plus_b_signed.clone().into();
                assert_eq!(a_plus_b, a_plus_b_from_repr);
            }
        }

        #[test]
        fn test_signed_representative_sub() {
            use ark_std::UniformRand;
            let rng = &mut ark_std::test_rng();
            for a in test_vec() {
                let b = <$t>::rand(rng);
                let a_signed: <$t as WithSignedRepresentative>::SignedRepresentative = a.into();
                let b_signed: <$t as WithSignedRepresentative>::SignedRepresentative = b.into();

                let a_minus_b: $t = a - b;
                let a_minus_b_signed: <$t as WithSignedRepresentative>::SignedRepresentative =
                    a_signed - b_signed;
                let a_minus_b_from_repr: $t = a_minus_b_signed.clone().into();
                assert_eq!(a_minus_b, a_minus_b_from_repr);
            }
        }

        #[test]
        fn test_signed_representative_mul() {
            use ark_std::UniformRand;
            let rng = &mut ark_std::test_rng();
            for a in test_vec() {
                let b = <$t>::rand(rng);
                let a_signed: <$t as WithSignedRepresentative>::SignedRepresentative = a.into();
                let b_signed: <$t as WithSignedRepresentative>::SignedRepresentative = b.into();

                let a_mul_b: $t = a * b;
                let a_mul_b_signed: <$t as WithSignedRepresentative>::SignedRepresentative =
                    a_signed * b_signed;
                let a_mul_b_from_repr: $t = a_mul_b_signed.clone().into();
                assert_eq!(a_mul_b, a_mul_b_from_repr);
            }
        }

        #[test]
        fn test_signed_representative_div_rem() {
            use num_traits::Zero;
            use ark_std::UniformRand;
            let rng = &mut ark_std::test_rng();
            for a in test_vec() {
                let mut b = <$t>::rand(rng);
                while b.is_zero() {
                    b = <$t>::rand(rng);
                }

                let a_signed: <$t as WithSignedRepresentative>::SignedRepresentative = a.into();
                let b_signed: <$t as WithSignedRepresentative>::SignedRepresentative = b.into();

                let div_result = a_signed.clone() / b_signed.clone();
                let rem_result = a_signed.clone() % b_signed.clone();

                let recomposed = b_signed.clone() * div_result.clone() + rem_result.clone();

                assert_eq!(recomposed, a_signed);
            }
        }

        #[test]
        fn test_signed_representative_rounded_div() {
            use num_traits::Zero;
            use rounded_div::RoundedDiv;
            use ark_std::UniformRand;
            let rng = &mut ark_std::test_rng();
            for a in test_vec() {
                let mut b = <$t>::rand(rng);
                while b.is_zero() {
                    b = <$t>::rand(rng);
                }

                let a_signed: <$t as WithSignedRepresentative>::SignedRepresentative = a.into();
                let b_signed: <$t as WithSignedRepresentative>::SignedRepresentative = b.into();

                let div_result = a_signed.clone().rounded_div(b_signed.clone());

                // Recompose approximately: (b * q) ≈ a
                let recomposed = b_signed.clone() * div_result.clone();

                // Since rounding may lose at most 0.5 * |b|, check that recomposed is "close enough" to a
                let diff = (a_signed.clone() - recomposed.clone()).0;

                // |diff| should be <= |b|/2
                let b_abs = b_signed.clone().0.abs();
                let b_half = b_abs.clone() / 2u32;

                assert!(
                    diff.clone().abs() <= b_half,
                    "RoundedDiv test failed: a = {:?}, b = {:?}, div = {:?}, recomposed = {:?}, diff = {:?}",
                    a_signed, b_signed, div_result, recomposed, diff
                );
            }
        }
    }
}
