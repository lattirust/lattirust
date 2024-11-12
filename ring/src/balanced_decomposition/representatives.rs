use ark_ff::{One, Zero};
use ark_serialize::{SerializationError, Valid};
use ark_std::ops::{Add, BitXor, Div, DivAssign, Mul, Neg, Rem, Sub};
use derive_more::{Add, AddAssign, Display, Mul, MulAssign, Product, Sub, SubAssign, Sum};
use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{sign::Signed, Num};

// Work-around to allow us implementing From traits
#[derive(
    Clone,
    Copy,
    Debug,
    Display,
    PartialEq,
    Eq,
    Add,
    AddAssign,
    Sub,
    SubAssign,
    Mul,
    MulAssign,
    Sum,
    Product,
)]
#[mul(forward)]
#[mul_assign(forward)]
#[repr(transparent)]
pub struct UnsignedRepresentative<T>(pub T);

#[derive(
    Clone,
    Copy,
    Debug,
    Display,
    PartialEq,
    Eq,
    Add,
    AddAssign,
    Sub,
    SubAssign,
    Mul,
    MulAssign,
    Sum,
    Product,
)]
#[mul(forward)]
#[mul_assign(forward)]
#[repr(transparent)]
pub struct SignedRepresentative<T>(pub T);

impl<T: Zero> Zero for SignedRepresentative<T> {
    fn zero() -> Self {
        SignedRepresentative(T::zero())
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl<T: One> One for SignedRepresentative<T> {
    fn one() -> Self {
        SignedRepresentative(T::one())
    }
}

impl<T: Signed> Signed for SignedRepresentative<T> {
    fn abs(&self) -> Self {
        SignedRepresentative(T::abs(&self.0))
    }

    fn abs_sub(&self, other: &Self) -> Self {
        SignedRepresentative(T::abs_sub(&self.0, &other.0))
    }

    fn signum(&self) -> Self {
        SignedRepresentative(T::signum(&self.0))
    }

    fn is_positive(&self) -> bool {
        T::is_positive(&self.0)
    }

    fn is_negative(&self) -> bool {
        T::is_negative(&self.0)
    }
}

impl<T: Neg<Output = T>> Neg for SignedRepresentative<T> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        SignedRepresentative(self.0.neg())
    }
}

impl<T: Zero> Zero for UnsignedRepresentative<T> {
    fn zero() -> Self {
        UnsignedRepresentative(T::zero())
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl<T: One> One for UnsignedRepresentative<T> {
    fn one() -> Self {
        UnsignedRepresentative(T::one())
    }
}

impl<'a, T: Clone + Mul<&'a T, Output = T>> Mul<&'a SignedRepresentative<T>>
    for &'a SignedRepresentative<T>
{
    type Output = SignedRepresentative<T>;

    fn mul(self, rhs: &'a SignedRepresentative<T>) -> Self::Output {
        SignedRepresentative(self.0.clone() * &rhs.0)
    }
}

impl<'a, T: Clone + Mul<&'a T, Output = T>> Mul<&'a UnsignedRepresentative<T>>
    for &'a UnsignedRepresentative<T>
{
    type Output = UnsignedRepresentative<T>;

    fn mul(self, rhs: &'a UnsignedRepresentative<T>) -> Self::Output {
        UnsignedRepresentative(self.0.clone() * &rhs.0)
    }
}

impl<T: Send + Sync> Valid for SignedRepresentative<T> {
    fn check(&self) -> Result<(), SerializationError> {
        Ok(())
    }
}

impl<T: Into<BigInt>> From<SignedRepresentative<T>> for BigInt {
    fn from(value: SignedRepresentative<T>) -> Self {
        value.0.into()
    }
}

impl<T: Into<BigUint>> From<UnsignedRepresentative<T>> for BigUint {
    fn from(value: UnsignedRepresentative<T>) -> Self {
        value.0.into()
    }
}

impl<T: PartialOrd> PartialOrd for SignedRepresentative<T> {
    fn partial_cmp(&self, other: &Self) -> Option<ark_std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl<T: PartialOrd> PartialOrd for UnsignedRepresentative<T> {
    fn partial_cmp(&self, other: &Self) -> Option<ark_std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl<T: Ord> Ord for SignedRepresentative<T> {
    fn cmp(&self, other: &Self) -> ark_std::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}
impl<T: Ord> Ord for UnsignedRepresentative<T> {
    fn cmp(&self, other: &Self) -> ark_std::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

impl<T: Num> Num for SignedRepresentative<T> {
    type FromStrRadixErr = T::FromStrRadixErr;

    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        T::from_str_radix(str, radix).map(|x| Self(x))
    }
}

impl<T: Num> Num for UnsignedRepresentative<T> {
    type FromStrRadixErr = T::FromStrRadixErr;

    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        T::from_str_radix(str, radix).map(|x| Self(x))
    }
}

impl<T: Div<Output = T>> Div for SignedRepresentative<T> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self(self.0.div(rhs.0))
    }
}

impl<T: Rem<Output = T>> Rem for SignedRepresentative<T> {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        Self(self.0.rem(rhs.0))
    }
}

impl<T: Rem<i128, Output = T>> Rem<i128> for SignedRepresentative<T> {
    type Output = Self;

    fn rem(self, rhs: i128) -> Self::Output {
        Self(self.0.rem(rhs))
    }
}

impl<T: Div<Output = T>> Div for UnsignedRepresentative<T> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self(self.0.div(rhs.0))
    }
}

impl<T: Rem<Output = T>> Rem for UnsignedRepresentative<T> {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        Self(self.0.rem(rhs.0))
    }
}

impl<T: Integer> Integer for SignedRepresentative<T> {
    fn div_floor(&self, other: &Self) -> Self {
        Self(self.0.div_floor(&other.0))
    }

    fn mod_floor(&self, other: &Self) -> Self {
        Self(self.0.mod_floor(&other.0))
    }

    fn gcd(&self, other: &Self) -> Self {
        Self(self.0.gcd(&other.0))
    }

    fn lcm(&self, other: &Self) -> Self {
        Self(self.0.lcm(&other.0))
    }

    fn is_multiple_of(&self, other: &Self) -> bool {
        self.0.is_multiple_of(&other.0)
    }

    fn is_even(&self) -> bool {
        self.0.is_even()
    }

    fn is_odd(&self) -> bool {
        self.0.is_odd()
    }

    fn div_rem(&self, other: &Self) -> (Self, Self) {
        let (q, r) = self.0.div_rem(&other.0);

        (Self(q), Self(r))
    }
}

impl<T: Integer> Integer for UnsignedRepresentative<T> {
    fn div_floor(&self, other: &Self) -> Self {
        Self(self.0.div_floor(&other.0))
    }

    fn mod_floor(&self, other: &Self) -> Self {
        Self(self.0.mod_floor(&other.0))
    }

    fn gcd(&self, other: &Self) -> Self {
        Self(self.0.gcd(&other.0))
    }

    fn lcm(&self, other: &Self) -> Self {
        Self(self.0.lcm(&other.0))
    }

    fn is_multiple_of(&self, other: &Self) -> bool {
        self.0.is_multiple_of(&other.0)
    }

    fn is_even(&self) -> bool {
        self.0.is_even()
    }

    fn is_odd(&self) -> bool {
        self.0.is_odd()
    }

    fn div_rem(&self, other: &Self) -> (Self, Self) {
        let (q, r) = self.0.div_rem(&other.0);

        (Self(q), Self(r))
    }
}

impl<T: Add<i128, Output = T>> Add<i128> for SignedRepresentative<T> {
    type Output = Self;
    fn add(self, value: i128) -> Self {
        Self(self.0 + value)
    }
}

impl<T: Sub<i128, Output = T>> Sub<i128> for SignedRepresentative<T> {
    type Output = Self;
    fn sub(self, value: i128) -> Self {
        Self(self.0 - value)
    }
}

impl<T: Div<i128, Output = T>> Div<i128> for SignedRepresentative<T> {
    type Output = Self;
    fn div(self, value: i128) -> Self {
        Self(self.0 / value)
    }
}

impl<T: DivAssign> DivAssign for SignedRepresentative<T> {
    fn div_assign(&mut self, rhs: Self) {
        self.0.div_assign(rhs.0);
    }
}

impl<T: DivAssign<i128>> DivAssign<i128> for SignedRepresentative<T> {
    fn div_assign(&mut self, rhs: i128) {
        self.0.div_assign(rhs);
    }
}

impl<T: DivAssign> DivAssign for UnsignedRepresentative<T> {
    fn div_assign(&mut self, rhs: Self) {
        self.0.div_assign(rhs.0);
    }
}

impl<T> From<T> for SignedRepresentative<T> {
    fn from(value: T) -> Self {
        Self(value)
    }
}

impl<T> From<T> for UnsignedRepresentative<T> {
    fn from(value: T) -> Self {
        Self(value)
    }
}

impl<T: BitXor<Output = T>> BitXor<SignedRepresentative<T>> for SignedRepresentative<T> {
    type Output = Self;

    fn bitxor(self, rhs: SignedRepresentative<T>) -> Self::Output {
        Self(self.0.bitxor(rhs.0))
    }
}

impl<T: BitXor<Output = T>> BitXor<UnsignedRepresentative<T>> for UnsignedRepresentative<T> {
    type Output = Self;

    fn bitxor(self, rhs: UnsignedRepresentative<T>) -> Self::Output {
        Self(self.0.bitxor(rhs.0))
    }
}

impl From<u128> for UnsignedRepresentative<BigUint> {
    fn from(value: u128) -> Self {
        Self(BigUint::from(value))
    }
}

impl From<i128> for SignedRepresentative<BigInt> {
    fn from(value: i128) -> Self {
        Self(BigInt::from(value))
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;

    use super::*;
    use crate::zn::z_q::Zq;

    const Q: u128 = 2u128.pow(61) - 1;
    const Q_HALF: i128 = (Q as i128 - 1) / 2;
    const TEST_STEP_SIZE: usize = (Q / 10) as usize;

    type F = Zq<{ Q as u64 }>;

    #[test]
    fn test_unsigned_representative() {
        for i in (0..Q).step_by(TEST_STEP_SIZE) {
            let f1 = F::from(i);
            let v1: UnsignedRepresentative<u128> = UnsignedRepresentative::from(f1);
            assert_eq!(i, v1.0);
        }
    }

    #[test]
    fn test_signed_representative() {
        assert_eq!(Q_HALF, F::MODULUS_MINUS_ONE_DIV_TWO.0[0] as i128);
        for i in (-Q_HALF..=Q_HALF).step_by(TEST_STEP_SIZE) {
            let v1 = SignedRepresentative::from(i);
            let f2 = F::from(v1);
            let v2 = SignedRepresentative::from(f2);
            assert_eq!(v1.0, v2.0);
        }
    }
}
