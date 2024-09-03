use derive_more::{Add, Div, From, Into, Mul, Neg, Rem, Sub};
use num_bigint::ToBigInt;
use num_traits::{Num, One, Signed, ToPrimitive, Zero};
use rounded_div::RoundedDiv;


#[derive(
    Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, From, Into, Add, Sub, Mul, Rem, Div, Neg,
)]
#[mul(forward)]
#[rem(forward)]
#[div(forward)]
pub struct ZqSignedRepresentative(pub(crate) num_bigint::BigInt);

impl Zero for ZqSignedRepresentative {
    fn zero() -> Self {
        ZqSignedRepresentative(num_bigint::BigInt::zero())
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl One for ZqSignedRepresentative {
    fn one() -> Self {
        ZqSignedRepresentative(num_bigint::BigInt::one())
    }
}

impl RoundedDiv<Self> for ZqSignedRepresentative {
    type Output = Self;

    fn rounded_div(self, other: Self) -> Self {
        f64::round(self.0.to_f64().unwrap() / other.0.to_f64().unwrap())
            .to_bigint()
            .unwrap()
            .into()
    }
}

impl Num for ZqSignedRepresentative {
    type FromStrRadixErr = <num_bigint::BigInt as Num>::FromStrRadixErr;

    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        num_bigint::BigInt::from_str_radix(str, radix).map(ZqSignedRepresentative)
    }
}

impl Signed for ZqSignedRepresentative {
    fn abs(&self) -> Self {
        self.0.abs().into()
    }

    fn abs_sub(&self, other: &Self) -> Self {
        self.0.abs_sub(&other.0).into()
    }

    fn signum(&self) -> Self {
        self.0.signum().into()
    }

    fn is_positive(&self) -> bool {
        self.0.is_positive()
    }

    fn is_negative(&self) -> bool {
        self.0.is_negative()
    }
}

impl From<i128> for ZqSignedRepresentative {
    fn from(value: i128) -> Self {
        value.to_bigint().unwrap().into()
    }
}
