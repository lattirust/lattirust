use std::convert::Into;

use ark_ff::{BigInt as ArkBigInt, SmallFp, SmallFpConfig};
use num_bigint::{BigInt, BigUint};
use num_traits::Signed;

use crate::ring::ntt::{generator, largest_power_of_two_dividing, two_adic_root_of_unity};
use crate::ring::representatives::{SignedRepresentative, WithSignedRepresentative};
use crate::traits::Modulus;

pub struct FqConfig<const Q: u64> {}

const fn assert_odd_prime<const Q: u64>() {
    assert!(
        Q > 2 && const_primes::is_prime(Q),
        "You tried to instantiate an FqConfig<Q> with either Q = 2 or Q not a prime"
    );
}

const fn r_mod_q<const Q: u64>() -> u64 {
    ((1u128 << 64) % (Q as u128)) as u64
}

const fn mod_mul_u64<const Q: u64>(a: u64, b: u64) -> u64 {
    (((a as u128) * (b as u128)) % (Q as u128)) as u64
}

const fn montgomery_inv_u64<const Q: u64>() -> u64 {
    let n = Q as u128;
    let mut inv = 1u128;
    let mut i = 0;
    while i < 7 {
        inv = inv.wrapping_mul(2u128.wrapping_sub(n.wrapping_mul(inv)));
        i += 1;
    }
    inv.wrapping_neg() as u64
}

const fn montgomery_mul_u64<const Q: u64>(a: u64, b: u64) -> u64 {
    let t = (a as u128) * (b as u128);
    let m = (t as u64).wrapping_mul(montgomery_inv_u64::<Q>());
    let mq = (m as u128) * (Q as u128);
    let (.., carry) = (t as u64).overflowing_add(mq as u64);
    let mut reduced = (t >> 64) + (mq >> 64) + (carry as u128);
    if reduced >= Q as u128 {
        reduced -= Q as u128;
    }
    reduced as u64
}

impl<const Q: u64> FqConfig<Q> {
    const R_SQUARED_MOD_Q: u64 = mod_mul_u64::<Q>(r_mod_q::<Q>(), r_mod_q::<Q>());

    /// Convert a canonical residue to Montgomery form without wide division.
    /// Returns `None` for inputs outside `0..Q`.
    pub const fn from_canonical_u64(value: u64) -> Option<SmallFp<Self>> {
        if value >= Q {
            None
        } else {
            Some(SmallFp::from_raw(montgomery_mul_u64::<Q>(
                value,
                Self::R_SQUARED_MOD_Q,
            )))
        }
    }

    pub const fn from_u64(value: u64) -> SmallFp<Self> {
        let reduced = value % Q;
        SmallFp::from_raw(mod_mul_u64::<Q>(reduced, r_mod_q::<Q>()))
    }
}

impl<const Q: u64> SmallFpConfig for FqConfig<Q> {
    type T = u64;

    const MODULUS: Self::T = {
        assert_odd_prime::<Q>();
        Q
    };
    const MODULUS_U128: u128 = Q as u128;
    const GENERATOR: SmallFp<Self> = Self::from_u64(generator::<Q>());
    const ZERO: SmallFp<Self> = SmallFp::from_raw(0);
    const ONE: SmallFp<Self> = SmallFp::from_raw(r_mod_q::<Q>());
    const NEG_ONE: SmallFp<Self> = Self::from_u64(Q - 1);
    const TWO_ADICITY: u32 = largest_power_of_two_dividing(Q - 1);
    const TWO_ADIC_ROOT_OF_UNITY: SmallFp<Self> = Self::from_u64(two_adic_root_of_unity::<Q>());
    const SQRT_PRECOMP: Option<ark_ff::SqrtPrecomputation<SmallFp<Self>>> = None;

    #[inline(always)]
    fn add_assign(a: &mut SmallFp<Self>, b: &SmallFp<Self>) {
        let (mut val, overflow) = a.value.overflowing_add(b.value);
        if overflow {
            val = val.wrapping_add(u64::MAX - Q + 1);
        }
        if val >= Q {
            val -= Q;
        }
        a.value = val;
    }

    #[inline(always)]
    fn sub_assign(a: &mut SmallFp<Self>, b: &SmallFp<Self>) {
        if a.value >= b.value {
            a.value -= b.value;
        } else {
            a.value = Q - (b.value - a.value);
        }
    }

    #[inline(always)]
    fn double_in_place(a: &mut SmallFp<Self>) {
        let tmp = *a;
        Self::add_assign(a, &tmp);
    }

    #[inline(always)]
    fn neg_in_place(a: &mut SmallFp<Self>) {
        if a.value != 0 {
            a.value = Q - a.value;
        }
    }

    #[inline(always)]
    fn mul_assign(a: &mut SmallFp<Self>, b: &SmallFp<Self>) {
        a.value = montgomery_mul_u64::<Q>(a.value, b.value);
    }

    #[inline(always)]
    fn sum_of_products<const T: usize>(
        a: &[SmallFp<Self>; T],
        b: &[SmallFp<Self>; T],
    ) -> SmallFp<Self> {
        let mut acc = Self::ZERO;
        let mut i = 0;
        while i < T {
            let mut prod = a[i];
            Self::mul_assign(&mut prod, &b[i]);
            Self::add_assign(&mut acc, &prod);
            i += 1;
        }
        acc
    }

    #[inline(always)]
    fn square_in_place(a: &mut SmallFp<Self>) {
        let tmp = *a;
        Self::mul_assign(a, &tmp);
    }

    #[inline]
    fn inverse(a: &SmallFp<Self>) -> Option<SmallFp<Self>> {
        if a.value == 0 {
            return None;
        }

        let canonical = Self::into_bigint(*a).0[0] as i128;
        let mut t = 0i128;
        let mut new_t = 1i128;
        let mut r = Q as i128;
        let mut new_r = canonical;

        while new_r != 0 {
            let quotient = r / new_r;
            let tmp_t = t - quotient * new_t;
            t = new_t;
            new_t = tmp_t;
            let tmp_r = r - quotient * new_r;
            r = new_r;
            new_r = tmp_r;
        }

        if r != 1 {
            return None;
        }
        if t < 0 {
            t += Q as i128;
        }

        Some(Self::new(t as u64))
    }

    #[inline(always)]
    fn new(value: Self::T) -> SmallFp<Self> {
        let reduced = value % Q;
        SmallFp::from_raw(mod_mul_u64::<Q>(reduced, r_mod_q::<Q>()))
    }

    #[inline(always)]
    fn from_bigint(other: ArkBigInt<1>) -> Option<SmallFp<Self>> {
        let value = other.0[0];
        if value >= Q {
            None
        } else {
            Some(Self::new(value))
        }
    }

    #[inline(always)]
    fn into_bigint(other: SmallFp<Self>) -> ArkBigInt<1> {
        ArkBigInt::<1>([montgomery_mul_u64::<Q>(other.value, 1)])
    }
}

/// `Fq<Q>` is a prime field with modulus `Q`, where `Q` is less than 64 bits.
pub type Fq<const Q: u64> = SmallFp<FqConfig<Q>>;

pub const fn fq_zero<const Q: u64>() -> Fq<Q> {
    SmallFp::from_raw(0)
}

impl<P: SmallFpConfig> Modulus for SmallFp<P> {
    fn modulus() -> BigUint {
        P::MODULUS_U128.into()
    }
}

impl<M: Modulus, P: SmallFpConfig> From<SignedRepresentative<M>> for SmallFp<P> {
    fn from(value: SignedRepresentative<M>) -> Self {
        if Self::modulus() < M::modulus() {
            panic!("Cannot convert signed representative to {}, as the signed representative modulus {} is larger than the modulus {}.", M::modulus(), std::any::type_name::<Self>(),  Self::modulus());
        }
        let mut bigint = value.0;
        if bigint.is_negative() {
            let modulus: BigInt = Self::modulus().into();
            bigint += modulus;
        }
        let biguint: BigUint = bigint.to_biguint().unwrap();
        Self::from(biguint)
    }
}

impl<P: SmallFpConfig> WithSignedRepresentative for SmallFp<P> {
    type SignedRepresentative = SignedRepresentative<SmallFp<P>>;

    fn as_signed_representative(&self) -> Self::SignedRepresentative {
        (*self).into()
    }
}

#[cfg(test)]
mod test {
    use crate::*;

    use super::*;

    // Some primes, big and small. In particular, this tests that the implementation does not rely on any special structure of the prime, nor on the primes being specified in any particular order.
    const Q1: u64 = (1 << 31) - (1 << 27) + 1; // BabyBear prime, NTT-friendly up to N=2^26
    const Q2: u64 = 274177; // LaBRADOR modulus factor 1, NTT-friendly up to N=2^7
    const Q3: u64 = 67280421310721; // LaBRADOR modulus factor 2, NTT-friendly up to N=2^7
    const Q4: u64 = ((1u128 << 64) - (1u128 << 32) + 1) as u64; // Goldilocks prime, NTT-friendly up to N=2^31
    const Q5: u64 = 3; // Not NTT-friendly
    const Q6: u64 = (1 << 31) - 1; // Mersenne prime, not NTT-friendly
    const Q7: u64 = (1 << 13) - 1; // Not NTT-friendly
    const Q8: u64 = 27644437; // Not NTT-friendly
    const Q9: u64 = 200560490131; // Not NTT-friendly
    const Q10: u64 = 7; // Not NTT-friendly

    #[cfg(test)]
    mod test_f1 {
        use super::*;

        type F = Fq<Q1>;
        test_field_ring!(F, 100);
    }

    #[cfg(test)]
    mod test_f2 {
        use super::*;

        type F = Fq<Q2>;
        test_field_ring!(F, 100);
    }

    #[cfg(test)]
    mod test_f3 {
        use super::*;

        type F = Fq<Q3>;
        test_field_ring!(F, 100);
    }

    #[cfg(test)]
    mod test_f4 {
        use super::*;

        type F = Fq<Q4>;
        test_field_ring!(F, 100);
    }

    #[cfg(test)]
    mod test_f5 {
        use super::*;

        type F = Fq<Q5>;
        test_field_ring!(F, 100);
    }

    #[cfg(test)]
    mod test_f6 {
        use super::*;

        type F = Fq<Q6>;
        test_field_ring!(F, 100);
    }

    #[cfg(test)]
    mod test_f7 {
        use super::*;

        type F = Fq<Q7>;
        test_field_ring!(F, 100);
    }

    #[cfg(test)]
    mod test_f8 {
        use super::*;

        type F = Fq<Q8>;
        test_field_ring!(F, 100);
    }

    #[cfg(test)]
    mod test_f9 {
        use super::*;

        type F = Fq<Q9>;
        test_field_ring!(F, 100);
    }

    #[cfg(test)]
    mod test_f10 {
        use super::*;

        type F = Fq<Q10>;
        test_field_ring!(F, 100);
    }

    #[cfg(test)]
    mod test_signed_repr_1 {
        use super::*;

        type F = Fq<Q1>;
        test_signed_representative!(F, 100);
    }

    #[cfg(test)]
    mod test_signed_repr_2 {
        use super::*;

        type F = Fq<Q2>;
        test_signed_representative!(F, 100);
    }

    #[cfg(test)]
    mod test_signed_repr_3 {
        use super::*;

        type F = Fq<Q3>;
        test_signed_representative!(F, 100);
    }

    #[cfg(test)]
    mod test_signed_repr_4 {
        use super::*;

        type F = Fq<Q4>;
        test_signed_representative!(F, 100);
    }

    #[cfg(test)]
    mod test_signed_repr_5 {
        use super::*;

        type F = Fq<Q5>;
        test_signed_representative!(F, 100);
    }

    #[cfg(test)]
    mod test_signed_repr_6 {
        use super::*;

        type F = Fq<Q6>;
        test_signed_representative!(F, 100);
    }

    #[cfg(test)]
    mod test_signed_repr_7 {
        use super::*;

        type F = Fq<Q7>;
        test_signed_representative!(F, 100);
    }

    #[cfg(test)]
    mod test_signed_repr_8 {
        use super::*;

        type F = Fq<Q8>;
        test_signed_representative!(F, 100);
    }

    #[cfg(test)]
    mod test_signed_repr_9 {
        use super::*;

        type F = Fq<Q9>;
        test_signed_representative!(F, 100);
    }

    #[cfg(test)]
    mod test_signed_repr_10 {
        use super::*;

        type F = Fq<Q10>;
        test_signed_representative!(F, 100);
    }
}

#[cfg(test)]
mod canonical_conversion_tests {
    use super::*;

    fn check<const Q: u64>() {
        for value in [0, 1, Q / 2, Q - 1, Q, u64::MAX] {
            let expected = (value < Q).then(|| FqConfig::<Q>::from_u64(value));
            assert_eq!(FqConfig::<Q>::from_canonical_u64(value), expected);
        }
    }

    #[test]
    fn canonical_conversion_matches_division_based_conversion() {
        check::<2013265921>();
        check::<2113929217>();
        check::<1811939329>();
        check::<469762049>();
        check::<274877906837>();
        check::<18446744069414584321>();
    }
}
