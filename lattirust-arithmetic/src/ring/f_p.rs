use std::convert::Into;

use ark_ff::{Fp, Fp64, FpConfig, MontBackend, MontConfig, PrimeField};
use num_bigint::{BigInt, BigUint};
use num_traits::Signed;

use crate::ring::ntt::{const_fq_from, generator, two_adic_root_of_unity};
use crate::ring::representatives::{SignedRepresentative, WithSignedRepresentative};
use crate::ring::{Zq, Zq1, ZqConfig};
use crate::traits::Modulus;

pub struct FqConfig<const Q: u64> {}

const fn to_bigint_assert_odd_prime<const Q: u64>() -> ark_ff::BigInt<1> {
    assert!(
        Q > 2 && const_primes::is_prime(Q),
        "You tried to instantiate an FqConfig<Q> with either Q = 2 or Q not a prime"
    );
    ark_ff::BigInt::<1>([Q])
}

impl<const Q: u64> MontConfig<1> for FqConfig<Q> {
    const MODULUS: ark_ff::BigInt<1> = to_bigint_assert_odd_prime::<Q>(); // Fails at compile time if Q is not an odd prime

    // TODO: As far as I can tell, `GENERATOR`and `TWO_ADIC_ROOT_OF_UNITY` are only used for FftField.
    const GENERATOR: Fp<MontBackend<Self, 1>, 1> = const_fq_from(generator::<Q>());
    const TWO_ADIC_ROOT_OF_UNITY: Fp<MontBackend<Self, 1>, 1> =
        const_fq_from(two_adic_root_of_unity::<Q>());
}

/// `Fq<Q>` is a prime field with modulus `Q`, where `Q` is less than 64 bits.
// TODO: The arkworks implementation is not correct for 64-bit primes Q. Not sure if this is due to some overflows, and if we should restrict Q to be e.g. <= 62 bits.
pub type Fq<const Q: u64> = Fp64<MontBackend<FqConfig<Q>, 1>>;

pub const fn fq_zero<const Q: u64>() -> Fq<Q> {
    MontBackend::<FqConfig<Q>, 1>::ZERO
}

impl<P: FpConfig<N>, const N: usize> Modulus for Fp<P, N> {
    fn modulus() -> BigUint {
        Self::MODULUS.into()
    }
}


pub trait FromZqSignedRepresentative<T> {
    fn from_zq_signed_representative(value: T) -> Self;
}

impl<C: ZqConfig<1>, P: FpConfig<N>, const N: usize> FromZqSignedRepresentative<SignedRepresentative<Zq<C, 1>>> for Fp<P, N> {
    fn from_zq_signed_representative(value: SignedRepresentative<Zq<C, 1>>) -> Self {
        let zq_modulus = Zq::<C, 1>::modulus();
        let fp_modulus = Self::modulus();
        
        if fp_modulus < zq_modulus {
            panic!("Cannot convert signed representative to {}, as the signed representative modulus {} is larger than the modulus {}.", 
                   zq_modulus, std::any::type_name::<Self>(), fp_modulus);
        }
        
        let mut bigint = value.0;
        if bigint.is_negative() {
            let modulus: BigInt = zq_modulus.into();
            bigint += modulus;
        }
        let biguint: BigUint = bigint.to_biguint().unwrap();
        Self::from(biguint)
    }
}

impl<P: FpConfig<N>, const N: usize> From<SignedRepresentative<Fp<P, N>>> for Fp<P, N> {
    fn from(value: SignedRepresentative<Fp<P, N>>) -> Self {
        let mut bigint = value.0;
        if bigint.is_negative() {
            let modulus: BigInt = Self::modulus().into();
            bigint += modulus;
        }
        let biguint: BigUint = bigint.to_biguint().unwrap();
        Self::from(biguint)
    }
}

impl<C: FpConfig<N>, const N: usize> WithSignedRepresentative for Fp<C, N> {
    type SignedRepresentative = SignedRepresentative<Fp<C, N>>;

    fn as_signed_representative(&self) -> Self::SignedRepresentative {
        (*self).into()
    }

    fn signed_representative_to_bigint(repr: &Self::SignedRepresentative) ->num_bigint::BigInt {
        repr.0.clone()  
    }
    
    fn signed_representative_from_bigint(value: num_bigint::BigInt) -> Option<Self::SignedRepresentative> {
        Some(SignedRepresentative::new(value))
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
