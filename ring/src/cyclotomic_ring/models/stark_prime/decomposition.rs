use ark_ff::PrimeField;
use num_bigint::{BigInt, BigUint, ToBigInt};
use num_traits::Euclid;

use crate::{ConvertibleRing, SignedRepresentative, UnsignedRepresentative};

use super::Fq;

impl ConvertibleRing for Fq {
    type UnsignedInt = UnsignedRepresentative<BigUint>;
    type SignedInt = SignedRepresentative<BigInt>;
}

impl From<SignedRepresentative<BigInt>> for Fq {
    fn from(value: SignedRepresentative<BigInt>) -> Self {
        let q: BigInt = UnsignedRepresentative::from(BigUint::from(Fq::MODULUS))
            .0
            .to_bigint()
            .unwrap();

        let r = value.0.rem_euclid(&q);

        Fq::from(r.to_biguint().unwrap())
    }
}

impl From<Fq> for SignedRepresentative<BigInt> {
    fn from(value: Fq) -> Self {
        let unsigned: BigUint = UnsignedRepresentative::from(value).0;
        let q_half: ark_ff::BigInt<4> = Fq::MODULUS_MINUS_ONE_DIV_TWO;
        let q: ark_ff::BigInt<4> = Fq::MODULUS;
        if unsigned > q_half.into() {
            SignedRepresentative(
                unsigned.to_bigint().unwrap() - BigUint::from(q).to_bigint().unwrap(),
            )
        } else {
            SignedRepresentative(unsigned.to_bigint().unwrap())
        }
    }
}

impl From<Fq> for UnsignedRepresentative<BigUint> {
    fn from(value: Fq) -> Self {
        UnsignedRepresentative::from(BigUint::from(value))
    }
}

impl From<UnsignedRepresentative<BigUint>> for Fq {
    fn from(value: UnsignedRepresentative<BigUint>) -> Self {
        Fq::from(ark_ff::BigInt::try_from(value.0).unwrap())
    }
}

#[cfg(test)]
mod test {
    use ark_ff::MontFp;

    use crate::{balanced_decomposition::Decompose, cyclotomic_ring::models::stark_prime::Fq};

    #[test]
    fn test_stark_prime_decomposition() {
        let x: Fq = MontFp!("253532532532352325");

        let decomposition = x.decompose(1 << 16, None);

        assert_eq!(
            decomposition,
            vec![
                -Fq::from(27323),
                -Fq::from(17255),
                -Fq::from(17793),
                Fq::from(901)
            ]
        )
    }
}
