use ark_ff::{Field, Fp3, Fp3Config, Fp64, MontBackend, MontFp};

mod ntt;

#[derive(ark_ff::MontConfig)]
#[modulus = "18446744069414584321"]
#[generator = "7"]
pub struct FqConfig;
pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

pub struct Goldilocks3Config;

impl Fp3Config for Goldilocks3Config {
    type Fp = Fq;

    const NONRESIDUE: Self::Fp = MontFp!("1099511627776");

    // Do we even need these?
    // NQR ^ (MODULUS^i - 1)/3, i=0,1,2 with NQR = u = (0,1,0)
    const FROBENIUS_COEFF_FP3_C1: &'static [Fq] = &[];

    // NQR ^ (2*MODULUS^i - 2)/3, i=0,1,2 with NQR = u = (0,1,0)
    const FROBENIUS_COEFF_FP3_C2: &'static [Fq] = &[];

    const TWO_ADICITY: u32 = 32;

    const TRACE_MINUS_ONE_DIV_TWO: &'static [u64] =
        &[9223372049739677694, 9223372049739677692, 2147483646];
    // 7 ^ t.
    const QUADRATIC_NONRESIDUE_TO_T: ark_ff::Fp3<Self> =
        Fp3::new(MontFp!("3607031617444012685"), Fq::ZERO, Fq::ZERO);
}

pub type Fq3 = Fp3<Goldilocks3Config>;

mod test {
    #[allow(unused_imports)]
    use crate::goldilocks::{Fq, Goldilocks3Config};
    #[allow(unused_imports)]
    use ark_ff::Fp3Config;
    #[allow(unused_imports)]
    use ark_std::One;

    // Note: if ord X = 24 then X can't be a cubic residue.
    #[test]
    fn test_nonresidue_is_order_24() {
        let nonresidue: Fq = Goldilocks3Config::NONRESIDUE;

        let mut pow = nonresidue;

        for _i in 0..22 {
            pow *= nonresidue;
            assert_ne!(pow, Fq::one());
        }

        assert_eq!(pow * nonresidue, Fq::one());
    }
}
