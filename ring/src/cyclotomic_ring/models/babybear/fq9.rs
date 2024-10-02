use ark_ff::{CubicExtConfig, CubicExtField, Fp3, Fp3Config, SqrtPrecomputation};
use ark_std::marker::PhantomData;

use super::{BabyBear3ExtConfig, Fq, Fq9};

pub trait Fp9Config: 'static + Send + Sync + Sized {
    type Fp3Config: Fp3Config;

    const NONRESIDUE: Fp3<Self::Fp3Config>;

    const SQRT_PRECOMP: Option<SqrtPrecomputation<Fp9<Self>>> = None;

    /// Coefficients for the Frobenius automorphism in Fp9.
    const FROBENIUS_COEFF_FP9_C1: &'static [<Self::Fp3Config as Fp3Config>::Fp];
    const FROBENIUS_COEFF_FP9_C2: &'static [<Self::Fp3Config as Fp3Config>::Fp];

    #[inline(always)]
    fn mul_fp3_by_nonresidue_in_place(fe: &mut Fp3<Self::Fp3Config>) -> &mut Fp3<Self::Fp3Config> {
        let old_c2 = fe.c2;
        fe.c2 = fe.c1;
        fe.c1 = fe.c0;
        fe.c0 = old_c2;
        <Self::Fp3Config as Fp3Config>::mul_fp_by_nonresidue_in_place(&mut fe.c0);
        fe
    }
}

pub type Fq3 = Fp3<BabyBear3ExtConfig>;

pub struct Fp9ConfigWrapper<P: Fp9Config>(PhantomData<P>);

impl<P: Fp9Config> CubicExtConfig for Fp9ConfigWrapper<P> {
    type BasePrimeField = <P::Fp3Config as Fp3Config>::Fp;
    type BaseField = Fp3<P::Fp3Config>;
    type FrobCoeff = Self::BasePrimeField;

    const SQRT_PRECOMP: Option<ark_ff::SqrtPrecomputation<CubicExtField<Self>>> = P::SQRT_PRECOMP;

    const DEGREE_OVER_BASE_PRIME_FIELD: usize = 9;

    const NONRESIDUE: Self::BaseField = P::NONRESIDUE;

    const FROBENIUS_COEFF_C1: &'static [Self::FrobCoeff] = P::FROBENIUS_COEFF_FP9_C1;
    const FROBENIUS_COEFF_C2: &'static [Self::FrobCoeff] = P::FROBENIUS_COEFF_FP9_C2;

    #[inline(always)]
    fn mul_base_field_by_nonresidue_in_place(fe: &mut Self::BaseField) -> &mut Self::BaseField {
        P::mul_fp3_by_nonresidue_in_place(fe)
    }

    fn mul_base_field_by_frob_coeff(
        _c1: &mut Self::BaseField,
        _c2: &mut Self::BaseField,
        _power: usize,
    ) {
        unimplemented!() // Is this needed?
    }
}

pub type Fp9<P> = CubicExtField<Fp9ConfigWrapper<P>>;

// utils
pub(super) fn fq_vec_to_fq9_vec(mut vec: Vec<Fq>) -> Vec<Fq9> {
    vec.shrink_to_fit();

    let (ptr, len, cap) = vec.into_raw_parts();

    assert_eq!(len, cap);

    unsafe { Vec::from_raw_parts(ptr as *mut Fq9, super::ntt::N, super::ntt::N) }
}

pub(super) fn fq9_vec_to_fq_vec(mut vec: Vec<Fq9>) -> Vec<Fq> {
    vec.shrink_to_fit();

    let (ptr, len, cap) = vec.into_raw_parts();

    assert_eq!(len, cap);

    unsafe { Vec::from_raw_parts(ptr as *mut Fq, super::ntt::D, super::ntt::D) }
}
