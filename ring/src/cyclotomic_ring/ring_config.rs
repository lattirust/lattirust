use ark_ff::{Field, Fp, FpConfig};
use ark_poly::EvaluationDomain;

use crate::Ring;

/// The trait for describing cyclotomic ring parameters.
/// It is used to specify:
/// * The field of coefficients of the cyclotomic ring.
/// * The field of CRT-components of the ring.
/// * Implementation of CRT/iCRT and reduction modulo the cyclotomic polynomial.
pub trait CyclotomicConfig<const N: usize>: Send + Sync + 'static + Sized {
    /// The base prime field configuration of the underlying polynomial ring.
    type BaseFieldConfig: FpConfig<N>;
    /// The field of the CRT components.
    type BaseCRTField: Field<BasePrimeField = Fp<Self::BaseFieldConfig, N>> + Ring;

    /// The field of the CRT components is a field extension of the base field.
    /// This constant specifies the degree.
    const CRT_FIELD_EXTENSION_DEGREE: usize;

    /// Given coefficients of a polynomial of degree 2 * phi(D)
    /// reduces it mod Phi_D(X).
    fn reduce_in_place(coefficients: &mut Vec<Fp<Self::BaseFieldConfig, N>>);

    /// Computes the evaluations of the polynomial with
    /// coefficients `coefficients` on the primitive Dth roots of unity.
    fn crt_in_place(coefficients: &mut [Fp<Self::BaseFieldConfig, N>]);

    fn crt(coefficients: Vec<Fp<Self::BaseFieldConfig, N>>) -> Vec<Self::BaseCRTField>;
    fn icrt(evaluations: Vec<Self::BaseCRTField>) -> Vec<Fp<Self::BaseFieldConfig, N>>;

    /// Given primitive Dth root of unity evaluations of a polynomial
    /// computes its coefficient form.
    fn icrt_in_place(evaluations: &mut [Fp<Self::BaseFieldConfig, N>]);
}

/// A configuration trait for fully splitting pow2 cyclotomic rings.
/// Non-efficient. Shouldn't be used for time-sensitive applications.
/// Piggybacks on the implementation of FFT in arkworks.
pub trait RpConfig<const N: usize>: Send + Sync + 'static + Sized {
    /// The FpConfig of the underlying Fp field.
    type FpConfig: FpConfig<N>;
    /// The type of evaluation domain to compute the NTT on.
    type EvaluationDomain: EvaluationDomain<Fp<Self::FpConfig, N>>;

    /// Constant instance of the evaluation domain.
    const DOMAIN: Self::EvaluationDomain;

    /// The parameter of the cyclotomic polynomial.
    /// The ring determined by the configuration is meant to be
    /// Fp[X]/(Phi_D(X)).
    const D: usize;
    /// Euler's totient of D: \varphi(D). I.e. number of components
    const PHI_D: usize;
    /// An element of Fp of order D.
    const PRIMITIVE_DTH_ROOT_OF_UNITY: Fp<Self::FpConfig, N>;

    /// Given an index i\in [0,..., D-1] tells if its coprime with D,
    /// which implies that \psi^i is again a primitive Dth root of unity
    /// if \psi is.
    fn is_primitive_index(i: usize) -> bool;

    /// Given coefficients of a polynomial of degree D
    /// reduces it mod Phi_D(X).
    fn reduce_in_place(coefficients: &mut Vec<Fp<Self::FpConfig, N>>);

    /// Computes the evaluations of the polynomial with
    /// coefficients `coefficients` on the primitive Dth roots of unity.
    fn crt_in_place(coefficients: &mut Vec<Fp<Self::FpConfig, N>>) {
        coefficients.resize(Self::PHI_D, <Fp<Self::FpConfig, N> as Field>::ZERO);

        // We resize the coefficient vector with D - PHI_D zeros
        // to have a polynomial of "degree" D.
        coefficients
            .extend((0..(Self::D - Self::PHI_D)).map(|_| <Fp<Self::FpConfig, N> as Field>::ZERO));
        Self::DOMAIN.fft_in_place(coefficients);

        // Once we've done the NTT we remove the evaluations at PRIMITIVE_INDICES.
        // Those evaluations correspond to the evaluations at non-primitive roots of unity.
        let evals: Vec<Fp<Self::FpConfig, N>> = coefficients
            .iter()
            .enumerate()
            .filter_map(|(i, &x)| Self::is_primitive_index(i).then_some(x))
            .collect();

        *coefficients = evals
    }

    /// Given primitive Dth root of unity evaluations of a polynomial
    /// computes its coefficient form.
    fn icrt_in_place(evaluations: &mut Vec<Fp<Self::FpConfig, N>>) {
        // Enlarge the evaluation vector to be of length D.
        // Add dummy zero evaluations at the non primitive indices.
        // let mut evals = evaluations.clone();
        let mut evals = Vec::new();
        evals.resize(Self::D, <Fp<Self::FpConfig, N> as Field>::ZERO);

        for (i, value) in evals.iter_mut().enumerate() {
            if Self::is_primitive_index(i) {
                *value = evaluations.remove(0);
            }
        }
        Self::DOMAIN.ifft_in_place(&mut evals);

        // Reduce the resulting polynomial mod Phi_D(X).
        Self::reduce_in_place(&mut evals);
        *evaluations = evals;
    }
}

impl<C: RpConfig<N>, const N: usize> CyclotomicConfig<N> for C {
    type BaseFieldConfig = C::FpConfig;
    type BaseCRTField = Fp<Self::BaseFieldConfig, N>;

    const CRT_FIELD_EXTENSION_DEGREE: usize = 1;

    #[inline(always)]
    fn reduce_in_place(coefficients: &mut Vec<Fp<Self::BaseFieldConfig, N>>) {
        <Self as RpConfig<N>>::reduce_in_place(coefficients);
    }

    #[inline(always)]
    fn crt_in_place(coefficients: &mut [Fp<Self::BaseFieldConfig, N>]) {
        let mut vec = coefficients.to_vec();

        <Self as RpConfig<N>>::crt_in_place(&mut vec);

        for (coeff_mut, vec_elem) in coefficients.iter_mut().zip(vec) {
            *coeff_mut = vec_elem;
        }
    }

    #[inline(always)]
    fn icrt_in_place(evaluations: &mut [Fp<Self::BaseFieldConfig, N>]) {
        let mut vec = evaluations.to_vec();

        <Self as RpConfig<N>>::icrt_in_place(&mut vec);

        for (eval_mut, vec_elem) in evaluations.iter_mut().zip(vec) {
            *eval_mut = vec_elem;
        }
    }

    fn crt(mut coefficients: Vec<Fp<Self::BaseFieldConfig, N>>) -> Vec<Self::BaseCRTField> {
        Self::crt_in_place(&mut coefficients);

        coefficients
    }

    fn icrt(mut evaluations: Vec<Self::BaseCRTField>) -> Vec<Fp<Self::BaseFieldConfig, N>> {
        Self::icrt_in_place(&mut evaluations);

        evaluations
    }
}
