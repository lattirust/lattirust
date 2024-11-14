use ark_std::{ptr, slice};

/// Enables Chinese Remainder Transform (CRT) conversion from coefficient representation to NTT representation.
pub trait CRT: Sized {
    type CRTForm;

    fn crt(self) -> Self::CRTForm;
    fn elementwise_crt(vec: Vec<Self>) -> Vec<Self::CRTForm> {
        let (ptr, len, cap) = vec.into_raw_parts();

        let slice = unsafe { slice::from_raw_parts_mut(ptr, len) };

        for (i, elem) in slice.iter_mut().enumerate() {
            unsafe {
                ptr::write(
                    &mut (*(ptr as *mut Self::CRTForm).add(i)),
                    ptr::read(elem).crt(),
                );
            }
        }

        unsafe { Vec::from_raw_parts(ptr as *mut Self::CRTForm, len, cap) }
    }
}

/// Enables Inverse Chinese Remainder Transform (ICRT) conversion from NTT back to coefficient representation.
pub trait ICRT: Sized {
    type ICRTForm;

    fn icrt(self) -> Self::ICRTForm;
    fn elementwise_icrt(vec: Vec<Self>) -> Vec<Self::ICRTForm> {
        let (ptr, len, cap) = vec.into_raw_parts();

        let slice = unsafe { slice::from_raw_parts_mut(ptr, len) };

        for (i, elem) in slice.iter_mut().enumerate() {
            unsafe {
                ptr::write(
                    &mut (*(ptr as *mut Self::ICRTForm).add(i)),
                    ptr::read(elem).icrt(),
                );
            }
        }

        unsafe { Vec::from_raw_parts(ptr as *mut Self::ICRTForm, len, cap) }
    }
}

#[macro_export]
macro_rules! impl_crt_icrt_for_a_ring {
    ($crt: ty, $icrt: ty, $ring_config: ty) => {
        impl CRT for $icrt {
            type CRTForm = $crt;

            fn crt(mut self) -> $crt {
                <$ring_config>::crt_in_place(&mut self.0);

                unsafe { ark_std::mem::transmute(self) }
            }
        }

        impl ICRT for $crt {
            type ICRTForm = $icrt;

            fn icrt(self) -> $icrt {
                let mut icrt_self: $icrt = unsafe { ark_std::mem::transmute(self) };

                <$ring_config>::icrt_in_place(&mut icrt_self.0);

                icrt_self
            }
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic_ring::models::{babybear, frog_ring, goldilocks, stark_prime};
    use ark_std::{test_rng, UniformRand};

    macro_rules! test_crt_conversion {
        ($poly_type:ty, $ntt_type:ty, $test_name:ident, $test_vec_name:ident) => {
            #[test]
            fn $test_name() {
                let mut rng = test_rng();
                for _ in 0..1000 {
                    let original_poly = <$poly_type>::rand(&mut rng);

                    // Convert to NTT form using CRT trait
                    let ntt_form: $ntt_type = original_poly.crt();

                    // Convert back to PolyRing using ICRT trait
                    let converted_back_poly: $poly_type = ntt_form.icrt();

                    // Check if the original and converted back polynomials are the same
                    assert_eq!(original_poly.0, converted_back_poly.0);
                }
            }

            #[test]
            fn $test_vec_name() {
                let mut rng = test_rng();
                for _ in 0..100 {
                    let original_vec: Vec<$poly_type> =
                        (0..100).map(|_| <$poly_type>::rand(&mut rng)).collect();

                    // Convert to NTT form using CRT trait
                    let ntt_form: Vec<$ntt_type> = CRT::elementwise_crt(original_vec.clone());

                    // Convert back to PolyRing using ICRT trait
                    let converted_back_vec: Vec<$poly_type> = ICRT::elementwise_icrt(ntt_form);

                    // Check if the original and converted back poly vectors are the same
                    assert_eq!(original_vec, converted_back_vec);
                }
            }
        };
    }

    test_crt_conversion!(
        babybear::RqPoly,
        babybear::RqNTT,
        test_crt_trait_conversion_babybear,
        test_crt_trait_vec_conversion_babybear
    );
    test_crt_conversion!(
        frog_ring::RqPoly,
        frog_ring::RqNTT,
        test_crt_trait_conversion_frog_ring,
        test_crt_trait_vec_conversion_frog_ring
    );
    test_crt_conversion!(
        goldilocks::RqPoly,
        goldilocks::RqNTT,
        test_crt_trait_conversion_goldilocks,
        test_crt_trait_vec_conversion_goldilocks
    );
    test_crt_conversion!(
        stark_prime::RqPoly,
        stark_prime::RqNTT,
        test_crt_trait_conversion_stark_prime,
        test_crt_trait_vec_conversion_stark_prime
    );
}
