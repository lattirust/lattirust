use ark_ff::{One, Zero};
use crate::{linear_algebra::Vector, ring::{PolyRing, Pow2CyclotomicPolyRing, Zq}, traits::FromRandomBytes};


pub struct PPKChallengeSet<R: PolyRing> {
    _marker: std::marker::PhantomData<R>,
}

impl<const M: u64, const N: usize> FromRandomBytes<Pow2CyclotomicPolyRing<Zq<M>, N>> 
    for PPKChallengeSet<Pow2CyclotomicPolyRing<Zq<M>, N>> {
        fn byte_size() -> usize {
            8
        }
    
        fn try_from_random_bytes(bytes: &[u8]) -> Option<Pow2CyclotomicPolyRing<Zq<M>, N>> {
            let mut coeffs: [Zq<M>; N] = [Zq::<M>::zero(); N];
            // TODO: fix, is it ok with the modulo?
            let index = usize::from_be_bytes(bytes.try_into().ok()?) % N;
            coeffs[index] = Zq::<M>::one();
            
            Some(Pow2CyclotomicPolyRing::from(coeffs))
        }
}

#[test]
fn generate_monomial() {
    let bytes = [0u8; 8];
    let result = PPKChallengeSet::<Pow2CyclotomicPolyRing<Zq<17>, 8>>::try_from_random_bytes(&bytes);
    assert!(result.is_some());

    let poly = result.unwrap();
    let coeffs = poly.coeffs().clone();
    assert_eq!(coeffs[0], Zq::<17>::one());

    for i in 1..8 {
        assert_eq!(coeffs[i], Zq::<17>::zero());
    }
}

#[test]
fn generate_monomial_with_non_zero_bytes() {
    let bytes = Vector::<u8>::rand(8, &mut rand::thread_rng());
    let bytes_array: [u8; 8] = bytes.clone().as_slice().try_into().expect("slice with incorrect length");
    let index = usize::from_be_bytes(bytes_array) % 8;
    
    let result = PPKChallengeSet::<Pow2CyclotomicPolyRing<Zq<17>, 8>>::try_from_random_bytes(bytes.as_slice());
    assert!(result.is_some());
    let poly = result.unwrap();
    let coeffs = poly.coeffs().clone();

    // println!("{index:?}");
    assert_eq!(coeffs[index], Zq::<17>::one());
    for i in 0..8 {
        if i != index {
            assert_eq!(coeffs[i], Zq::<17>::zero());
        }
    }
}
