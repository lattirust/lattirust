use super::{Fq, Fq3};
use ark_std::vec::*;

pub(super) fn fq_vec_to_fq3_vec(mut vec: Vec<Fq>) -> Vec<Fq3> {
    vec.shrink_to_fit();

    let (ptr, len, cap) = vec.into_raw_parts();

    assert_eq!(len, cap);

    unsafe { Vec::from_raw_parts(ptr as *mut Fq3, super::ntt::N, super::ntt::N) }
}

pub(super) fn fq3_vec_to_fq_vec(mut vec: Vec<Fq3>) -> Vec<Fq> {
    vec.shrink_to_fit();

    let (ptr, len, cap) = vec.into_raw_parts();

    assert_eq!(len, cap);

    unsafe { Vec::from_raw_parts(ptr as *mut Fq, super::ntt::D, super::ntt::D) }
}
