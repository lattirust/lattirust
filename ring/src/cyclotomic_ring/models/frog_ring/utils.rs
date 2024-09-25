use super::{Fq, Fq4};

pub(super) fn fq_vec_to_fq4_vec(mut vec: Vec<Fq>) -> Vec<Fq4> {
    vec.shrink_to_fit();

    let (ptr, len, cap) = vec.into_raw_parts();

    assert_eq!(len, cap);

    unsafe { Vec::from_raw_parts(ptr as *mut Fq4, super::ntt::N, super::ntt::N) }
}

pub(super) fn fq4_vec_to_fq_vec(mut vec: Vec<Fq4>) -> Vec<Fq> {
    vec.shrink_to_fit();

    let (ptr, len, cap) = vec.into_raw_parts();

    assert_eq!(len, cap);

    unsafe { Vec::from_raw_parts(ptr as *mut Fq, super::ntt::D, super::ntt::D) }
}
