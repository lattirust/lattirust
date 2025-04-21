use num_traits::Signed;

use crate::ring::representatives::WithSignedRepresentative;
use crate::ring::Ring;

pub fn powers_of_basis<R: Ring + WithSignedRepresentative>(basis: R, length: usize) -> Vec<R> {
    powers_of_basis_int(basis, length)
        .into_iter()
        .map(|x| Into::<R>::into(x))
        .collect()
}

fn powers_of_basis_int<R: Ring + WithSignedRepresentative>(basis: R, length: usize) -> Vec<R> {
    assert!(basis.as_signed_representative().is_positive());
    let mut pows = Vec::<R>::with_capacity(length);
    pows.push(R::one());
    for i in 1..length {
        pows.push(pows[i - 1] * basis);
    }
    pows
}
