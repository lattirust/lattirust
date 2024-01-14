use ark_ff::{Field, Fp, MontBackend};
use ark_ff::fields::MontConfig;
use ark_relations::r1cs::ConstraintSystem;

use crate::labrador::setup::CommonReferenceString;
use crate::lattice_arithmetic::poly_ring::PolyRing;

#[derive(MontConfig)]
#[modulus = "18446744073709551617"] // 2^64+1
#[generator = "3"] // TODO
pub struct FqConfig;

pub type F64b = Fp<MontBackend<FqConfig, 2>, 2>;

pub fn prove_r1cs<R: PolyRing>(cs: ConstraintSystem<F64b>,crs: CommonReferenceString<R>) {
    assert_eq!(F64b::extension_degree(), 0);
    assert_eq!(F64b::characteristic(), [1, 1]);
    todo!("{:?} {:?}", cs, crs);
}