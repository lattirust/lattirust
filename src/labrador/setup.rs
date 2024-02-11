use crate::lattice_arithmetic::poly_ring::PolyRing;
pub use crate::labrador::common_reference_string::CommonReferenceString;

pub fn setup<R: PolyRing>(r: usize, n: usize, beta_sq: f64, num_constraints: usize, num_constant_constraints: usize) -> CommonReferenceString<R> {
    CommonReferenceString::<R>::new(r, n, beta_sq, num_constraints, num_constant_constraints)
}