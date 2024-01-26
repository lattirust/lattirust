use crate::lattice_arithmetic::poly_ring::PolyRing;
pub use crate::labrador::common_reference_string::CommonReferenceString;

pub fn setup<R: PolyRing>(r: usize, n: usize, d: usize, beta: f64, kappa: usize, kappa1: usize, kappa2: usize, num_constraints: usize, num_constant_constraints: usize, decomposition_basis: R::BaseRing) -> CommonReferenceString<R> {
    CommonReferenceString::<R>::new(r, n, d, beta, kappa, kappa1, kappa2, num_constraints, num_constant_constraints, decomposition_basis)
}