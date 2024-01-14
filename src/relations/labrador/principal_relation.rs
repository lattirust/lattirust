use rand::thread_rng;

use crate::lattice_arithmetic::matrix::{Matrix, sample_uniform_mat_symmetric, sample_uniform_vec, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;

#[derive(Clone)]
#[allow(non_snake_case)]
pub struct QuadDotProdFunction<R: PolyRing> {
    pub A: Matrix<R>,
    pub phi: Vec<Vector<R>>,
    pub b: R,
}

#[allow(non_snake_case)]
impl<R: PolyRing> QuadDotProdFunction<R> {
    pub fn new(A: Matrix<R>, phi: Vec<Vector<R>>, b: R) -> Self {
        let r = A.nrows();
        assert_eq!(r, A.ncols());
        for i in 0..r {
            for j in i..r {
                assert_eq!(A.index((i, j)), A.index((j, i)));
            }
        }

        assert_eq!(r, phi.len());
        let n = phi[0].len();
        for i in 1..r {
            assert_eq!(n, phi[i].len());
        }

        Self { A, phi, b }
    }

    pub fn new_dummy(r: usize, n: usize) -> Self {
        Self::new(sample_uniform_mat_symmetric(r, r), vec![sample_uniform_vec(n); r], R::rand(&mut thread_rng()))
    }
}

#[derive(Clone)]
pub struct PrincipalRelation<R: PolyRing> {
    pub r: usize,
    // multiplicity
    pub n: usize,
    // rank
    pub norm_bound: f64,
    // beta
    pub quad_dot_prod_funcs: Vec<QuadDotProdFunction<R>>,
    pub ct_quad_dot_prod_funcs: Vec<QuadDotProdFunction<R>>,
}

impl<R: PolyRing> PrincipalRelation<R> {
    pub fn new_dummy(r: usize, n: usize, norm_bound: f64, num_constraints: usize, num_ct_constraints: usize) -> PrincipalRelation<R> {
        Self {
            r,
            n,
            norm_bound,
            quad_dot_prod_funcs: vec![QuadDotProdFunction::new_dummy(r, n); num_constraints],
            ct_quad_dot_prod_funcs: vec![QuadDotProdFunction::new_dummy(r, n); num_ct_constraints],
        }
    }
}