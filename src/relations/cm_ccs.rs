use num_traits::Pow;

use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::relations::{ajtai_cm, ccs};

pub struct Crs<R: PolyRing> {
    pub crs_ccs: ccs::Crs<R>,
    pub crs_cm: ajtai_cm::Crs<R>,
    pub gadget_matrix: Matrix<R>,
}

impl<R: PolyRing> Crs<R> {
    pub fn new(crs_ccs: ccs::Crs<R>, crs_cm: ajtai_cm::Crs<R>) -> Crs<R> {
        let l = crs_cm.m / crs_ccs.n_cols;
        let b = crs_cm.norm_bound;
        debug_assert!(b.pow(l as u32) >= (R::modulus() / 2) as usize);
        let id_nc = Matrix::<R>::identity(crs_ccs.n_cols, crs_ccs.n_cols);
        let b_pows = Vector::<R>::from(
            (0..crs_ccs.n_cols).map(|i| R::from(b as u128).pow([i as u64])).collect::<Vec<_>>()
        );
        let gadget_matrix = id_nc.kronecker(&b_pows); // n_cols * l
        Crs { crs_ccs, crs_cm, gadget_matrix }
    }
}

pub struct Instance<R: PolyRing> {
    cm: ajtai_cm::Instance<R>,
    pub x_ccs: ccs::Instance<R>,
}

pub struct Witness<R: PolyRing> {
    pub f: ajtai_cm::Witness<R>,
    pub w_ccs: ccs::Witness<R>,
}

pub fn is_satisfied<R: PolyRing>(crs: &Crs<R>, x: &Instance<R>, w: &Witness<R>) -> bool {
    ccs::is_satisfied(&crs.crs_ccs, &x.x_ccs, &w.w_ccs) &&
        ajtai_cm::is_satisfied(&crs.crs_cm, &x.cm, &w.f) &&
        w.w_ccs == &crs.gadget_matrix * &w.f
}