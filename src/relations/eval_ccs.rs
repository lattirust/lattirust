use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::relations::{ajtai_cm, ccs, cm_ccs, eval};

pub struct Crs<R: PolyRing> {
    pub crs_cm_ccs: cm_ccs::Crs<R>,
}

impl<R: PolyRing> Crs<R> {
    pub fn new(crs_cm_ccs: cm_ccs::Crs<R>) -> Crs<R> {
        Crs { crs_cm_ccs }
    }
}

pub struct Instance<R: PolyRing> {
    x_cm_ccs: cm_ccs::Instance<R>,
    r: Vector<R>,
    v: R,
}

pub struct Witness<R: PolyRing> {
    w_cm_ccs: cm_ccs::Witness<R>,
}

pub fn is_satisfied<R: PolyRing>(crs: &Crs<R>, x: &Instance<R>, w: &Witness<R>) -> bool {
    cm_ccs::is_satisfied(&crs.crs_cm_ccs, &x.x_cm_ccs, &w.w_cm_ccs)
// TODO
}