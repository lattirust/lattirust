use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::relations::ajtai_cm;

pub struct Crs<R: PolyRing> {
    pub crs_cm: ajtai_cm::Crs<R>,
}

impl<R: PolyRing> Crs<R> {
    pub fn new(crs_cm: ajtai_cm::Crs<R>) -> Crs<R> {
        Crs { crs_cm }
    }
}

pub struct Instance<R: PolyRing> {
    r: Vector<R>,
    v: R,
    cm: ajtai_cm::Instance<R>,
}

pub struct Witness<R: PolyRing> {
    f: ajtai_cm::Witness<R>,
}

pub fn is_satisfied<R: PolyRing>(crs: &Crs<R>, x: &Instance<R>, w: &Witness<R>) -> bool {
    ajtai_cm::is_satisfied(&crs.crs_cm, &x.cm, &w.f)
    // TODO: && mle[f^](r) = v
}