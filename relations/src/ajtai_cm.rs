use ark_std::rand;
use lattirust_arithmetic::linear_algebra::Matrix;
use lattirust_arithmetic::linear_algebra::Vector;
use lattirust_arithmetic::ring::PolyRing;

pub struct Crs<R: PolyRing> {
    pub n: usize,
    pub m: usize,
    pub norm_bound: usize,
    // Use usize, since we use the L-infinity norm
    pub ck: Matrix<R>,
}

impl<R: PolyRing> Crs<R> {
    pub fn new<Rng: rand::Rng + ?Sized>(
        n: usize,
        m: usize,
        norm_bound: usize,
        rng: &mut Rng,
    ) -> Crs<R> {
        let ck = Matrix::<R>::rand(n, m, rng);
        Crs {
            n,
            m,
            norm_bound,
            ck,
        }
    }
}

pub type Instance<R> = Vector<R>;

pub type Witness<R> = Vector<R>;

pub fn is_satisfied<R: PolyRing>(crs: &Crs<R>, x: &Instance<R>, w: &Witness<R>) -> bool {
    *x == &crs.ck * w && w.iter().all(|a| a.l2_norm() <= crs.norm_bound as f64)
}
