use multiset::HashMultiSet;

use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::relations::traits::Relation;

pub struct Crs<R: PolyRing> {
    pub n_rows: usize,
    pub n_cols: usize,
    pub n_matrices: usize,
    pub n_private_inputs: usize,
    pub deg: usize,
    pub matrices: Vec<Matrix<R>>,
    pub multisets: Vec<HashMultiSet<usize>>,
    pub scalars: Vec<R>,
}

impl<R: PolyRing> Crs<R> {
    pub fn new(matrices: &[Matrix<R>], multisets: &Vec<HashMultiSet<usize>>, scalars: &[R], n_private_inputs: usize) -> Crs<R> {
        let n_matrices = matrices.len();
        let n_rows = matrices[0].nrows();
        let n_cols = matrices[0].ncols();
        let deg = multisets.iter().map(|m| m.len()).max().unwrap();
        assert_eq!(scalars.len(), n_matrices);
        assert!(matrices.iter().all(|m| m.nrows() == n_rows && m.ncols() == n_cols));
        assert!(n_private_inputs < n_cols);

        Crs {
            n_rows,
            n_cols,
            n_matrices,
            deg,
            n_private_inputs,
            matrices: matrices.to_vec(),
            multisets: multisets.clone(),
            scalars: scalars.to_vec(),
        }
    }
}

pub type Instance<R> = Vector<R>;

pub type Witness<R> = Vector<R>;

pub fn concat<R: PolyRing>(x: &Instance<R>, w: &Witness<R>) -> Vector<R> {
    Vector::<R>::from_vec([x.as_slice(), &[R::one()], w.as_slice()].concat())
}

pub struct Ccs<R: PolyRing> {
    _marker: std::marker::PhantomData<R>,
}

impl<R: PolyRing> Relation for Ccs<R> {
    type Instance = Instance<R>;
    type Witness = Witness<R>;
    type Crs = Crs<R>;

    fn is_satisfied(crs: &Crs<R>, x: &Instance<R>, w: &Witness<R>) -> bool {
        let z = &concat(x, w);
        assert_eq!(z.len(), crs.n_cols);
        let mut y: Vector<R> = Vector::zeros(crs.n_rows);
        for i in 0..crs.n_matrices {
            for j in crs.multisets[i].iter() {
                y += (&crs.matrices[*j] * z) * crs.scalars[i];
            }
        }
        y == Vector::zeros(crs.n_rows)
    }
}

pub fn is_satisfied<R: PolyRing>(crs: &Crs<R>, x: &Instance<R>, w: &Witness<R>) -> bool {
    Ccs::<R>::is_satisfied(crs, x, w)
}