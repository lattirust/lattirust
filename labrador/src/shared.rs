#![allow(non_snake_case)]

use lattirust_arithmetic::linear_algebra::{Matrix, SymmetricMatrix};
use lattirust_arithmetic::linear_algebra::Vector;
use lattirust_arithmetic::ring::PolyRing;
use relations::principal_relation::PrincipalRelation;

use crate::common_reference_string::CommonReferenceString;

/// A (subset of a) transcript of one execution of the core Labrador protocol, as needed to compute the next instance/crs/witness for recursion
pub struct BaseTranscript<'a, R: PolyRing> {
    pub instance: &'a PrincipalRelation<R>,
    pub crs: &'a CommonReferenceString<R>,
    pub(crate) u_1: Option<Vector<R>>,
    pub(crate) Pi: Option<Vec<Matrix<R>>>,
    pub(crate) p: Option<Vector<R::BaseRing>>,
    pub(crate) psi: Option<Vec<Vector<R::BaseRing>>>,
    pub(crate) omega: Option<Vec<Vector<R::BaseRing>>>,
    pub(crate) b__: Option<Vec<R>>,
    pub(crate) alpha: Option<Vector<R>>,
    pub(crate) beta: Option<Vector<R>>,
    pub(crate) u_2: Option<Vector<R>>,
    pub(crate) c: Option<Vec<R>>,
    pub(crate) z: Option<Vector<R>>,
    pub(crate) t: Option<Vec<Vector<R>>>,
    pub(crate) G: Option<SymmetricMatrix<R>>,
    pub(crate) H: Option<SymmetricMatrix<R>>,

    // Note: phi is not part of the transcript, but it is a linear combination of values in the transcript
    pub(crate) phi: Option<Vec<Vector<R>>>,
}

impl<'a, R: PolyRing> BaseTranscript<'a, R> {
    pub fn init(
        crs: &'a CommonReferenceString<R>,
        instance: &'a PrincipalRelation<R>,
    ) -> BaseTranscript<'a, R> {
        BaseTranscript {
            instance: &instance,
            crs: &crs,
            u_1: None,
            Pi: None,
            p: None,
            psi: None,
            omega: None,
            b__: None,
            alpha: None,
            beta: None,
            u_2: None,
            c: None,
            z: None,
            t: None,
            G: None,
            H: None,
            phi: None,
        }
    }

    pub fn new_core(
        crs: &'a CommonReferenceString<R>,
        instance: &'a PrincipalRelation<R>,
        u_1: Vector<R>,
        Pi: Vec<Matrix<R>>,
        psi: Vec<Vector<R::BaseRing>>,
        omega: Vec<Vector<R::BaseRing>>,
        b__: Vec<R>,
        alpha: Vector<R>,
        beta: Vector<R>,
        u_2: Vector<R>,
        c: Vec<R>,
        phi: Vec<Vector<R>>,
    ) -> BaseTranscript<'a, R> {
        BaseTranscript {
            instance: &instance,
            crs: &crs,
            u_1: Some(u_1),
            Pi: Some(Pi),
            p: None,
            psi: Some(psi),
            omega: Some(omega),
            b__: Some(b__),
            alpha: Some(alpha),
            beta: Some(beta),
            u_2: Some(u_2),
            c: Some(c),
            z: None,
            t: None,
            G: None,
            H: None,
            phi: Some(phi),
        }
    }

    pub fn new(
        crs: &'a CommonReferenceString<R>,
        instance: &'a PrincipalRelation<R>,
        u_1: Vector<R>,
        Pi: Vec<Matrix<R>>,
        p: Vector<R::BaseRing>,
        psi: Vec<Vector<R::BaseRing>>,
        omega: Vec<Vector<R::BaseRing>>,
        b__: Vec<R>,
        alpha: Vector<R>,
        beta: Vector<R>,
        u_2: Vector<R>,
        c: Vec<R>,
        z: Vector<R>,
        t: Vec<Vector<R>>,
        G: SymmetricMatrix<R>,
        H: SymmetricMatrix<R>,
        phi: Vec<Vector<R>>,
    ) -> BaseTranscript<'a, R> {
        BaseTranscript {
            instance: &instance,
            crs: &crs,
            u_1: Some(u_1),
            Pi: Some(Pi),
            p: Some(p),
            psi: Some(psi),
            omega: Some(omega),
            b__: Some(b__),
            alpha: Some(alpha),
            beta: Some(beta),
            u_2: Some(u_2),
            c: Some(c),
            z: Some(z),
            t: Some(t),
            G: Some(G),
            H: Some(H),
            phi: Some(phi),
        }
    }
}
