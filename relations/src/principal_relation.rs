#![allow(non_snake_case)]
use ark_std::rand;
use ark_std::rand::thread_rng;
use derive_more::Display;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;

use lattirust_arithmetic::linear_algebra::inner_products::inner_products;
use lattirust_arithmetic::linear_algebra::Matrix;
use lattirust_arithmetic::linear_algebra::Vector;
use lattirust_arithmetic::ring::PolyRing;
use lattirust_arithmetic::traits::WithL2Norm;

use crate::Relation;

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Index<R: PolyRing> {
    /// Number of witness vectors
    pub r: usize,
    /// Number of entries in a witness vector
    pub n: usize,
    /// Square of the L2-norm bound on the concatenation of witness vectors
    pub norm_bound_squared: f64,
    /// Number of quadratic-linear constraints
    pub num_constraints: usize,
    /// Number of quadratic-linear constraints on constant coefficients
    pub num_constant_constraints: usize,
    _marker: std::marker::PhantomData<R>,
}

#[derive(Clone, Eq, PartialEq, Debug, Display)]
#[display(
    "PrincipalRelation::Instance: \nquad_dot_prod_funcs: {:?}\nct_quad_dot_prod_funcs: {:?}",
    quad_dot_prod_funcs,
    ct_quad_dot_prod_funcs
)]
pub struct Instance<R: PolyRing> {
    pub quad_dot_prod_funcs: Vec<QuadDotProdFunction<R>>,
    pub ct_quad_dot_prod_funcs: Vec<ConstantQuadDotProdFunction<R>>,
}

#[derive(Clone, Debug, Eq, PartialEq, Display)]
#[display("QuadDotProdConstraint: A: {:?}, phi: {:?}, b: {:?}", A, phi, b)]
pub struct QuadDotProdFunction<R: PolyRing> {
    // TODO: A is always symmetric, so we could at least use a symmetric matrix type. A is also very sparse in some cases.
    pub A: Option<Matrix<R>>,
    // TODO: phi can be quite sparse
    pub phi: Vec<Vector<R>>,
    pub b: R,
    _private: (), // Forbid direct initialization, force users to use new(), which does some basis debug_asserts
}

impl<R: PolyRing> QuadDotProdFunction<R> {
    pub fn new(A: Matrix<R>, phi: Vec<Vector<R>>, b: R) -> Self {
        let (r, n) = (A.nrows(), phi[0].len());
        debug_assert_eq!(A.ncols(), r, "A should be square");
        debug_assert_eq!(A.transpose(), A, "A should be symmetric");

        debug_assert_eq!(
            phi.len(),
            r,
            "phi should have the same length as the dimensions of A"
        );
        debug_assert!(
            phi.iter().all(|phi_i| phi_i.len() == n),
            "each phi_i should have the same length"
        );
        Self {
            A: Some(A),
            phi,
            b,
            _private: (),
        }
    }

    pub fn new_linear(phi: Vec<Vector<R>>, b: R) -> Self {
        let n = phi[0].len();
        debug_assert!(
            phi.iter().all(|phi_i| phi_i.len() == n),
            "each phi_i should have the same length"
        );
        Self {
            A: None,
            phi,
            b,
            _private: (),
        }
    }

    pub fn new_dummy<Rng: rand::Rng + ?Sized>(r: usize, n: usize, rng: &mut Rng) -> Self {
        Self::new(
            Matrix::<R>::rand_symmetric(r, rng),
            vec![Vector::<R>::rand(n, rng); r],
            R::zero(),
        )
    }

    pub fn new_empty(r: usize, n: usize) -> Self {
        Self::new(
            Matrix::<R>::zeros(r, r),
            vec![Vector::<R>::zeros(n); r],
            R::zero(),
        )
    }

    pub fn is_valid_witness(&self, witness: &Witness<R>) -> bool {
        let inner_prods = inner_products(&witness.s);

        let mut res = R::zero();
        if let Some(A) = &self.A {
            let r = A.nrows();
            for i in 0..r {
                for j in 0..r {
                    res += A[(i, j)] * inner_prods[(i, j)];
                }
            }
        }

        for i in 0..self.phi.len() {
            res += self.phi[i].dot(&witness.s[i]);
        }

        res == self.b
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Display)]
#[display(
    "ConstantQuadDotProdConstraint: A: {:?}, phi: {:?}, b: {:?}",
    A,
    phi,
    b
)]
pub struct ConstantQuadDotProdFunction<R: PolyRing> {
    pub A: Option<Matrix<R>>,
    pub phi: Vec<Vector<R>>,
    pub b: R::BaseRing,
    _private: (), // Forbid direct initialization, force users to use new(), which does some basis debug_asserts
}

impl<R: PolyRing> ConstantQuadDotProdFunction<R> {
    pub fn new(A: Matrix<R>, phi: Vec<Vector<R>>, b: R::BaseRing) -> Self {
        let (r, n) = (A.nrows(), phi[0].len());
        debug_assert_eq!(A.ncols(), r, "A should be square");
        debug_assert_eq!(A.transpose(), A, "A should be symmetric");

        debug_assert_eq!(
            phi.len(),
            r,
            "phi should have the same length as the dimensions of A"
        );
        debug_assert!(
            phi.iter().all(|phi_i| phi_i.len() == n),
            "each phi_i should have the same length"
        );
        Self {
            A: Some(A),
            phi,
            b,
            _private: (),
        }
    }

    pub fn new_linear(phi: Vec<Vector<R>>, b: R::BaseRing) -> Self {
        let n = phi[0].len();
        debug_assert!(
            phi.iter().all(|phi_i| phi_i.len() == n),
            "each phi_i should have the same length"
        );
        Self {
            A: None,
            phi,
            b,
            _private: (),
        }
    }

    pub fn new_dummy<Rng: rand::Rng + ?Sized>(r: usize, n: usize, rng: &mut Rng) -> Self {
        Self::new(
            Matrix::<R>::rand_symmetric(r, rng),
            vec![Vector::<R>::rand(n, rng); r],
            R::BaseRing::zero(),
        )
    }

    pub fn is_valid_witness(&self, witness: &Witness<R>) -> bool {
        let inner_prods = inner_products(&witness.s);

        let mut res = R::zero();
        if let Some(A) = &self.A {
            let r = A.nrows();
            for i in 0..r {
                for j in 0..i + 1 {
                    res += A[(i, j)] * inner_prods[(i, j)];
                }
                for j in i + 1..r {
                    res += A[(i, j)] * inner_prods[(j, i)];
                }
            }
        }

        for i in 0..self.phi.len() {
            res += self.phi[i].dot(&witness.s[i]);
        }

        res.coefficients()[0] == self.b
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Witness<R: PolyRing> {
    pub s: Vec<Vector<R>>,
}

impl<R: PolyRing> Witness<R> {
    pub fn zero(rank: usize, multiplicity: usize) -> Self {
        Self {
            s: vec![Vector::<R>::from_fn(multiplicity, |_, _| R::zero()); rank],
        }
    }
}

impl<R: PolyRing> Instance<R> {
    pub fn is_valid_witness(&self, witness: &Witness<R>) -> anyhow::Result<()> {
        let unsat_indices = self
            .quad_dot_prod_funcs
            .iter()
            .enumerate()
            .filter_map(|(idx, c)| (!c.is_valid_witness(witness)).then(|| idx))
            .collect::<Vec<_>>();
        let ct_unsat_indices = self
            .ct_quad_dot_prod_funcs
            .iter()
            .enumerate()
            .filter_map(|(idx, c)| (!c.is_valid_witness(witness)).then(|| idx))
            .collect::<Vec<_>>();

        if unsat_indices.len() > 0 && ct_unsat_indices.len() > 0 {
            return Err(anyhow::anyhow!(
                "Standard constraints {:?} are not satisfied",
                unsat_indices
            ));
        } else if unsat_indices.len() > 0 {
            return Err(anyhow::anyhow!(
                "Standard constraints {:?} are not satisfied",
                unsat_indices
            ));
        } else if ct_unsat_indices.len() > 0 {
            return Err(anyhow::anyhow!(
                "Constant constraints {:?} are not satisfied",
                ct_unsat_indices
            ));
        }
        Ok(())
    }
}

pub struct Size {
    pub num_witnesses: usize,
    pub witness_size: usize,
    pub norm_bound_sq: f64,
    pub num_constraints: usize,
    pub num_constant_constraints: usize,
}

impl<R: PolyRing> Index<R> {
    pub fn new(size: &Size) -> Self {
        Self {
            r: size.num_witnesses,
            n: size.witness_size,
            norm_bound_squared: size.norm_bound_sq,
            num_constraints: size.num_constraints,
            num_constant_constraints: size.num_constant_constraints,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn is_wellformed_constraint(&self, c: &QuadDotProdFunction<R>) -> bool {
        match c.A {
            Some(ref A) => A.nrows() == self.r && A.ncols() == self.r && A.transpose() == *A,
            None => true,
        }
    }

    pub fn is_wellformed_const_constraint(&self, c: &ConstantQuadDotProdFunction<R>) -> bool {
        match c.A {
            Some(ref A) => A.nrows() == self.r && A.ncols() == self.r && A.transpose() == *A,
            None => true,
        }
    }

    pub fn is_wellformed_instance(&self, instance: &Instance<R>) -> anyhow::Result<()> {
        if instance.quad_dot_prod_funcs.len() != self.num_constraints {
            return Err(anyhow::anyhow!(
                "Number of quadratic-linear constraints does not match num_constraints"
            ));
        }
        if instance.ct_quad_dot_prod_funcs.len() != self.num_constant_constraints {
            return Err(anyhow::anyhow!(
                "Number of constant quadratic-linear constraints does not match num_constant_constraints"
            ));
        }
        let malformed_constraint_indices: Vec<_> = instance
            .quad_dot_prod_funcs
            .iter()
            .enumerate()
            .filter_map(|(idx, c)| (!self.is_wellformed_constraint(c)).then(|| idx))
            .collect();
        if malformed_constraint_indices.len() > 0 {
            return Err(anyhow::anyhow!(
                "Constraints {:?} are not well-formed",
                malformed_constraint_indices
            ));
        }

        let malformed_const_constraint_indices: Vec<_> = instance
            .ct_quad_dot_prod_funcs
            .iter()
            .enumerate()
            .filter_map(|(idx, c)| (!self.is_wellformed_const_constraint(c)).then(|| idx))
            .collect();
        if malformed_const_constraint_indices.len() > 0 {
            return Err(anyhow::anyhow!(
                "Constant constraints {:?} are not well-formed",
                malformed_const_constraint_indices
            ));
        }
        return Ok(());
    }

    pub fn is_wellformed_witness(&self, witness: &Witness<R>) -> anyhow::Result<()> {
        if witness.s.len() != self.r {
            return Err(anyhow::anyhow!(
                "Number of witness vectors does not match r"
            ));
        }
        let malformed_witness_indices: Vec<_> = witness
            .s
            .iter()
            .enumerate()
            .filter_map(|(idx, s)| (s.len() != self.n).then(|| idx))
            .collect();
        if malformed_witness_indices.len() > 0 {
            return Err(anyhow::anyhow!(
                "Witness vectors {:?} have the wrong length",
                malformed_witness_indices
            ));
        }
        if witness
            .s
            .iter()
            .map(|s_i| s_i.l2_norm_squared())
            .sum::<BigUint>()
            .to_f64()
            .unwrap()
            > self.norm_bound_squared
        {
            return Err(anyhow::anyhow!(
                "L2 norm of witness vectors exceeds norm_bound_squared"
            ));
        }
        Ok(())
    }
}

pub struct PrincipalRelation<R: PolyRing> {
    _marker: std::marker::PhantomData<R>,
}

impl<R: PolyRing> Relation for PrincipalRelation<R> {
    type Size = Size;
    type Index = Index<R>;
    type Instance = Instance<R>;
    type Witness = Witness<R>;

    fn is_well_defined(pp: &Self::Index, x: &Self::Instance, w: Option<&Self::Witness>) -> bool {
        Self::is_well_defined_err(pp, x, w).is_ok()
    }

    fn is_well_defined_err(
        pp: &Self::Index,
        x: &Self::Instance,
        w: Option<&Self::Witness>,
    ) -> anyhow::Result<()> {
        pp.is_wellformed_instance(x)?;
        match w {
            Some(w) => pp.is_wellformed_witness(w),
            None => Ok(()),
        }
    }

    fn is_satisfied(pp: &Self::Index, x: &Self::Instance, w: &Self::Witness) -> bool {
        Self::is_satisfied_err(pp, x, w).is_ok()
    }

    fn is_satisfied_err(
        pp: &Self::Index,
        x: &Self::Instance,
        w: &Self::Witness,
    ) -> anyhow::Result<()> {
        Self::is_well_defined_err(pp, x, Some(w)).and(x.is_valid_witness(w))
    }

    fn generate_satisfied_instance(
        size: &Self::Size,
    ) -> (Self::Index, Self::Instance, Self::Witness) {
        assert!(size.num_witnesses > 0, "Need at least one witness");
        assert!(size.witness_size > 0, "Need positive witness size");
        assert!(
            size.num_constraints > 0 || size.num_constant_constraints > 0,
            "Need at least one constraint"
        );

        // TODO: implement a more interesting satisfied relation
        let index = Index::<R>::new(size);
        let instance = Instance::<R> {
            quad_dot_prod_funcs: (0..index.num_constraints)
                .map(|_| {
                    QuadDotProdFunction::<R>::new_dummy(
                        size.num_witnesses,
                        size.witness_size,
                        &mut thread_rng(),
                    )
                })
                .collect(),
            ct_quad_dot_prod_funcs: (0..index.num_constant_constraints)
                .map(|_| {
                    ConstantQuadDotProdFunction::<R>::new_dummy(
                        size.num_witnesses,
                        size.witness_size,
                        &mut thread_rng(),
                    )
                })
                .collect(),
        };

        let witness = Witness::zero(size.num_witnesses, size.witness_size);
        (index, instance, witness)
    }

    fn generate_unsatisfied_instance(
        size: &Self::Size,
    ) -> (Self::Index, Self::Instance, Self::Witness) {
        assert!(size.num_witnesses > 0, "Need at least one witness");
        assert!(size.witness_size > 0, "Need positive witness size");
        assert!(
            size.num_constraints > 0 || size.num_constant_constraints > 0,
            "Need at least one constraint"
        );

        // TODO: implement a more interesting satisfied relation
        let index = Index::<R>::new(size);
        let mut instance = Instance::<R> {
            quad_dot_prod_funcs: (0..index.num_constraints)
                .map(|_| {
                    QuadDotProdFunction::<R>::new_dummy(
                        size.num_witnesses,
                        size.witness_size,
                        &mut thread_rng(),
                    )
                })
                .collect(),
            ct_quad_dot_prod_funcs: (0..index.num_constant_constraints)
                .map(|_| {
                    ConstantQuadDotProdFunction::<R>::new_dummy(
                        size.num_witnesses,
                        size.witness_size,
                        &mut thread_rng(),
                    )
                })
                .collect(),
        };

        // Add single unsatisfied constraint
        instance.ct_quad_dot_prod_funcs.last_mut().unwrap().b = R::BaseRing::one();

        let witness = Witness::zero(size.num_witnesses, size.witness_size);
        (index, instance, witness)
    }
}

#[cfg(test)]
mod test {
    use lattirust_arithmetic::ring::ntt::ntt_prime;
    use lattirust_arithmetic::ring::{Pow2CyclotomicPolyRingNTT, Zq1};

    use crate::principal_relation::{PrincipalRelation, Size};
    use crate::Relation;
    use crate::{test_generate_satisfied_instance, test_generate_unsatisfied_instance};

    const Q: u64 = ntt_prime::<64>(32);
    const D: usize = 64;

    type BaseRing = Zq1<Q>;
    type R = Pow2CyclotomicPolyRingNTT<BaseRing, D>;
    type RELATION = PrincipalRelation<R>;

    const TEST_SIZE: Size = Size {
        num_witnesses: 10,
        witness_size: 128,
        norm_bound_sq: (Q as f64) / 10.,
        num_constraints: 32,
        num_constant_constraints: 64,
    };

    test_generate_satisfied_instance!(RELATION, TEST_SIZE);

    test_generate_unsatisfied_instance!(RELATION, TEST_SIZE);
}
