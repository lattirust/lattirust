#![allow(non_snake_case)]

use ark_std::rand;
use ark_std::rand::Rng;
use ark_std::rand::thread_rng;
use derive_more::Display;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;

use lattirust_arithmetic::linear_algebra::inner_products::inner_products;
use lattirust_arithmetic::linear_algebra::SymmetricMatrix;
use lattirust_arithmetic::linear_algebra::Vector;
use lattirust_arithmetic::ring::PolyRing;
use lattirust_arithmetic::ring::representatives::WithSignedRepresentative;
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

#[derive(Clone, PartialEq, Debug, Display)]
#[display(
    "PrincipalRelation::Instance: \nquad_dot_prod_funcs: {:?}\nct_quad_dot_prod_funcs: {:?}",
    quad_dot_prod_funcs,
    ct_quad_dot_prod_funcs
)]
pub struct Instance<R: PolyRing> {
    pub quad_dot_prod_funcs: Vec<QuadraticConstraint<R>>,
    pub ct_quad_dot_prod_funcs: Vec<ConstantQuadraticConstraint<R>>,
}

#[derive(Clone, Debug, PartialEq, Display)]
#[display("QuadraticConstraint: A: {:?}, phi: {:?}, b: {:?}", A, phi, b)]
pub struct QuadraticConstraint<R: PolyRing> {
    // TODO: A is always symmetric, so we could at least use a symmetric matrix type. A is also very sparse in some cases.
    pub A: Option<SymmetricMatrix<R>>,
    // TODO: phi can be quite sparse
    pub phi: Vec<Vector<R>>,
    pub b: R,
    _private: (), // Forbid direct initialization, force users to use new(), which does some basis debug_asserts
}

impl<R: PolyRing> QuadraticConstraint<R> {
    pub fn new(A: SymmetricMatrix<R>, phi: Vec<Vector<R>>, b: R) -> Self {
        let (r, n) = (A.size(), phi[0].len());
        debug_assert_eq!(
            phi.len(),
            r,
            "phi should have the same length as the dimensions of A, got {} and {}",
            phi.len(),
            r
        );
        debug_assert!(
            phi.iter().all(|phi_i| phi_i.len() == n),
            "each phi_i should have the same length, got lengths {:?}",
            phi.iter().map(|v| v.len()).collect::<Vec<_>>()
        );
        Self {
            A: Some(A),
            phi,
            b,
            _private: (),
        }
    }

    pub fn eval(&self, witness: &Witness<R>) -> R {
        let inner_prods = inner_products(&witness.s); // TODO: don't recompute every time

        let mut res = -self.b;
        if let Some(A) = &self.A {
            let r = A.size();
            for i in 0..r {
                for j in 0..r {
                    res += A[(i, j)] * inner_prods[(i, j)];
                }
            }
        }
        for i in 0..self.phi.len() {
            res += self.phi[i].dot(&witness.s[i]);
        }

        res
    }

    pub fn new_linear(phi: Vec<Vector<R>>, b: R) -> Self {
        let n = phi[0].len();
        debug_assert!(
            phi.iter().all(|phi_i| phi_i.len() == n),
            "each phi_i should have the same length, got lengths {:?}",
            phi.iter().map(|v| v.len()).collect::<Vec<_>>()
        );
        Self {
            A: None,
            phi,
            b,
            _private: (),
        }
    }

    pub fn new_homogeneous_linear(phi: Vec<Vector<R>>) -> Self {
        Self::new_linear(phi, R::zero())
    }

    pub fn new_dummy_homogeneous<Rng: rand::Rng + ?Sized>(
        r: usize,
        n: usize,
        rng: &mut Rng,
    ) -> Self {
        Self::new(
            SymmetricMatrix::<R>::rand(r, rng),
            vec![Vector::<R>::rand(n, rng); r],
            R::zero(),
        )
    }

    pub fn new_empty(r: usize, n: usize) -> Self {
        Self {
            A: None,
            phi: vec![Vector::<R>::zeros(n); r],
            b: R::zero(),
            _private: (),
        }
    }

    fn is_valid_witness(&self, witness: &Witness<R>) -> bool {
        self.eval(witness).is_zero()
    }
}

#[derive(Clone, Debug, PartialEq, Display)]
#[display("ConstantQuadraticConstraint: A: {:?}, phi: {:?}, b: {:?}", A, phi, b)]
pub struct ConstantQuadraticConstraint<R: PolyRing> {
    pub A: Option<SymmetricMatrix<R>>,
    pub phi: Vec<Vector<R>>,
    pub b: R::BaseRing,
    _private: (), // Forbid direct initialization, force users to use new(), which does some basis debug_asserts
}

impl<R: PolyRing> ConstantQuadraticConstraint<R> {
    pub fn new(A: SymmetricMatrix<R>, phi: Vec<Vector<R>>, b: R::BaseRing) -> Self {
        let (r, n) = (A.size(), phi[0].len());

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

    pub fn eval(&self, witness: &Witness<R>) -> R::BaseRing {
        let inner_prods = inner_products(&witness.s); // TODO: don't recompute every time

        let mut res = -self.b;
        if let Some(A) = &self.A {
            let r = A.size();
            for i in 0..r {
                for j in 0..r {
                    res += (A[(i, j)] * inner_prods[(i, j)]).coefficients()[0];
                }
            }
        }
        for i in 0..self.phi.len() {
            res += (self.phi[i].dot(&witness.s[i])).coefficients()[0];
        }

        res
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

    pub fn is_valid_witness(&self, witness: &Witness<R>) -> bool {
        self.eval(witness).is_zero()
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Witness<R: PolyRing> {
    pub s: Vec<Vector<R>>,
    _private: (), // Forbid direct initialization, force users to use new(), which does some basis debug_asserts
}

impl<R: PolyRing> Witness<R>
where
    R::BaseRing: WithSignedRepresentative,
{
    pub fn rand<Rng: rand::Rng + ?Sized>(
        num_witnesses: usize,
        witness_len: usize,
        norm_bound_sq: f64,
        rng: &mut Rng,
    ) -> Self {
        // Restrict the L∞ norm
        // TODO: sample uniformly from the ball of radius norm_bound instead?
        let linf_bound = (norm_bound_sq.sqrt()
            / f64::sqrt((num_witnesses * witness_len * R::dimension()) as f64))
        .floor() as i128;
        let mut s = vec![];
        for _ in 0..num_witnesses {
            s.push(Vector::<R>::from(
                (0..witness_len)
                    .map(|_| {
                        R::try_from_coefficients(
                            (0..R::dimension())
                                .map(|_| {
                                    let rand_i128 = rng.sample(rand::distributions::Uniform::<i128>::new(
                                        -linf_bound,
                                        linf_bound,
                                    ));
                                    <R::BaseRing as WithSignedRepresentative>::SignedRepresentative::try_from(
                                        rand_i128
                                    ).unwrap()
                                    .into()
                                })
                                .collect::<Vec<R::BaseRing>>()
                                .as_slice(),
                        )
                        .unwrap()
                    })
                    .collect::<Vec<R>>(),
            ));
        }
        Self::new(s)
    }
}

impl<R: PolyRing> Witness<R> {
    pub fn new(s: Vec<Vector<R>>) -> Self {
        debug_assert!(s.len() > 0, "Witness must have at least one vector");
        debug_assert!(
            s.iter().all(|v| v.len() == s[0].len()),
            "All witness vectors must have the same length, got {:?}",
            s.iter().map(|v| v.len()).collect::<Vec<_>>()
        );
        Self { s, _private: () }
    }

    pub fn zero(rank: usize, multiplicity: usize) -> Self {
        Self::new(vec![
            Vector::<R>::from_fn(multiplicity, |_, _| R::zero());
            rank
        ])
    }
}

impl<R: PolyRing> Instance<R> {
    pub fn num_constraints(&self) -> usize {
        self.quad_dot_prod_funcs.len()
    }

    pub fn num_const_constraints(&self) -> usize {
        self.ct_quad_dot_prod_funcs.len()
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Size {
    pub num_witnesses: usize,
    pub witness_len: usize,
    pub norm_bound_sq: f64,
    pub num_constraints: usize,
    pub num_constant_constraints: usize,
}

impl<R: PolyRing> Index<R> {
    pub fn new(size: &Size) -> Self {
        Self {
            r: size.num_witnesses,
            n: size.witness_len,
            norm_bound_squared: size.norm_bound_sq,
            num_constraints: size.num_constraints,
            num_constant_constraints: size.num_constant_constraints,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn is_wellformed_constraint(&self, c: &QuadraticConstraint<R>) -> bool {
        c.phi.len() == self.r
            && c.phi.iter().all(|phi_i| phi_i.len() == self.n)
            && match c.A {
                Some(ref A) => A.size() == self.r,
                None => true,
            }
    }

    pub fn is_wellformed_const_constraint(&self, c: &ConstantQuadraticConstraint<R>) -> bool {
        c.phi.len() == self.r
            && c.phi.iter().all(|phi_i| phi_i.len() == self.n)
            && match c.A {
                Some(ref A) => A.size() == self.r,
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
            .filter_map(|(idx, c)| (!self.is_wellformed_constraint(c)).then_some(idx))
            .collect();
        if !malformed_constraint_indices.is_empty() {
            return Err(anyhow::anyhow!(
                "Constraints {:?} (/{}) are not well-formed",
                malformed_constraint_indices,
                self.num_constraints
            ));
        }

        let malformed_const_constraint_indices: Vec<_> = instance
            .ct_quad_dot_prod_funcs
            .iter()
            .enumerate()
            .filter_map(|(idx, c)| (!self.is_wellformed_const_constraint(c)).then_some(idx))
            .collect();
        if !malformed_const_constraint_indices.is_empty() {
            return Err(anyhow::anyhow!(
                "Constant constraints {:?} are not well-formed",
                malformed_const_constraint_indices
            ));
        }
        Ok(())
    }

    pub fn is_wellformed_witness(&self, witness: &Witness<R>) -> anyhow::Result<()> {
        if witness.s.len() != self.r {
            return Err(anyhow::anyhow!(
                "Number of witness vectors {} does not match r={}",
                self.r,
                witness.s.len(),
            ));
        }
        let malformed_witness_indices: Vec<_> = witness
            .s
            .iter()
            .enumerate()
            .filter_map(|(idx, s)| (s.len() != self.n).then_some(idx))
            .collect();
        if !malformed_witness_indices.is_empty() {
            return Err(anyhow::anyhow!(
                "Witness vectors {:?} have the wrong lengths ({:?} instead of {})",
                malformed_witness_indices,
                witness.s.iter().map(|s| s.len()).collect::<Vec<_>>(),
                self.n
            ));
        }

        let l2_norm_sq = witness
            .s
            .iter()
            .map(|s_i| s_i.l2_norm_squared())
            .sum::<BigUint>()
            .to_f64()
            .unwrap();
        if l2_norm_sq > self.norm_bound_squared {
            return Err(anyhow::anyhow!(
                "Squared L2 norm of witness vectors was {l2_norm_sq}, which is larger than norm_bound_squared = {}",
                self.norm_bound_squared
            ));
        }
        Ok(())
    }
}

pub struct PrincipalRelation<R: PolyRing> {
    _marker: std::marker::PhantomData<R>,
}

impl<R: PolyRing> Relation for PrincipalRelation<R>
where
    R::BaseRing: WithSignedRepresentative,
{
    type Size = Size;
    type Index = Index<R>;
    type Instance = Instance<R>;
    type Witness = Witness<R>;

    fn is_well_defined(i: &Self::Index, x: &Self::Instance, w: Option<&Self::Witness>) -> bool {
        Self::is_well_defined_err(i, x, w).is_ok()
    }

    fn is_well_defined_err(
        i: &Self::Index,
        x: &Self::Instance,
        w: Option<&Self::Witness>,
    ) -> anyhow::Result<()> {
        i.is_wellformed_instance(x)?;
        match w {
            Some(w) => i.is_wellformed_witness(w),
            None => Ok(()),
        }
    }

    fn is_satisfied(i: &Self::Index, x: &Self::Instance, w: &Self::Witness) -> bool {
        Self::is_satisfied_err(i, x, w).is_ok()
    }

    fn is_satisfied_err(
        i: &Self::Index,
        x: &Self::Instance,
        w: &Self::Witness,
    ) -> anyhow::Result<()> {
        Self::is_well_defined_err(i, x, Some(w))?;

        let unsat_indices = x
            .quad_dot_prod_funcs
            .iter()
            .enumerate()
            .filter_map(|(idx, c)| (!c.is_valid_witness(w)).then_some(idx))
            .collect::<Vec<_>>();
        let ct_unsat_indices = x
            .ct_quad_dot_prod_funcs
            .iter()
            .enumerate()
            .filter_map(|(idx, c)| (!c.is_valid_witness(w)).then_some(idx))
            .collect::<Vec<_>>();

        if !unsat_indices.is_empty() && !ct_unsat_indices.is_empty() {
            return Err(anyhow::anyhow!(
                "Standard constraints {:?} (/{}) and constant-coeff constraints {:?} (/{}) are not satisfied",
                unsat_indices,
                x.num_constraints(),
                ct_unsat_indices,
                x.num_const_constraints()
            ));
        } else if !unsat_indices.is_empty() {
            return Err(anyhow::anyhow!(
                "Standard constraints {:?} (/{}) are not satisfied",
                unsat_indices,
                x.num_constraints()
            ));
        } else if !ct_unsat_indices.is_empty() {
            return Err(anyhow::anyhow!(
                "Constant-coeff constraints {:?} (/{}) are not satisfied",
                ct_unsat_indices,
                x.num_const_constraints()
            ));
        }
        Ok(())
    }

    fn generate_satisfied_instance(
        size: &Self::Size,
    ) -> (Self::Index, Self::Instance, Self::Witness) {
        assert!(size.num_witnesses > 0, "Need at least one witness");
        assert!(size.witness_len > 0, "Need positive witness size");
        assert!(
            size.num_constraints > 0 || size.num_constant_constraints > 0,
            "Need at least one constraint"
        );

        let rng = &mut thread_rng();
        let index = Index::<R>::new(size);
        let witness = Witness::rand(
            size.num_witnesses,
            size.witness_len,
            size.norm_bound_sq,
            rng,
        );
        let instance = Instance::<R> {
            quad_dot_prod_funcs: (0..index.num_constraints)
                .map(|_| {
                    let A = SymmetricMatrix::<R>::rand(index.r, rng);
                    let phi = (0..size.num_witnesses)
                        .map(|_| Vector::<R>::rand(size.witness_len, rng))
                        .collect();
                    let mut constraint = QuadraticConstraint::new(A, phi, R::zero());
                    let b = constraint.eval(&witness);
                    constraint.b = b;
                    constraint
                })
                .collect(),
            ct_quad_dot_prod_funcs: (0..index.num_constant_constraints)
                .map(|_| {
                    let A = SymmetricMatrix::<R>::rand(index.r, rng);
                    let phi = (0..size.num_witnesses)
                        .map(|_| Vector::<R>::rand(size.witness_len, rng))
                        .collect();
                    let mut constraint =
                        ConstantQuadraticConstraint::new(A, phi, R::BaseRing::zero());
                    let b = constraint.eval(&witness);
                    constraint.b = b;
                    constraint
                })
                .collect(),
        };

        (index, instance, witness)
    }

    fn generate_unsatisfied_instance(
        size: &Self::Size,
    ) -> (Self::Index, Self::Instance, Self::Witness) {
        assert!(size.num_witnesses > 0, "Need at least one witness");
        assert!(size.witness_len > 0, "Need positive witness size");
        assert!(
            size.num_constraints > 0 || size.num_constant_constraints > 0,
            "Need at least one constraint"
        );

        let sat = Self::generate_satisfied_instance(size);

        let index = sat.0;
        let mut instance = sat.1;
        let witness = sat.2;

        let rng = &mut thread_rng();

        instance.quad_dot_prod_funcs[rng.gen_range(0..index.num_constraints)].b += R::one();
        instance.ct_quad_dot_prod_funcs[rng.gen_range(0..index.num_constant_constraints)].b +=
            R::BaseRing::one();

        (index, instance, witness)
    }
}

#[cfg(test)]
mod test {
    use lattirust_arithmetic::ring::{Pow2CyclotomicPolyRingNTT, Zq1};
    use lattirust_arithmetic::ring::ntt::ntt_prime;

    use crate::{test_generate_satisfied_instance, test_generate_unsatisfied_instance};
    use crate::principal_relation::{PrincipalRelation, Size};
    use crate::Relation;

    const Q: u64 = ntt_prime::<64>(32);
    const D: usize = 64;

    type BaseRing = Zq1<Q>;
    type R = Pow2CyclotomicPolyRingNTT<BaseRing, D>;
    type RELATION = PrincipalRelation<R>;

    const TEST_SIZE: Size = Size {
        num_witnesses: 10,
        witness_len: 128,
        norm_bound_sq: (Q as f64) / 100.,
        num_constraints: 16,
        num_constant_constraints: 8,
    };

    test_generate_satisfied_instance!(RELATION, TEST_SIZE);

    test_generate_unsatisfied_instance!(RELATION, TEST_SIZE);
}
