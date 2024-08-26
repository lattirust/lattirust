use ark_std::rand::rngs::OsRng;

use lattirust_arithmetic::linear_algebra::{SparseMatrix, Vector};
use lattirust_arithmetic::ring::Ring;

use crate::Relation;

pub struct R1CS<R: Ring> {
    _marker: std::marker::PhantomData<R>,
}

pub struct Index<R: Ring> {
    a: SparseMatrix<R>,
    b: SparseMatrix<R>,
    c: SparseMatrix<R>,
}

pub struct Instance<R: Ring>(Vec<R>);

pub struct Witness<R: Ring>(Vec<R>);

pub struct Size {
    num_constraints: usize,
    num_public_variables: usize,
    num_private_variables: usize,
}

impl<R: Ring> Relation for R1CS<R> {
    type Size = Size;
    type Index = Index<R>;
    type Instance = Instance<R>;
    type Witness = Witness<R>;

    fn is_well_defined(i: &Self::Index, x: &Self::Instance, w: Option<&Self::Witness>) -> bool {
        let matrices_same_dim = i.a.nrows() == i.b.nrows()
            && i.b.nrows() == i.c.nrows()
            && i.a.ncols() == i.b.ncols()
            && i.b.ncols() == i.c.ncols();
        match w {
            Some(w) => matrices_same_dim && i.a.ncols() == x.0.len() + w.0.len(),
            None => matrices_same_dim && i.a.ncols() <= x.0.len(),
        }
    }

    fn is_satisfied(i: &Self::Index, x: &Self::Instance, w: &Self::Witness) -> bool {
        let z = Vector::<R>::from_vec(
            x.0.clone()
                .into_iter()
                .chain(w.0.clone().into_iter())
                .collect::<Vec<R>>(),
        );
        Self::is_well_defined(i, x, Some(w)) && {
            let a_z = &i.a * &z;
            let b_z = &i.b * &z;
            let c_z = &(i.c) * &z;
            a_z.component_mul(&b_z) == c_z
        }
    }

    fn generate_satisfied_instance(
        size: &Self::Size,
    ) -> (Self::Index, Self::Instance, Self::Witness) {
        assert!(size.num_private_variables > 0 || size.num_public_variables > 0);
        
        let mut csprng: OsRng = OsRng;
        let num_variables = size.num_public_variables + size.num_private_variables;

        let mut z: Vec<R> = (0..num_variables).map(|_| R::rand(&mut csprng)).collect();
        z[0] = R::one(); // set the constant term to 1

        let mut a_triplets = Vec::with_capacity(size.num_constraints);
        let mut b_triplets = Vec::with_capacity(size.num_constraints);
        let mut c_triplets = Vec::with_capacity(size.num_constraints);

        for i in 0..size.num_constraints {
            let a_idx = i % num_variables;
            let b_idx = (i + 1) % num_variables;
            let c_idx = (i + 2) % num_variables;

            a_triplets.push((i, a_idx, R::one()));
            b_triplets.push((i, b_idx, R::one()));

            let ab_val = z[a_idx] * z[b_idx];
            match z[c_idx].inverse() {
                Some(c_val_inv) => c_triplets.push((i, c_idx, ab_val * c_val_inv)),
                None => c_triplets.push((i, 0, ab_val)),
            };
        }

        let index = Index {
            a: SparseMatrix::try_from_triplets(size.num_constraints, num_variables, a_triplets)
                .unwrap(),
            b: SparseMatrix::try_from_triplets(size.num_constraints, num_variables, b_triplets)
                .unwrap(),
            c: SparseMatrix::try_from_triplets(size.num_constraints, num_variables, c_triplets)
                .unwrap(),
        };
        let instance = Instance(z[..size.num_public_variables].to_vec());
        let witness = Witness(z[size.num_public_variables..].to_vec());

        debug_assert!(Self::is_well_defined(&index, &instance, Some(&witness)));
        debug_assert!(Self::is_satisfied(&index, &instance, &witness));
        (index, instance, witness)
    }

    fn generate_unsatisfied_instance(
        size: &Self::Size,
    ) -> (Self::Index, Self::Instance, Self::Witness) {
        assert!(size.num_private_variables > 0 || size.num_public_variables > 0);
        
        let mut csprng: OsRng = OsRng;
        let num_variables = size.num_public_variables + size.num_private_variables;

        let mut z: Vec<R> = (0..num_variables).map(|_| R::rand(&mut csprng)).collect();
        z[0] = R::one(); // set the constant term to 1

        let mut a_triplets = Vec::with_capacity(size.num_constraints);
        let mut b_triplets = Vec::with_capacity(size.num_constraints);
        let mut c_triplets = Vec::with_capacity(size.num_constraints);

        for i in 0..size.num_constraints - 1 {
            let a_idx = i % num_variables;
            let b_idx = (i + 1) % num_variables;
            let c_idx = (i + 2) % num_variables;

            a_triplets.push((i, a_idx, R::one()));
            b_triplets.push((i, b_idx, R::one()));

            let ab_val = z[a_idx] * z[b_idx];
            match z[c_idx].inverse() {
                Some(c_val_inv) => c_triplets.push((i, c_idx, ab_val * c_val_inv)),
                None => c_triplets.push((i, 0, ab_val)),
            };
        }

        // Insert single unsatisfiable constraint; 0x * 1x == 1
        let unsatisfied_idx = if size.num_private_variables > 0 {
            size.num_public_variables
        } else {
            0
        };
        a_triplets.push((size.num_constraints - 1, unsatisfied_idx, R::zero()));
        b_triplets.push((size.num_constraints - 1, unsatisfied_idx, R::one()));
        c_triplets.push((size.num_constraints - 1, 0, R::one()));

        let index = Index {
            a: SparseMatrix::try_from_triplets(size.num_constraints, num_variables, a_triplets)
                .unwrap(),
            b: SparseMatrix::try_from_triplets(size.num_constraints, num_variables, b_triplets)
                .unwrap(),
            c: SparseMatrix::try_from_triplets(size.num_constraints, num_variables, c_triplets)
                .unwrap(),
        };
        let instance = Instance(z[..size.num_public_variables].to_vec());
        let witness = Witness(z[size.num_public_variables..].to_vec());

        debug_assert!(Self::is_well_defined(&index, &instance, Some(&witness)));
        debug_assert!(!Self::is_satisfied(&index, &instance, &witness));
        (index, instance, witness)
    }
}

#[cfg(test)]
mod test {
    use lattirust_arithmetic::ring::Zq;

    use crate::{test_generate_satisfied_instance, test_generate_unsatisfied_instance};
    use crate::Relation;

    use super::*;

    const Q: u64 = 4294967297; // fifth fermat number, non-prime
    type R = Zq<Q>;
    type RELATION = R1CS<R>;

    const TEST_SIZE: Size = Size {
        num_constraints: 512,
        num_public_variables: 64,
        num_private_variables: 128,
    };

    test_generate_satisfied_instance!(RELATION, TEST_SIZE);

    test_generate_unsatisfied_instance!(RELATION, TEST_SIZE);
}
