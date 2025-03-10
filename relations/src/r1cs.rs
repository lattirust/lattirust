use std::ops::AddAssign;

use ark_std::rand::rngs::OsRng;
use num_traits::Zero;

use lattirust_arithmetic::linear_algebra::{Scalar, SparseMatrix, Vector};
use lattirust_arithmetic::ring::Ring;

use crate::Relation;

pub struct R1CS<R: Ring> {
    _marker: std::marker::PhantomData<R>,
}

pub struct Index<R: Ring> {
    pub a: SparseMatrix<R>,
    pub b: SparseMatrix<R>,
    pub c: SparseMatrix<R>,
}

pub struct Instance<R: Ring>(pub Vec<R>);

pub struct Witness<R: Ring>(pub Vec<R>);

pub struct Size {
    pub num_constraints: usize,
    pub num_instance_variables: usize,
    pub num_witness_variables: usize,
}

impl<R: Ring> Relation for R1CS<R> {
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
        let matrices_same_dim = i.a.nrows() == i.b.nrows()
            && i.b.nrows() == i.c.nrows()
            && i.a.ncols() == i.b.ncols()
            && i.b.ncols() == i.c.ncols();
        if !matrices_same_dim {
            anyhow::bail!("Matrices A, B, and C must have the same dimensions");
        }
        match w {
            Some(w) => {
                if i.a.ncols() != x.0.len() + w.0.len() {
                    anyhow::bail!("The number of columns in A must be equal to the sum of the number of instance and witness variables");
                } else {
                    Ok(())
                }
            }
            None => {
                if i.a.ncols() != x.0.len() {
                    anyhow::bail!("The number of columns in A must be equal to the number of instance variables");
                } else {
                    Ok(())
                }
            }
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

        let z = Vector::<R>::from_vec(
            x.0.clone()
                .into_iter()
                .chain(w.0.clone())
                .collect::<Vec<R>>(),
        );

        let a_z = &i.a * &z;
        let b_z = &i.b * &z;
        let c_z = &i.c * &z;

        if a_z.component_mul(&b_z) != c_z {
            anyhow::bail!("The R1CS constraint is not satisfied");
        } else {
            Ok(())
        }
    }

    fn generate_satisfied_instance(
        size: &Self::Size,
    ) -> (Self::Index, Self::Instance, Self::Witness) {
        assert!(size.num_witness_variables > 0 || size.num_instance_variables > 0);

        let mut csprng: OsRng = OsRng;
        let num_variables = size.num_instance_variables + size.num_witness_variables;

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
        let instance = Instance(z[..size.num_instance_variables].to_vec());
        let witness = Witness(z[size.num_instance_variables..].to_vec());

        debug_assert!(Self::is_well_defined(&index, &instance, Some(&witness)));
        debug_assert!(Self::is_satisfied(&index, &instance, &witness));
        (index, instance, witness)
    }

    fn generate_unsatisfied_instance(
        size: &Self::Size,
    ) -> (Self::Index, Self::Instance, Self::Witness) {
        assert!(size.num_witness_variables > 0 || size.num_instance_variables > 0);

        let mut csprng: OsRng = OsRng;
        let num_variables = size.num_instance_variables + size.num_witness_variables;

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

        // Insert single unsatisfiable constraint; 0 * 1 == 1
        a_triplets.push((size.num_constraints - 1, 0, R::zero()));
        b_triplets.push((size.num_constraints - 1, 0, R::one()));
        c_triplets.push((size.num_constraints - 1, 0, R::one()));

        let index = Index {
            a: SparseMatrix::try_from_triplets(size.num_constraints, num_variables, a_triplets)
                .unwrap(),
            b: SparseMatrix::try_from_triplets(size.num_constraints, num_variables, b_triplets)
                .unwrap(),
            c: SparseMatrix::try_from_triplets(size.num_constraints, num_variables, c_triplets)
                .unwrap(),
        };
        let instance = Instance(z[..size.num_instance_variables].to_vec());
        let witness = Witness(z[size.num_instance_variables..].to_vec());

        debug_assert!(Self::is_well_defined(&index, &instance, Some(&witness)));
        debug_assert!(!Self::is_satisfied(&index, &instance, &witness));
        (index, instance, witness)
    }
}

pub fn sparse_matrix_from_ark_matrix<R: Scalar + Copy + Zero + AddAssign>(
    matrix: ark_relations::r1cs::Matrix<R>,
    nrows: usize,
    ncols: usize,
) -> SparseMatrix<R> {
    assert_eq!(nrows, matrix.len());
    let triplets = matrix
        .iter()
        .enumerate()
        .flat_map(|(row_index, row)| {
            row.iter()
                .map(move |(elem, col_index)| (row_index, *col_index, *elem))
        })
        .collect();
    SparseMatrix::<R>::try_from_triplets(nrows, ncols, triplets).unwrap()
}

#[cfg(test)]
mod test {
    use lattirust_arithmetic::ring::Zq1;

    use crate::test_generate_unsatisfied_instance;
    use crate::{test_generate_satisfied_instance, Relation};

    use super::*;

    const Q: u64 = 294967297; // fifth fermat number, non-prime
    type R = Zq1<Q>;
    type RELATION = R1CS<R>;

    const TEST_SIZE: Size = Size {
        num_constraints: 1,
        num_instance_variables: 1,
        num_witness_variables: 1,
    };

    test_generate_satisfied_instance!(RELATION, TEST_SIZE);

    test_generate_unsatisfied_instance!(RELATION, TEST_SIZE);
}
