#![allow(non_snake_case)]

use lattirust_arithmetic::decomposition::DecompositionFriendlySignedRepresentative;
use lattirust_arithmetic::linear_algebra::{Matrix, SymmetricMatrix, Vector};
use lattirust_arithmetic::ring::representatives::WithSignedRepresentative;
use lattirust_arithmetic::ring::PolyRing;
use relations::principal_relation::{Index, Instance, QuadraticConstraint, Size};

use crate::common_reference_string::{CommonReferenceString, FoldedSize};
use crate::util::{flatten_symmetric_matrix, mul_basescalar_vector};

/// A view of the transcript of one execution of the core Labrador protocol
pub struct TranscriptView<R: PolyRing> {
    pub(crate) u_1: Vector<R>,
    pub(crate) b__: Vec<R>,
    pub(crate) alpha: Vector<R>,
    pub(crate) beta: Vector<R>,
    pub(crate) u_2: Vector<R>,
    pub(crate) c: Vec<R>,
    pub(crate) phi: Vec<Vector<R>>,
    pub(crate) a__: Vec<SymmetricMatrix<R>>,
}

pub struct Layouter<R: PolyRing> {
    pub folded_size: FoldedSize,
    pub vec: Vec<R>,
}

impl<R: PolyRing> Layouter<R> {
    #[inline(always)]
    fn z1_offset(folded_size: FoldedSize) -> usize {
        folded_size.nu * folded_size.size.witness_len
    }

    #[inline(always)]
    fn zi_len(folded_size: FoldedSize) -> usize {
        folded_size.size_z / 2
    }

    #[inline(always)]
    fn t_offset(folded_size: FoldedSize) -> usize {
        // t starts at 2*nu * n_next
        2 * folded_size.nu * folded_size.size.witness_len
    }

    #[inline(always)]
    fn g_offset(folded_size: FoldedSize) -> usize {
        Self::t_offset(folded_size) + folded_size.size_t
    }

    #[inline(always)]
    fn h_offset(folded_size: FoldedSize) -> usize {
        Self::g_offset(folded_size) + folded_size.size_g
    }

    #[inline(always)]
    pub fn new(folded_size: FoldedSize) -> Self {
        Self {
            folded_size,
            vec: vec![
                R::zero();
                (2 * folded_size.nu + folded_size.mu) * folded_size.size.witness_len
            ],
        }
    }

    #[inline(always)]
    pub fn set_z0(&mut self, v: &[R]) {
        self.vec[0..Self::zi_len(self.folded_size)].copy_from_slice(&v);
    }

    #[inline(always)]
    pub fn set_z1(&mut self, v: &[R]) {
        self.vec[Self::z1_offset(self.folded_size)
            ..Self::z1_offset(self.folded_size) + Self::zi_len(self.folded_size)]
            .copy_from_slice(&v);
    }

    #[inline(always)]
    pub fn set_t(&mut self, v: &[R]) {
        self.vec[Self::t_offset(self.folded_size)
            ..Self::t_offset(self.folded_size) + self.folded_size.size_t]
            .copy_from_slice(&v);
    }

    #[inline(always)]
    pub fn set_g(&mut self, v: &[R]) {
        self.vec[Self::g_offset(self.folded_size)
            ..Self::g_offset(self.folded_size) + self.folded_size.size_g]
            .copy_from_slice(&v);
    }

    pub fn set_h(&mut self, v: &[R]) {
        self.vec[Self::h_offset(self.folded_size)
            ..Self::h_offset(self.folded_size) + self.folded_size.size_h]
            .copy_from_slice(&v);
    }

    pub fn split(&self) -> Vec<Vector<R>> {
        self.vec
            .chunks_exact(self.folded_size.size.witness_len)
            .map(|chunk| Vector::<R>::from_slice(chunk))
            .collect()
    }
}

// TODO: add tracing info with size of padding to enable efficiency fine-tuning
pub fn fold_instance<R: PolyRing>(
    crs: &CommonReferenceString<R>,
    instance: &Instance<R>,
    transcript: &TranscriptView<R>,
) -> (Index<R>, Instance<R>)
where
    <R as PolyRing>::BaseRing: WithSignedRepresentative,
    <<R as PolyRing>::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
{
    // Generate instance for next iteration of the protocol
    let next_size = crs.next_size();

    let r_next = next_size.size.num_witnesses;
    let n_next = next_size.size.witness_len;
    let nu = next_size.nu;

    let mut quad_dot_prod_funcs_next =
        Vec::<QuadraticConstraint<R>>::with_capacity(next_size.size.num_constraints);

    let b_ring = R::try_from(crs.b).unwrap();

    let b1_ring = R::try_from(crs.b1).unwrap();
    let mut b1_pows = Vec::<R>::with_capacity(crs.t1);
    b1_pows.push(R::one());
    for i in 1..crs.t1 {
        b1_pows.push(b1_ring * b1_pows[i - 1]);
    }

    let b2_ring = R::try_from(crs.b2).unwrap();
    let mut b2_pows = Vec::<R>::with_capacity(crs.t2);
    b2_pows.push(R::one());
    for i in 1..crs.t2 {
        b2_pows.push(b2_ring * b2_pows[i - 1]);
    }

    // Constraints for Az = sum_{i in [r]} c_i t_i
    {
        for l in 0..crs.k {
            let mut layouter = Layouter::<R>::new(next_size);
            let A_l = crs.A.row(l).transpose();
            layouter.set_z0(A_l.as_slice()); // <A_l, z_0>
            layouter.set_z1((&A_l * b_ring).as_slice()); // <A_l * b, z_1>

            // <(t_i^(k)_l)_{k in [t1], i in [r]}, (c_i * b^k)_{k in [t1], i in [r]})> = sum_{k in [t1]} sum_{i in [r]} t_i^(k)_l * c_i * b^k
            let mut c_vec = vec![R::zero(); next_size.size_t];
            for k in 0..crs.t1 {
                for i in 0..crs.r {
                    c_vec[k * crs.r * crs.k + i * crs.k + l] = -transcript.c[i] * b1_pows[k];
                }
            }
            layouter.set_t(c_vec.as_slice());

            quad_dot_prod_funcs_next.push(QuadraticConstraint::<R>::new_homogeneous_linear(
                layouter.split(),
            ));
        }
    }

    let c_prods = SymmetricMatrix::<R>::from_fn(crs.r, |i, j| {
        if i == j {
            transcript.c[i] * transcript.c[j]
        } else {
            transcript.c[i] * transcript.c[j] + transcript.c[j] * transcript.c[i]
        }
    });
    let c_prods_vec = Vector::<R>::from(flatten_symmetric_matrix(&c_prods)); // r * (r+1) / 2
    let c_prods_b1: Vec<R> = (0..crs.t1)
        .flat_map(|k| (c_prods_vec.clone() * b1_pows[k]).as_slice().to_vec())
        .collect();
    let c_prods_b2: Vec<R> = (0..crs.t2)
        .flat_map(|k| (c_prods_vec.clone() * b2_pows[k]).as_slice().to_vec())
        .collect();

    // Constraint for <z, z> = sum_{i, j in [r]} g_ij c_i c_j
    {
        // <z, z> = <[z^(0) + b*z^(1)], [z^(0) + b*z^(1)]>
        let mut A = SymmetricMatrix::<R>::zero(r_next);
        let b_sq = R::try_from(crs.b * crs.b).expect(
            format!(
                "Ring R with modulus {} is too small to compute b^2 = {}^2",
                R::modulus(),
                crs.b
            )
            .as_str(),
        );
        for i in 0..nu {
            A[(i, i)] = -R::one();
            A[(i, i + nu)] = -b_ring; // sets A[(i+nu, i)] as well
            A[(i + nu, i + nu)] = -b_sq;
        }

        let mut layouter = Layouter::<R>::new(next_size);
        layouter.set_g(&c_prods_b2);

        quad_dot_prod_funcs_next.push(QuadraticConstraint::<R>::new(
            A,
            layouter.split(),
            R::zero(),
        ));
    }

    // Constraint for sum_{i in [r]} <phi_i, z> c_i = sum_{i, j in [r]} h_ij c_i c_j
    {
        let mut layouter = Layouter::<R>::new(next_size);

        let phi_lc_0: Vector<R> = transcript
            .phi
            .clone()
            .into_iter()
            .zip(transcript.c.iter())
            .map(|(phi_i, c_i)| (phi_i * -c_i.clone()).into())
            .sum();

        layouter.set_z0(phi_lc_0.as_slice()); // <sum_{i in [r]} phi_i * c_i, z_0>
        layouter.set_z1((phi_lc_0 * b_ring).as_slice()); // <sum_{i in [r]} phi_i * c_i * b, z_1>

        // <(c_i * c_j * b^k)_{k, i, j}, (h_ij^(k))_{k, i, j}>
        layouter.set_h(&c_prods_b1);

        quad_dot_prod_funcs_next.push(QuadraticConstraint::<R>::new_homogeneous_linear(
            layouter.split(),
        ));
    }

    // Constraint for sum_{i,j in [r]} a_ij * g_ij + sum_{i in [r]} h_ii = b
    {
        // Compute a_ij = sum_{k in [K]} alpha_k * a_ij^(k) + sum_{k in [K']} beta_k * a_ij''^(k), where a_ij''^(k) = sum_{l in [L]} psi_l^(k) * a_ij'^(l)
        let a__ = &transcript.a__;

        let mut A = Matrix::<R>::zeros(crs.r, crs.r);
        for k in 0..crs.num_constraints {
            if let Some(ref A_k) = instance.quad_dot_prod_funcs[k].A {
                let A_k_mat = Matrix::<R>::from(A_k.clone());
                A += A_k_mat * transcript.alpha[k];
            }
        }
        for k in 0..crs.num_aggregs {
            let A_k_mat = Matrix::<R>::from(a__[k].clone());
            A += &A_k_mat * transcript.beta[k];
        }

        let mut layouter = Layouter::<R>::new(next_size);
        // Set phis for a_ij * g_ij
        let mut A_vec = Vec::<R>::with_capacity(next_size.size_g);
        for k in 0..crs.t2 {
            for i in 0..crs.r {
                for j in 0..i {
                    A_vec.push((A[(i, j)] + A[(j, i)]) * b2_pows[k]);
                }
                A_vec.push(A[(i, i)] * b2_pows[k]);
            }
        }
        layouter.set_g(A_vec.as_slice());

        // Set phis for 1 * h_ii = sum_{k in [t1]} h_ii^(k) * b1^k, i.e.
        let mut indic = vec![R::zero(); next_size.size_h];
        let h_num_entries = (crs.r * (crs.r + 1)) / 2;
        for k in 0..crs.t1 {
            let mut idx = 0;
            for i in 0..crs.r {
                indic[k * h_num_entries + idx] = b1_pows[k];
                idx += crs.r - i;
            }
        }
        layouter.set_h(indic.as_slice());

        // Compute b = sum_{k in [K]} alpha_k * b^(k) + sum_{k in [K']} beta_k * b''^(k)
        let mut b = R::zero();
        for k in 0..crs.num_constraints {
            b += transcript.alpha[k] * instance.quad_dot_prod_funcs[k].b;
        }
        for k in 0..crs.num_aggregs {
            b += transcript.beta[k] * transcript.b__[k];
        }

        quad_dot_prod_funcs_next.push(QuadraticConstraint::<R>::new_linear(layouter.split(), b));
    }

    // Constraints for u_1
    {
        for l in 0..crs.k1 {
            let mut layouter = Layouter::<R>::new(next_size);
            layouter.set_t(&crs.B.row(l).transpose().as_slice());
            layouter.set_g(crs.C.row(l).transpose().as_slice());

            quad_dot_prod_funcs_next.push(QuadraticConstraint::<R>::new_linear(
                layouter.split(),
                transcript.u_1[l],
            ));
        }
    }

    // Constraints for u_2
    {
        for l in 0..crs.k2 {
            let mut layouter = Layouter::<R>::new(next_size);
            layouter.set_h(crs.D.row(l).transpose().as_slice());

            quad_dot_prod_funcs_next.push(QuadraticConstraint::<R>::new_linear(
                layouter.split(),
                transcript.u_2[l],
            ));
        }
    }

    let instance_next = Instance::<R> {
        quad_dot_prod_funcs: quad_dot_prod_funcs_next,
        ct_quad_dot_prod_funcs: vec![],
    };

    let size_next = Size {
        num_witnesses: r_next,
        witness_len: n_next,
        norm_bound_sq: CommonReferenceString::<R>::next_norm_bound_sq(
            crs.r,
            crs.n,
            crs.norm_bound_squared,
            crs.k,
            crs.b,
        ),
        num_constraints: instance_next.num_constraints(),
        num_constant_constraints: instance_next.num_const_constraints(),
    };
    let index_next = Index::<R>::new(&size_next);
    (index_next, instance_next)
}

pub fn compute_phi__<R: PolyRing>(
    crs: &CommonReferenceString<R>,
    index: &Index<R>,
    instance: &Instance<R>,
    Pi: &Vec<Matrix<R>>,
    psi: &Vec<Vector<R::BaseRing>>,
    omega: &Vec<Vector<R::BaseRing>>,
) -> Vec<Vec<Vector<R>>> {
    let mut phi__ = vec![vec![Vector::<R>::zeros(index.n); index.r]; crs.num_aggregs];
    for k in 0..crs.num_aggregs {
        for i in 0..index.r {
            // Compute vec{phi}''_i^{(k)}
            for l in 0..instance.ct_quad_dot_prod_funcs.len() {
                phi__[k][i] +=
                    mul_basescalar_vector(psi[k][l], &instance.ct_quad_dot_prod_funcs[l].phi[i]);
            }
            for j in 0..256 {
                let pi_ij = &Pi[i].row(j).transpose(); // Vector of n elements in R
                phi__[k][i] +=
                    mul_basescalar_vector(omega[k][j], &R::apply_automorphism_vec(&pi_ij));
            }
        }
    }
    phi__
}

pub fn compute_phi<R: PolyRing>(
    crs: &CommonReferenceString<R>,
    instance: &Instance<R>,
    alpha: &Vector<R>,
    beta: &Vector<R>,
    phi__: &Vec<Vec<Vector<R>>>,
) -> Vec<Vector<R>> {
    let phi = (0..crs.r)
        .into_iter()
        .map(|i| {
            let mut phi_i = Vector::<R>::zeros(crs.n);
            for k in 0..instance.quad_dot_prod_funcs.len() {
                phi_i += &instance.quad_dot_prod_funcs[k].phi[i] * alpha[k];
            }
            for k in 0..crs.num_aggregs {
                phi_i += &phi__[k][i] * beta[k];
            }
            phi_i
        })
        .collect::<Vec<_>>();
    phi
}

pub fn compute_a__<R: PolyRing>(
    _crs: &CommonReferenceString<R>,
    instance: &Instance<R>,
    psi: &Vec<Vector<R::BaseRing>>,
) -> Vec<SymmetricMatrix<R>> {
    psi.into_iter()
        .map(|psi_k| {
            psi_k
                .into_iter()
                .zip(instance.ct_quad_dot_prod_funcs.to_owned().into_iter())
                .filter_map(|(psi_k_j, constraint)| {
                    constraint.A.and_then(|A_j| Some((psi_k_j, A_j)))
                })
                .fold(None, |acc, (psi_k_j, A_j)| {
                    let prod = A_j * psi_k_j.clone();
                    match acc {
                        None => Some(prod),
                        Some(a) => Some(a + prod),
                    }
                })
                .unwrap()
        })
        .collect()
}
