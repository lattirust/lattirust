#![allow(non_snake_case)]

use ark_linear_sumcheck::ml_sumcheck::MLSumcheck;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::iterable::Iterable;
use nimue::IOPatternError;
use num_traits::One;

use crate::lattice_arithmetic::balanced_decomposition::decompose_balanced_vec;
use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::ntt::NTT;
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::traits::FromRandomBytes;
use crate::latticefold::util::*;
use crate::nimue::lattice_arthur::LatticeArthur;
use crate::relations::{cm_ccs, eval_ccs};
use crate::sumcheck::boolean_hypercube::BooleanHypercube;

// R_cmccs^B -> R_evalccs^B
pub fn linearize(arthur: &mut LatticeArthur<R>, crs: &cm_ccs::Crs<R>, x: &cm_ccs::Instance<R>, w: &cm_ccs::Witness<R>) -> Result<(), IOPatternError>
{
    // Params: sampling set C = Z_q
    // Input:
    // - cm in R_q^k
    // - x_css in R_q^(l_in)
    // - f in R_q^m
    // - w_ccs: R_q^(n_c - l_in - 1)
    // Output:
    // - r_o in R_q^log(m)
    // - v in R_q
    // - cm
    // - u in R_q^t
    // - w_o = (f, w_css)

    // Get random beta in C^log(m), C = Z_q
    let log_m = crs.crs_cm.m.next_power_of_two().ilog2() as usize;
    let beta = arthur.squeeze_vector::<F, F>(log_m)?;

    // Construct sumcheck polynomial
    let gs = linearization_sumcheck_poly(crs, &x, &w, &beta);

    // Run sumcheck
    let r_o = gs.iter().map(|g| {
        // let mut fs_rng = Blake2s512Rng::setup();
        let mut fs_rng = arthur;
        let (proof, prover_state) = MLSumcheck::prove_as_subprotocol(&mut fs_rng, g)?;
        // TODO: feed prover messages to transcript?
        prover_state.randomness
    }).collect::<Vec<_>>();
    debug_assert!(r_o.len() == R::dimension());
    debug_assert!(r_o.iter().all(|r_i| r_i.len() == log_m));

    // Send v = mle[f^](r_o)
    let mle_f = mle_vec(&w.f);
    let v = mle_f.iter().zip(r_o).map(|(f_i, r_o_i)| f_i.evaluate(r_o_i.as_slice()).unwrap()).collect();
    arthur.absorb_vec(&v)?;

    // Send u_i = sum_{b in {0,1}^log(n_c)} mle[M_i](r_o, b) * mle[z_ccs](b)
    let t = crs.crs_ccs.n_matrices;
    let log_nc = crs.crs_ccs.n_cols.next_power_of_two().ilog2() as u8;
    let mle_z_css = mle_vec(&w.w_ccs); // TODO: or z_ccs?
    let u = (0..t).map(|i| {
        let mle_Mi = mle_mat(&crs.crs_ccs.matrices[i], false).iter().zip(&r_o).map(
            |(mle_Mi_j, r_o_j)| mle_Mi_j.fix_variables(r_o_j)
        );

        mle_Mi.iter().zip(mle_z_css).map(
            |(mle_Mi_j, mle_z_css_j)|
                BooleanHypercube::new(log_nc).iter().map(
                    |b| mle_Mi_j.evaluate(&b).unwrap() * mle_z_css_j.evaluate(&b).unwrap()
                ).sum()
        ).collect::<Vec<_>>()
    }
    ).collect::<Vec<_>>();

    arthur.absorb_vectors(&u)?;

    Ok(())
}

pub fn split<R: PolyRing>(b: usize, k: usize, f: &Vector<R>) -> Vec<Vector<R>> {
    // TODO: assert b^k = B
    // TODO: assert B < q/2
    // TODO: ||f||_inf < B
    decompose_balanced_vec(&f, b as u128, Some(k))
}

pub fn Lw<R: PolyRing>(f: &Vector<R>) -> Vector<R> {
    // Gf
    todo!()
}

// R_ccshom^B -> (R_ccshom^b)^k
pub fn decompose<R: PolyRing>(arthur: &mut LatticeArthur<R>, crs: &eval_ccs::Crs<crate::latticefold::util::R>, r: Vector<R>, v: R, y: R, u: R, x_w: &Vector<R>, f: &Vector<R>, w: &Vector<R>) {
    // x = (r, v, y, u, x_w)
    // w = (f, w)
    // out: x_i = (r, v_i, y_i, u_i, x_{w,i})
    // out: w_i = (f_i, w_i)

    // (f_0, ..., f_{k-1}) = split_{b,k}(f)
    let b = 0; //TODO
    let k = 0; //TODO
    let f_decomp = split(b, k, &f);

    // Send y_i = L(f_i)
    let y = f_decomp.iter().map(|f_i| &crs.crs_cm_ccs.crs_cm.ck * f_i).collect::<Vec<Vector<R>>>();
    arthur.absorb_vectors(&y)?;

    // Send v_i = mle[f_i](r)
    let v = f_decomp.iter().map(|f_i| {
        let mle_f_i = mle_vec(f_i);
        mle_f_i.iter().map(
            // TODO: split r by slice
            |f_i_j| f_i_j.evaluate(r.as_slice()).unwrap()
        ).collect()
    }).collect::<Vec<Vec<Vec<F>>>>();

    let M = Matrix::<R>::zeros();
    let r_tensor = tensor(r);

    let Lw_f = f_decomp.iter().map(|f_i|
                                       &crs.crs_cm_ccs.gadget_matrix * f_i // Lw(f_i)
    ).collect::<Vec<_>>();
    let u = Lw_f.iter().map(|Lw_f_i|
        (&M * &r_tensor).dot(Lw_f_i)
    ).collect::<Vec<_>>();

    let n_in = 0; // TODO
    let x_w = Lw_f.iter().map(|Lw_f_i|
        Vector::<R>::from_row_slice(&Lw_f_i.as_slice().to_vec()[0..n_in])
    ).collect::<Vec<_>>();

    // arthur.absorb(y, v, u, x_w)
}

fn sumcheck_poly_fold<R: PolyRing>(alpha: &Vec<R::BaseRing>, zeta: &Vec<R::BaseRing>, mu: &Vec<R::BaseRing>, beta: &Vec<R::BaseRing>, u: &Vector<R>, v: &Vector<R>) -> Vec<DenseMultilinearExtension<F>> {
    // TODO
    // Compute claimed sum
    let claimed_sum: Vec<R> = alpha.iter().zip(v).map(|(x, y)| x * y) +
        zeta.iter().zip(u).map(|(x, y)| x * y);

    // g_1,i[j](X_1, ..., X_logm) =eq(\vec{r_i,j}, \vec{X}) * mle[f_i](\vec{X})
    let g_1: Vec<F> = todo!();
}



// (R_evalccs^b)^2k -> R_evalccs^B
pub fn fold<R: PolyRing>(crs: &eval_ccs::Crs<R>, arthur: &mut LatticeArthur<R>)
    // where R::BaseRing: FromRandomBytes<C_small<R::BaseRing>>
{
    let n_in = 0; // TODO

    // Get verifier challenges
    let k = 0; // TODO
    let alpha = arthur.squeeze_vec::<R::BaseRing, R::BaseRing>(2 * k)?;
    let zeta = arthur.squeeze_vec::<R::BaseRing, R::BaseRing>(2 * k)?;
    let mu = arthur.squeeze_vec::<R::BaseRing, R::BaseRing>(2 * k - 1)?;
    let beta = arthur.squeeze_vec::<R::BaseRing, R::BaseRing>(crs.log_m())?;

    // Build sumcheck polynomial and claim
    let g = todo!();

    // Run sumcheck

    // Send v = mle[f^](r_o)
    // Send eta_i

    // Get verifier challenges
    // let mut rho = arthur.squeeze_vec::<C_small, C_small>(2 * crs.k - 1)?;
    // rho.splice(0..0, R::BaseRing::one());

    // Linearly combine
    // let f_o = f_decomp.iter().zip(rho.iter()).map(|(f_i, rho_i)| f_i * rho_i).fold(Vector::<R>::zeros(f_decomp[0].len()), |acc, f_i| acc + f_i);
    // let w_o = Lw(f_o).as_slice().to_vec()[n_in + 1..n_in + 1]; // TODO: or just Lw applied to f_o, since we're not actually committing to x
}