mod errors;
mod multilinear_polynomial;
mod util;
mod virtual_polynomial;

//TODO add exports
pub use errors::ArithErrors;
pub use multilinear_polynomial::{
    evaluate_no_par, evaluate_opt, fix_last_variables, fix_last_variables_no_par, fix_variables,
    identity_permutation, identity_permutation_mles, merge_polynomials, random_mle_list,
    random_permutation, random_permutation_mles, random_zero_mle_list, DenseMultilinearExtension,
};
pub use util::{bit_decompose, gen_eval_point, get_batched_nv, get_index};
pub use virtual_polynomial::{
    build_eq_x_r, build_eq_x_r_vec, eq_eval, VPAuxInfo, VirtualPolynomial,
};

extern crate alloc;

// ark-std v0.5 should re-export alloc/std::sync.
// While already released on crates.io, related versioning changes were reverted on the GitHub repo.
// Let's wait for this issue to stabilize.
#[cfg(target_has_atomic = "ptr")]
pub use alloc::sync::Arc as RefCounter;
#[cfg(not(target_has_atomic = "ptr"))]
pub use ark_std::rc::Rc as RefCounter;
