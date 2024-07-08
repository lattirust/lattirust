mod errors;
mod multilinear_polynomial;
mod univariate_polynomial;
mod util;
mod virtual_polynomial;

//TODO add exports
pub use errors::ArithErrors;
pub use util::{ bit_decompose, gen_eval_point, get_batched_nv, get_index };
