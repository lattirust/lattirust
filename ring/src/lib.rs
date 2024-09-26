#![allow(incomplete_features)]
#![allow(non_snake_case)]
#![feature(generic_const_exprs)]
#![feature(vec_into_raw_parts)]
// Exports
pub use poly_ring::*;
pub use representatives::{SignedRepresentative, UnsignedRepresentative};
pub use ring::*;

pub mod balanced_decomposition;
pub mod cyclotomic_ring;
pub mod traits;
pub mod zn;

mod poly_ring;
mod representatives;
mod ring;

extern crate core;
