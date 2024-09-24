#![feature(const_trait_impl)]
#![allow(incomplete_features)]
#![feature(inherent_associated_types)]
#![feature(int_roundings)]
#![feature(const_option)]
#![allow(non_snake_case)]
#![feature(generic_const_exprs)]
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
