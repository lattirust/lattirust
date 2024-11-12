// #![allow(dead_code)]
// #![allow(non_snake_case)]
// #![allow(unused_imports)]

// use cxx::{CxxVector, let_cxx_string};
// pub use cxx;

// #[cxx::bridge]
// pub mod ffi {
//     // Shared structs
//     struct GaussianParameters {
//         std_dev: f64,
//         base: i32,
//         smoothing_param: f64,
//     }

//     unsafe extern "C++" {
//         include!("openfhe/discretegaussiangenerator-impl.h");
//         include!("openfhe/discretegaussiangenerator.h");
//         include!("openfhe/discretegaussiangeneratorgeneric.h");
        
//         // Wrap DiscreteGaussianGenerator
//         type DiscreteGaussianGenerator;
        
//         // Constructor
//         fn new_gaussian_sampler(params: &GaussianParameters) -> UniquePtr<DiscreteGaussianGenerator>;
        
//         // Sample method
//         fn generate_sample(self: &DiscreteGaussianGenerator, center: f64) -> i64;
//     }
// }
