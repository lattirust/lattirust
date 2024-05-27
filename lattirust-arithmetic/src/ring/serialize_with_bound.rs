
pub struct WithBoundedBytes<T> {
    value: T,
    discard_trailing_bytes: usize,
}

// impl WithBoundedBytes<Z2_64> {
//     #[inline]
//     pub fn new(value: Z2_64, abs_max: i64) -> Self {
//         debug_assert!(value.0 .0.abs() as u64 <= abs_max);
//         Self {
//             value,
//             discard_trailing_bytes: abs_m,
//         }
//     }
// }
