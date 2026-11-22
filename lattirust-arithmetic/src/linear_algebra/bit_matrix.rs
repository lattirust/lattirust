//! Small packed-bit matrix operations.

/// Transpose eight little-endian bit rows. Output bit j of byte i equals
/// input bit i of byte j.
pub fn transpose8x8(rows: [u8; 8]) -> [u8; 8] {
    let mut word = u64::from_le_bytes(rows);
    for (shift, mask) in [
        (7, 0x00aa00aa00aa00aa),
        (14, 0x0000cccc0000cccc),
        (28, 0x00000000f0f0f0f0),
    ] {
        let swap = (word ^ (word >> shift)) & mask;
        word ^= swap ^ (swap << shift);
    }
    word.to_le_bytes()
}

/// Tabulate the sum of the entries selected by the eight-bit index.
pub fn subset_sums8<T: Copy + std::ops::Add<Output = T>>(weights: &[T; 8], zero: T) -> [T; 256] {
    let mut table = [zero; 256];
    for mask in 1usize..256 {
        let bit = mask.trailing_zeros() as usize;
        table[mask] = table[mask & (mask - 1)] + weights[bit];
    }
    table
}
