/// Return x^k mod Q
pub(crate) const fn const_pow_mod<const Q: u64>(x: u64, k: u64) -> u64 {
    let mut res: u128 = 1;
    let mut x: u128 = x as u128;
    let mut k: u128 = k as u128;
    while k > 0 {
        if k % 2 == 1 {
            res = (res * x) % (Q as u128);
        }
        x = (x * x) % (Q as u128);
        k /= 2;
    }
    res as u64
}
