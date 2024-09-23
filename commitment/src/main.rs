use commitments::ppk::get_gaussian_vec;
use rand::SeedableRng;

fn main() {
    // TODO: not suitable for crypto use, fix it
    let mut rng = rand::rngs::StdRng::seed_from_u64(20);
    let t = 12;         // Plaintext modulus
    let q = 65536;      // Ciphertext modulus
    let std_dev = 1.0;  // Standard deviation for generating the error
    let degree = 4;     // Degree of polynomials used for encoding and encrypting messages
    let rlk_base = (q as f64).log2() as i64; // The base for decomposition during relinearization

    let v = get_gaussian_vec(std_dev, 50, &mut rng);
    println!("{:?}", v);
}