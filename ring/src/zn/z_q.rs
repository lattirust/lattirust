use ark_ff::{BigInt, Fp, Fp64, MontBackend, MontConfig};

pub struct FqConfig<const Q: u64> {}

impl<const Q: u64> MontConfig<1> for FqConfig<Q> {
    const MODULUS: BigInt<1> = BigInt::<1>([Q]);
    // This works for primes of the form 2^x + 1
    const GENERATOR: Fp<MontBackend<Self, 1>, 1> = Fp::new(BigInt([3]));
    // This works for primes of the form 2^x + 1
    const TWO_ADIC_ROOT_OF_UNITY: Fp<MontBackend<Self, 1>, 1> = Fp::new(BigInt([3])); // TODO: implement this? Not using todo!() here to generate the docs
}

pub type Zq<const Q: u64> = Fp64<MontBackend<FqConfig<Q>, 1>>;

pub const fn const_fq_from<const Q: u64>(val: u64) -> Zq<Q> {
    Zq::new(BigInt::<1>([val]))
}
