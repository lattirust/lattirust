use std::ops::{Add, Mul};

use num_bigint::BigUint;

use crate::{linear_algebra::SVector, partial_ntt::PartialNTT, traits::Modulus};

use super::{Zq, Ring};

/// For an m-th cyclotomic ring its minimal polynomail has the form
/// `Phi_m = Prod_{j=1}^phi(m)(x - r_j)`
/// Where `ord_p(r_j) = m` and `r_j = w^i` where `i \in {i: gcd(i,m) = 1}`
/// and `j` is its enumaration
/// Note that this requires that `p = 1 mod m`
///
/// Another option is to split the minimal polynomial as
/// `Phi_m = Prod_{j = 1}^phi(z) (x^{m/z} - r_j)`
/// Where `ord_p(r_j) = z` and `r_j = w^i` where `i \in {i: gcd(i,z) = 1}`
/// this requires that `p = q mod z` so `m` doesn't have to divide `p-1` but
/// z has to share the prime factorizations with diferent powers
/// For more see:
///     Short invertible elements in partially splitting cyclotomic rings
///
/// We have that:
/// R = DirProd R_j and R_j \in Z_p[X]/ (X^{m/z} - r_j)
/// This struct can be thought as the concatenation of R_j elements
/// where the operations over each R_j are element-wise
/// let `D` = m/z so we have phi(z) R_j components each with D components
/// also phi(z) = N/D
pub struct CyclotomicPolyRingSplittedNTT<
    const Q: u64,
    const N: usize,
    const D: usize,
    const Z: usize,
    const PHI_Z: usize,
>(SVector<Zq<Q>, N>);

impl<const Q: u64, const N: usize, const D: usize, const Z: usize, const PHI_Z: usize>
    CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>
{
    #[allow(dead_code)]
    pub(crate) type Inner = SVector<Zq<Q>, N>;

    /// Constructs a polynomial from an array of coefficients in NTT form.
    /// This function is private since we can't enforce coeffs being in NTT form if called from outside.
    fn from_array(coeffs_ntt: [Zq<Q>; N]) -> Self {
        Self(Self::Inner::const_from_array(coeffs_ntt))
    }

    pub fn from_fn<F>(f: F, rou: Zq<Q>) -> Self
    where
        F: FnMut(usize) -> Zq<Q>,
    {
        let mut coeffs = core::array::from_fn(f);
        Self::ntt(&mut coeffs, rou);
        Self::from_array(coeffs)
    }

    pub fn ntt_mul(&self, rhs: &Self, rou: Zq<Q>) -> Self {
        let binding_lhs = self.0.iter().collect::<Vec<_>>();
        let concat_lhs = binding_lhs.chunks_exact(D);
        let binding_rhs = rhs.0.iter().collect::<Vec<_>>();
        let concat_rhs = binding_rhs.chunks_exact(D);
        let mut temp = [Zq::<Q>::from(0); N];
        let components = coprimes_set::<Z, PHI_Z>();
        for (k, (lhs, rhs)) in concat_lhs.zip(concat_rhs).enumerate() {
            // This is the multiplication of factors that doesn't need to be reduced
            for mu in 0..D {
                for nu in 0..D - mu {
                    temp[k * D + mu + nu] = lhs[mu] * rhs[nu];
                }
            }
            let rj_power = components[k] as u64;
            let rj = rou.pow([rj_power]);
            for mu in 1..D {
                for nu in D - mu..D {
                    temp[k * D + mu + nu - D] = lhs[mu] * rhs[nu] * rj;
                }
            }
        }
        Self::from_array(temp)
    }
}

// TODO: impl SplittedNTT
impl<const Q: u64, const N: usize, const D: usize, const Z: usize, const PHI_Z: usize>
    PartialNTT<Q, N, D, Z, PHI_Z> for CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>
{
}

impl<const Q: u64, const N: usize, const D: usize, const Z: usize, const PHI_Z: usize> Modulus
    for CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>
{
    fn modulus() -> BigUint {
        Zq::<Q>::modulus()
    }
}

const fn vec_from_element<const Q: u64, const N: usize>(elem: Zq<Q>) -> SVector<Zq<Q>, N> {
    SVector::<Zq<Q>, N>::const_from_array([elem; N])
}

impl<const Q: u64, const N: usize, const D: usize, const Z: usize, const PHI_Z: usize> Add
    for CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let mut res = [Zq::<Q>::from(0); N];
        for i in 0..N {
            res[i] = self.0[i] + rhs.0[i];
        }
        Self::from_array(res)
    }
}

const fn coprimes_set<const Z: usize, const PHI_Z: usize>() -> [usize; PHI_Z] {
    let mut i = 0;
    let mut j = 0;
    let mut set = [0; PHI_Z];
    while i < Z {
        if gcd::<Z>(i) {
            set[j] = i;
            j += 1;
        }
        i += 1;
    }
    set
}

const fn gcd<const Z: usize>(i: usize) -> bool {
    let mut a = i;
    let mut b = Z;
    while b != 0 {
        let temp = b;
        b = a % b;
        a = temp
    }
    a == 1
}
