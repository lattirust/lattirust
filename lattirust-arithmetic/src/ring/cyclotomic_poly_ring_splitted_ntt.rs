use std::{fmt::Display, ops::Add};

use num_bigint::BigUint;

use crate::{linear_algebra::SVector, partial_ntt::PartialNTT, traits::Modulus};

use super::{Ring, Zq};

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
#[derive(PartialEq, Eq, Debug)]
pub struct CyclotomicPolyRingSplittedNTT<
    const Q: u64,
    const N: usize,
    const D: usize,
    const Z: usize,
    const PHI_Z: usize,
>(SVector<Zq<Q>, N>);

impl<const Q: u64, const N: usize, const D: usize, const Z: usize, const PHI_Z: usize> Display
    for CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.0, f)
    }
}

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
        let mut temp = [Zq::<Q>::from(0); N];
        let components = coprimes_set::<Z, PHI_Z>();
        let num_chunks = N / D;
        for k in 0..num_chunks {
            for i in 0..D {
                for j in 0..D - i {
                    temp[k * D + i + j] += self.0[i + k * D] * rhs.0[j + k * D];
                }
            }
            let rj_power = components[k] as u64;
            let rj = rou.pow([rj_power]);
            for i in 1..D {
                for j in D - i..D {
                    temp[k * D + i + j - D] += self.0[i + k * D] * rhs.0[j + k * D] * rj;
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

#[cfg(test)]
mod tests {
    use crate::ring::{Ring, Zq};

    use super::CyclotomicPolyRingSplittedNTT;

    const Q: u64 = 15 * (1 << 27) + 1;
    const N: usize = 1 << 4;
    const Z: usize = 1 << 2;
    const D: usize = 1 << 5 - 2;
    const PHI_Z: usize = 1 << 1;
    #[test]
    fn mul_test() {
        assert!(Q % Z as u64 == 1);
        let resize_power = Q - 1 / Z as u64;
        let rou = Zq::<Q>::from(76160998 as u64).pow([resize_power]);
        let rou3 = rou.pow([3]);
        let poly1 = CyclotomicPolyRingSplittedNTT::<Q, N, D, Z, PHI_Z>::from_fn(
            |i| Zq::<Q>::from(i as u64),
            rou,
        );
        let mut red_poly = [Zq::<Q>::from(0); N];
        for i in 0..N / 2 {
            red_poly[i] = Zq::<Q>::from(i as u64) + Zq::<Q>::from(8 + i as u64) * rou;
        }
        for i in 0..N / 2 {
            red_poly[N / 2 + i] = Zq::<Q>::from(i as u64) + Zq::<Q>::from(8 + i as u64) * rou3;
        }
        let poly_reduced = CyclotomicPolyRingSplittedNTT::<Q, N, D, Z, PHI_Z>::from_array(red_poly);
        assert_eq!(poly_reduced, poly1, "Poly NTT incorrect");
        let poly2 = CyclotomicPolyRingSplittedNTT::<Q, N, D, Z, PHI_Z>::from_fn(
            |i| Zq::<Q>::from((N - 1 - i) as u64),
            rou,
        );
        let final_poly = poly1.ntt_mul(&poly2, rou);

        let mut poly1_split1 = [Zq::<Q>::from(0); D];
        let mut poly2_split1 = [Zq::<Q>::from(0); D];
        let mut poly1_split2 = [Zq::<Q>::from(0); D];
        let mut poly2_split2 = [Zq::<Q>::from(0); D];
        for i in 0..D {
            poly1_split1[i] = poly1.0[i];
            poly2_split1[i] = poly2.0[i];

            poly1_split2[i] = poly1.0[i + D];
            poly2_split2[i] = poly2.0[i + D];
        }

        let naive_mul_poly_split1 = poly_mul(&poly1_split1, &poly2_split1);
        let naive_mul_poly_split2 = poly_mul(&poly1_split2, &poly2_split2);
        let split1 = reduce(naive_mul_poly_split1, rou);
        let split2 = reduce(naive_mul_poly_split2, rou3);
        let mut naive_coeffs = [Zq::<Q>::from(0); N];
        for i in 0..D {
            naive_coeffs[i] = split1[i];
            naive_coeffs[i + D] = split2[i];
        }

        let naive_poly =
            CyclotomicPolyRingSplittedNTT::<Q, N, D, Z, PHI_Z>::from_array(naive_coeffs);
        assert_eq!(
            naive_poly, final_poly,
            "Splitted NTT multiplication incorrect"
        );
    }

    fn poly_mul(a: &[Zq<Q>; D], b: &[Zq<Q>; D]) -> [Zq<Q>; 2 * D] {
        let mut temp = [Zq::<Q>::from(0); 2 * D];
        for i in 0..D {
            for j in 0..D {
                temp[i + j] += a[i] * b[j];
            }
        }
        temp
    }

    fn reduce(a: [Zq<Q>; 2 * D], rou: Zq<Q>) -> [Zq<Q>; D] {
        let mut temp = [Zq::<Q>::from(0); D];
        for i in 0..D {
            temp[i] += a[i];
        }
        for i in 0..D {
            temp[i] += a[i + D] * rou;
        }
        temp
    }
}
