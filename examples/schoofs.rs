use ark_ff::One;
use ark_ff::PrimeField;
use ark_ff::Zero;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use ecfft::ec::Curve;
use ecfft::m31::Fp;
// use ecfft::secp256k1::Fp;
use ecfft::utils::div_rem;
use ecfft::utils::gcd;
use ecfft::utils::pow_mod;
use ecfft::utils::xgcd;
use num_bigint::BigInt;
use num_bigint::BigUint;
use num_integer::ExtendedGcd;
use num_integer::Integer;
use std::ops::Add;
use std::ops::Mul;

fn main() {
    let curve = Curve::new(Fp::from(12328), Fp::from(273812783));
    println!("Cardinality is: {}", cardinality(curve));
}

/// Returns the cardinality of a curve using Schoofs Algorithm
/// Implementation based on Algorithm 4.5 from "Elliptic Curves" book by LCW and
/// https://math.mit.edu/classes/18.783/2015/LectureNotes9.pdf.
pub fn cardinality<F: PrimeField>(curve: Curve<F>) -> BigUint {
    let p = F::MODULUS.into();
    let hasse_interval_len = BigInt::from(p.sqrt() * 4u32);
    let mut small_primes = SMALL_PRIMES.into_iter();
    let mut a_mod_m = BigInt::zero();
    let mut m = BigInt::one();
    while m < hasse_interval_len {
        let l = small_primes.next().unwrap();
        let a_mod_l = if l == 2 {
            // p + 1 + a ≡ 0 (mod 2) then a = 0 (mod 2) otherwise a = 1 (mod 2)
            if has_even_order(curve) {
                0
            } else {
                1
            }
        } else {
            // handle odd primes using Schoofs algorithm
            let psi = division_polynomial(l, &curve);
            frobenius_trace_mod_l(l, curve, &psi.into())
        };

        // Incremental Chinese Remainder Theorem. From algorithm 9.1 in:
        // https://math.mit.edu/classes/18.783/2015/LectureNotes9.pdf
        // calculates the frobenius trace mod `m` where m = l_0 * l_1 * ... * l_n
        let ExtendedGcd {
            x: m_inv_l,
            y: l_inv_m,
            ..
        } = m.extended_gcd(&l.into());
        a_mod_m = (m_inv_l * &m * a_mod_l + l_inv_m * l * a_mod_m) % (&m * l);
        m *= l;
    }

    let a_mod_m = BigUint::try_from((&m + a_mod_m) % &m).unwrap();
    let m = BigUint::try_from(m).unwrap();
    if a_mod_m > &m / 2u8 {
        p + 1u8 + a_mod_m - m
    } else {
        p + 1u8 + a_mod_m
    }
}

/// Calculates the trace of the frobenius mod `l` using Schoof's algorithm
/// `modulus` can only be made up of factors of the l-th division polynomial
/// Panics if `l` is not an odd prime or if modulus has degree 0
fn frobenius_trace_mod_l<F: PrimeField>(
    l: usize,
    curve: Curve<F>,
    ring: &QuotientRing<F>,
) -> usize {
    assert!(l.is_odd(), "odd primes only");
    assert!(!ring.modulus.degree().is_zero());
    let p = F::MODULUS.into();
    let p_mod_l = usize::try_from(&p % l).unwrap();

    if ring.modulus.degree() == 1 {
        // Modulus is only made up of factors of the l-th division polynomial.
        // If modulus is a degree 1 polynomial i.e. (x - c) for some constant c
        // then we've found a roots of the l-th division polynomial equal to c.
        // The roots of the l-th division polynomial correspond to the
        // x-coordinates of the points of the l-torsion group. Ok, so we found a
        // root and now know there is a point of order l. Therefore #E(F) = q + 1 + a =
        // 0 (mod l) and a = -q - 1 (mod l)
        return l - p_mod_l - 1;
    }

    let x = DensePolynomial::from_coefficients_vec(vec![F::zero(), F::one()]);
    let x3_ax_b = curve.x3_ax_b();

    // Compute π = (x^p, y^p) and π^2 = π ∘ π - note endomorphisms are multiplied by
    // composition. Also note y^p = (x^3 + A*x + B)^((p - 1) / 2) * y = b(x)*y
    // therefore π = (x^p, y^p) = (a(x), b(x)*y). Now we want to compute the
    // composition  of π_1 ∘ π_2. To do this, we apply π_1 to the input
    // coordinates (x, y), and then apply π_2 to the result. (1) Apply π_2:
    // (a_2(x), b_2(x)*y) (2) Apply π_1 to the result of π_2:
    // * the x-coord transform: a_1(a_2(x))
    // * the y-coord transform: b_1(x') * y' (where x' = a_2(x) and y' = b_2(x)*y)
    let pi_a = ring.pow(&x, p.clone());
    let pi_b = ring.pow(&x3_ax_b, (p - 1u8) / 2u8);
    let pi = Endomorphism::new(pi_a, pi_b, curve, ring);
    let pi_pi = pi.clone() * pi.clone();

    // (a_pl, b_pl * y) = (p mod l) * (x, y)
    let pl = match Endomorphism::identity(curve, ring) * BigUint::from(p_mod_l) {
        Err(Uninvertable(gcd)) => return frobenius_trace_mod_l(l, curve, &gcd.into()),
        Ok(res) => res,
    };

    // (x', y') = (a', b' * y) = π^2 + pl
    let e_prime = match pi_pi + pl {
        Err(Uninvertable(gcd)) => return frobenius_trace_mod_l(l, curve, &gcd.into()),
        Ok(res) => res,
    };

    for j in 1..l {
        // (a_j, b_j * y) = j * π
        let pi_j = match pi.clone() * BigUint::from(j) {
            Err(Uninvertable(gcd)) => return frobenius_trace_mod_l(l, curve, &gcd.into()),
            Ok(res) => res,
        };

        if pi_j == e_prime {
            return j;
        }
    }

    0
}

/// Stores an elliptic curve endomorphism in a quotient ring
/// (x, y) = (a(x), b(x) * y)
#[derive(Debug, Clone, PartialEq, Eq)]
struct Endomorphism<'a, F: PrimeField> {
    a: DensePolynomial<F>,
    b: DensePolynomial<F>,
    curve: Curve<F>,
    ring: &'a QuotientRing<F>,
}

impl<'a, F: PrimeField> Endomorphism<'a, F> {
    fn new(
        a: DensePolynomial<F>,
        b: DensePolynomial<F>,
        curve: Curve<F>,
        ring: &'a QuotientRing<F>,
    ) -> Self {
        Self {
            a: ring.reduce_poly(a),
            b: ring.reduce_poly(b),
            curve,
            ring,
        }
    }

    fn identity(curve: Curve<F>, ring: &'a QuotientRing<F>) -> Self {
        // TODO: identity endomorphism
        Self {
            a: DensePolynomial::from_coefficients_vec(vec![F::zero(), F::one()]),
            b: DensePolynomial::from_coefficients_vec(vec![F::one()]),
            curve,
            ring,
        }
    }
}

impl<F: PrimeField> Mul<BigUint> for Endomorphism<'_, F> {
    type Output = Result<Self, Uninvertable<F>>;

    /// Computes scalar multiplication of an endomorphism using double-and-add
    fn mul(self, mut n: BigUint) -> Self::Output {
        let mut res = None;
        let mut acc = self;
        while !n.is_zero() {
            if n.is_odd() {
                res = Some(match res.clone() {
                    Some(res) => (res + acc.clone())?,
                    None => acc.clone(),
                });
            }
            acc = (acc.clone() + acc)?;
            n >>= 1;
        }
        Ok(res.unwrap())
    }
}

impl<F: PrimeField> Mul<Self> for Endomorphism<'_, F> {
    type Output = Self;

    /// Multiplies two endomorphisms by composition
    fn mul(self, rhs: Self) -> Self {
        assert_eq!(self.ring, rhs.ring);
        assert_eq!(self.curve, rhs.curve);
        let ring = self.ring;
        let curve = self.curve;
        let Self { a: a1, b: b1, .. } = self;
        let Self { a: a2, b: b2, .. } = rhs;
        // compose a_1 with a_2 to obtain a_3 = (a_1(a_2(x)), b_1(a_2(x)) * b_2(x) * y)
        let mut a1_a2 = DensePolynomial::zero();
        let mut b1_a2 = DensePolynomial::zero();
        let mut acc = QuotientRing::one();
        for i in 0..a1.len().max(b1.len()) {
            if let Some(&coeff) = a1.get(i) {
                a1_a2 += &(&acc * coeff);
            }
            if let Some(&coeff) = b1.get(i) {
                b1_a2 += &(&acc * coeff);
            }
            acc = ring.mul(&acc, &a2);
        }
        // TODO: remove https://github.com/arkworks-rs/algebra/pull/638
        while a1_a2.coeffs.last().map_or(false, |c| c.is_zero()) {
            a1_a2.coeffs.pop();
        }
        while b1_a2.coeffs.last().map_or(false, |c| c.is_zero()) {
            b1_a2.coeffs.pop();
        }
        Self {
            a: a1_a2,
            b: ring.mul(&b1_a2, &b2),
            ring,
            curve,
        }
    }
}

impl<F: PrimeField> Add for Endomorphism<'_, F> {
    type Output = Result<Self, Uninvertable<F>>;

    /// Applies elliptic curve addition operation to endomorphisms
    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.ring, rhs.ring);
        assert_eq!(self.curve, rhs.curve);
        let ring = self.ring;
        let curve = self.curve;
        let x3_ax_b = curve.x3_ax_b();

        let Self { a: a1, b: b1, .. } = self;
        let Self { a: a2, b: b2, .. } = rhs;

        // m = y*(b_1 - b_2)/(a_1 - a_2) if a_1 != a_2
        // m = (3*a_1^2 + A)/(2*b_1*y)   if a_1 == a_2
        // m = r(x) * y
        let r = if a1 == a2 {
            // calculate tangent line
            let a = DensePolynomial::from_coefficients_vec(vec![curve.a]);
            let numerator = &(&ring.mul(&a1, &a1) * F::from(3u8)) + &a;
            let denominator = &x3_ax_b.naive_mul(&b1) * F::from(2u8);
            ring.div(&numerator, &denominator)?
        } else {
            // calculate slope
            ring.div(&(&b1 - &b2), &(&a1 - &a2))?
        };

        // note that y^2 = x^3 + A*x + B
        // a_3 = r^2 * y^2 - a_1 - a_2
        // b_3 = r * (a_1 - a_3) - b_1
        let rr = ring.pow(&r, 2u8.into());
        let a3 = &ring.mul(&rr, &x3_ax_b) - &(&a1 + &a2);
        let b3 = &ring.mul(&r, &(&a1 - &a3)) - &b1;
        Ok(Self::new(a3, b3, curve, ring))
    }
}

/// Quotient ring with modulus f(x)
#[derive(Debug, Clone, PartialEq, Eq)]
struct QuotientRing<F: PrimeField> {
    modulus: DensePolynomial<F>,
}

impl<F: PrimeField> QuotientRing<F> {
    fn new(modulus: DensePolynomial<F>) -> Self {
        assert!(!modulus.is_zero());
        Self { modulus }
    }

    fn one() -> DensePolynomial<F> {
        DensePolynomial::from_coefficients_vec(vec![F::one()])
    }

    fn mul(&self, a: &DensePolynomial<F>, b: &DensePolynomial<F>) -> DensePolynomial<F> {
        self.reduce_poly(a.naive_mul(b))
    }

    fn pow(&self, a: &DensePolynomial<F>, exp: BigUint) -> DensePolynomial<F> {
        if self.modulus.is_zero() {
            unimplemented!()
        }

        pow_mod(a, exp, &self.modulus)
    }

    /// If the denominator is invertible, the method returns a `Result`
    /// containing the division polynomial. Otherwise, it returns an error
    /// indicating that the denominator is not invertible.
    fn div(
        &self,
        numerator: &DensePolynomial<F>,
        denominator: &DensePolynomial<F>,
    ) -> Result<DensePolynomial<F>, Uninvertable<F>> {
        let denominator_inv = self.try_inverse(denominator)?;
        Ok(self.reduce_poly(numerator.naive_mul(&denominator_inv)))
    }

    /// Reduces a polynomial into the quotient ring, returning its canonical
    /// representative.
    fn reduce_poly(&self, a: DensePolynomial<F>) -> DensePolynomial<F> {
        if self.modulus.is_zero() {
            a
        } else {
            div_rem(&a, &self.modulus)
        }
    }

    fn try_inverse(&self, a: &DensePolynomial<F>) -> Result<DensePolynomial<F>, Uninvertable<F>> {
        let (a_inverse, _, gcd) = xgcd(a, &self.modulus);
        if gcd == Self::one() {
            Ok(a_inverse)
        } else {
            Err(Uninvertable(gcd))
        }
    }
}

impl<F: PrimeField> From<DensePolynomial<F>> for QuotientRing<F> {
    fn from(modulus: DensePolynomial<F>) -> Self {
        Self::new(modulus)
    }
}

/// Holds the GCD of the polynomial and the modulus
struct Uninvertable<F: PrimeField>(DensePolynomial<F>);

/// Returns true if there is an even number of points on the curve.
fn has_even_order<F: PrimeField>(curve: Curve<F>) -> bool {
    // "if x^3 + A*x + B has a root e ∈ Fp, then (e, 0) ∈ E[2] and (e, 0) ∈ E(Fp),
    // so E(Fp) has even order ... if x^3 + A*x + B has no roots in Fp, then E(Fp)
    // has no points of order 2" - Elliptic Curves, LCW.
    //
    // This can be determined by trying all values of x in the field but there is a
    // faster way. Note that x^p − x = (x - 1)(x - 2)...(x - (p - 1)) therefore
    // if gcd(x^p − x, x^3 + A*x + B) != 1 then a root exists and therefore
    // there is a point with order 2.
    //
    // Since "p" can be large we first compute xp ≡ x^p (mod x^3 + A*x + B) by
    // successive squaring. Then we just need to compute gcd(xp − x, x^3 + Ax + B).
    let p: BigUint = F::MODULUS.into();
    let x3_ax_b = curve.x3_ax_b();
    let one = DensePolynomial::from_coefficients_vec(vec![F::one()]);
    let x = DensePolynomial::from_coefficients_vec(vec![F::zero(), F::one()]);
    let xp = pow_mod(&x, p, &x3_ax_b);

    // Compute gcd(xp - x, x^3 + A*x + B). If the gcd is 1, then there's no root and
    // the order is odd.
    gcd(&(&xp - &x), &x3_ax_b) != one
}

/// Returns the nth division polynomial
/// The division polynomial for even values of n are missing a factor of 2y
pub fn division_polynomial<F: PrimeField>(n: usize, curve: &Curve<F>) -> DensePolynomial<F> {
    if n == 0 {
        // ψ_0 = 0
        DensePolynomial::default()
    } else if n == 1 || n == 2 {
        // ψ_1 = 1; ψ_2 = 2*y (2*y factor removed see above)
        DensePolynomial::from_coefficients_vec(vec![F::one()])
    } else if n == 3 {
        // ψ_3 = 3*x^4 + 6*A*x^2 + 12*B*x - A^2
        DensePolynomial::from_coefficients_vec(vec![
            -curve.a.square(),
            curve.b * F::from(12u8),
            curve.a * F::from(6u8),
            F::zero(),
            F::from(3u8),
        ])
    } else if n == 4 {
        // ψ_4 = 4*y*(x^6 + 5*A*x^4 + 20*B*x^3 − 5*A*A*x^2 − 4*A*B*x − 8*B^2 − A^3)
        //     = 2*y*(2*x^6 + 10*A*x^4 + 40*B*x^3 − 10*AA*x^2 − 8*AB*x − 16*B^2 − 2*A^3)
        DensePolynomial::from_coefficients_vec(vec![
            -curve.a.pow([3]) * F::from(2u8) - curve.b.square() * F::from(16u8),
            -curve.a * curve.b * F::from(8u8),
            -curve.a.square() * F::from(10u8),
            curve.b * F::from(40u8),
            curve.a * F::from(10u8),
            F::zero(),
            F::from(2u32),
        ])
    } else {
        let m = n / 2;
        let psi_m = division_polynomial(m, curve);
        let psi_ms1 = division_polynomial(m - 1, curve);
        let psi_mp1 = division_polynomial(m + 1, curve);
        let psi_mp2 = division_polynomial(m + 2, curve);
        if n.is_even() {
            // t0 = ψ_(m+2) * ψ_(m-1)^2
            // t1 = ψ_(m-2) * ψ_(m+1)^2
            let psi_ms2 = division_polynomial(m - 2, curve);
            let t0 = psi_ms1.naive_mul(&psi_ms1).naive_mul(&psi_mp2);
            let t1 = psi_mp1.naive_mul(&psi_mp1).naive_mul(&psi_ms2);
            (&t0 - &t1).naive_mul(&psi_m)
        } else {
            // = (2 * y)^4
            // = 16 * (y^2)^2
            // = 16 * (x^3 + A*x + B)^2
            let yy = curve.x3_ax_b();
            let yyyy16 = &yy.naive_mul(&yy) * F::from(16u8);
            // t0 = ψ_(m+2) * ψ_m^3
            // t1 = ψ_(m-1) * ψ_(m+1)^3
            let t0 = psi_mp2.naive_mul(&psi_m.naive_mul(&psi_m).naive_mul(&psi_m));
            let t1 = psi_ms1.naive_mul(&psi_mp1.naive_mul(&psi_mp1).naive_mul(&psi_mp1));
            if m.is_even() {
                &yyyy16.naive_mul(&t0) - &t1
            } else {
                &t0 - &yyyy16.naive_mul(&t1)
            }
        }
    }
}

// first 46 primes
const SMALL_PRIMES: [usize; 46] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199,
];

#[cfg(test)]
mod tests {
    use crate::ec::Curve;
    use crate::ecfft::FFTree;
    use crate::secp256k1::Fp;
    use crate::utils;
    use ark_ff::Fp256;
    use ark_ff::MontBackend;
    use ark_ff::MontConfig;
    use ark_ff::One;
    use ark_ff::PrimeField;
    use ark_ff::UniformRand;
    use ark_ff::Zero;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;

    #[test]
    fn roots_test() {
        let mut rng = rand::thread_rng();
        let a = Fp::rand(&mut rng);
        let b = Fp::rand(&mut rng);
        let x3_ax_b = DensePolynomial::from_coefficients_vec(vec![b, a, Fp::zero(), Fp::one()]);
        let curve = Curve::new(a, b);
        // println!("Roots: {:?}", utils::find_roots(&x3_ax_b));
    }
}
