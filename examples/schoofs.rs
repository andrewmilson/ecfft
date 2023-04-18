use ark_ff::One;
use ark_ff::PrimeField;
use ark_ff::Zero;
use ark_poly::univariate::DenseOrSparsePolynomial;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use ecfft::ec::Curve;
use ecfft::ec::Point;
use ecfft::m31::Fp;
use ecfft::utils::div_rem;
use ecfft::utils::gcd;
use ecfft::utils::pow_mod;
use ecfft::utils::xgcd;
use num_bigint::BigUint;
use num_integer::Integer;
use std::borrow::Borrow;
use std::collections::BTreeMap;
use std::ops::Range;

fn main() {
    let curve = Curve::new(Fp::one(), Fp::zero());
    let cardinality = schoofs_algorithm(&curve);
}

/// Returns the cardinality of a curve using Schoofs Algorithm
/// Implementation based on Algorithm 4.5 from "Elliptic Curves" book by LCW and
/// https://math.mit.edu/classes/18.783/2015/LectureNotes9.pdf.
pub fn schoofs_algorithm<F: PrimeField>(curve: &Curve<F>) -> BigUint {
    let Curve { a, b } = *curve;
    let x3_ax_b = DensePolynomial::from_coefficients_vec(vec![b, a, F::zero(), F::one()]);
    let p: BigUint = F::MODULUS.into();
    let mut congruences = BTreeMap::new();

    for l in schoof_primes::<F>() {
        if l == 2 {
            congruences.insert(l, if has_even_order(curve) { 0 } else { 1 });
            continue;
        }

        // handle odd primes
        let psi = division_polynomial(l, curve);
        congruences.insert(l, frobenius_trace_mod_l(l, &psi, curve));
    }

    // TODO: CRT

    todo!()
}

/// `modulus` can only be made up of factors of the l-th division polynomial
fn frobenius_trace_mod_l<F: PrimeField>(
    l: usize,
    modulus: &DensePolynomial<F>,
    curve: &Curve<F>,
) -> usize {
    assert!(l.is_odd(), "odd primes only");
    let p = F::MODULUS.into();
    let pl = usize::try_from(&p % l).unwrap();

    if modulus.degree() == 1 {
        // Modulus is only made up of factors of the l-th division polynomial.
        // If modulus is a degree 1 polynomial i.e. (x - c) for some constant c
        // then we've found a roots of the l-th division polynomial equal to c.
        // The roots of the l-th division polynomial correspond to the
        // x-coordinates of the points of the l-torsion group. Ok, so we found a
        // root and now know there is a point of order l. Therefore #E(F) = q + 1 + a =
        // 0 (mod l) and a = -q - 1 (mod l)
        return l - pl - 1;
    }

    let x = DensePolynomial::from_coefficients_vec(vec![F::zero(), F::one()]);
    let one = DensePolynomial::from_coefficients_vec(vec![F::one()]);
    let x3_ax_b =
        DensePolynomial::from_coefficients_vec(vec![curve.b, curve.a, F::zero(), F::one()]);

    // Compute π = (x^p, y^p) and π^2 = π ∘ π - note endomorphisms are multiplied by
    // composition. Also note y^p = (x^3 + A*x + B)^((p - 1) / 2) * y = b(x)*y
    // therefore π = (x^p, y^p) = (a(x), b(x)*y). Now we want to compute the
    // composition  of π_1 ∘ π_2. To do this, we apply π_1 to the input
    // coordinates (x, y), and then apply π_2 to the result. (1) Apply π_2:
    // (a_2(x), b_2(x)*y) (2) Apply π_1 to the result of π_2:
    //   * the x-coord transform: a_1(a_2(x))
    //   * the y-coord transform: b_1(x') * y' (where x' = a_2(x) and y' = b_2(x)*y)
    let pi_a = pow_mod(&x, p.clone(), modulus);
    let pi_b = pow_mod(&x3_ax_b, (p.clone() - 1u8) / 2u8, modulus);
    let pi_pi_a = pow_mod(&pi_a, p.clone(), modulus);
    let pi_b_pi_a = pow_mod(&pi_a, (p.clone() - 1u8) / 2u8, modulus);
    let pi_pi_b = div_rem(&pi_b_pi_a.naive_mul(&pi_b), modulus);
    let pi = (&pi_a, &pi_b);
    let pi_pi = (&pi_pi_a, &pi_pi_b);

    // (a_pl, b_pl * y) = pl * (x, y), pl = p (mod l)
    let (a_pl, b_pl) = match mul(pl.into(), (&x, &one), modulus, curve) {
        Err(UninvertablePolynomial(gcd)) => return frobenius_trace_mod_l(l, &gcd, curve),
        Ok(res) => res,
    };

    // (x', y') = (a', b' * y) = π^2 + (p mod l) * (x, y) (mod modulus)
    let (a_prime, b_prime) = match add(pi_pi, (&a_pl, &b_pl), modulus, curve) {
        Err(UninvertablePolynomial(gcd)) => return frobenius_trace_mod_l(l, &gcd, curve),
        Ok(res) => res,
    };

    for j in 1..l / 2 {
        // (a_j, b_j * y) = j * (x, y)
        let (a_j, b_j) = match mul(j.into(), (&x, &one), modulus, curve) {
            Err(UninvertablePolynomial(gcd)) => return frobenius_trace_mod_l(l, &gcd, curve),
            Ok(res) => res,
        };

        // (a_j)^p
        let pi_a_j = pow_mod(&a_j, p.clone(), modulus);

        // check if x' − x_j^p ≡ 0
        if (&a_prime - &pi_a_j).is_zero() {
            // compute (y' − y_j^p)/y = b' - b_j^((p - 1) / 2)
            let term1 = pow_mod(&b_j, (p - 1u8) / 2u8, modulus);

            // check if (y' − y_j^p)/y ≡ 0
            if (&b_prime - &term1).is_zero() {
                return j;
            } else {
                return l - j;
            }
        }
    }

    // All values 1 ≤ j ≤ (l − 1)/2 have been tried without success.
    // Now let ω^2 ≡ p (mod l). If ω does not exist, then a ≡ 0 (mod l).
    for omega in 0..l {
        if (omega * omega) % l == pl {}
    }

    // ω does not exist so a = 0
    0
}

// /// Computes the scalar multiplication map [n]
// /// Return value is of the form (r1, r2) where n(x, y) = (r1, r2 * y)
// fn mul<F: PrimeField>(
//     n: usize,
//     r1: &DensePolynomial<F>,
//     r2: &DensePolynomial<F>,
//     modulus: &DensePolynomial<F>,
// ) -> DensePolynomial<F> {
//     todo!()
// }

/// Computes scalar multiplication n * p using double-and-add
/// Note: p = (a(x), b(x) * y)
fn mul<F: PrimeField>(
    mut n: BigUint,
    (a, b): (&DensePolynomial<F>, &DensePolynomial<F>),
    modulus: &DensePolynomial<F>,
    curve: &Curve<F>,
) -> Result<(DensePolynomial<F>, DensePolynomial<F>), UninvertablePolynomial<F>> {
    let one = DensePolynomial::from_coefficients_slice(&[F::one()]);
    let (mut a_prime, mut b_prime) = (one.clone(), one);
    let (mut a_acc, mut b_acc) = (a.clone(), b.clone());
    while !n.is_zero() {
        if n.is_odd() {
            (a_prime, b_prime) = add((&a_prime, &b_prime), (&a_acc, &b_acc), modulus, curve)?;
        }
        (a_acc, b_acc) = add((&a_acc, &b_acc), (&a_acc, &b_acc), modulus, curve)?;
        n >>= 1;
    }
    Ok((a_prime, b_prime))
}

/// Computes α_1 + α_2 = α_3
/// Note: α_i = (a_i(x), b_i(x) * y)
fn add<F: PrimeField>(
    (a1, b1): (&DensePolynomial<F>, &DensePolynomial<F>),
    (a2, b2): (&DensePolynomial<F>, &DensePolynomial<F>),
    modulus: &DensePolynomial<F>,
    curve: &Curve<F>,
) -> Result<(DensePolynomial<F>, DensePolynomial<F>), UninvertablePolynomial<F>> {
    let one = DensePolynomial::from_coefficients_vec(vec![F::one()]);
    let x3_ax_b =
        DensePolynomial::from_coefficients_vec(vec![curve.b, curve.a, F::zero(), F::one()]);

    // m = y*(b_1 - b_2)/(a_1 - a_2) if a_1 != a_2
    // m = (3*a_1^2 + A)/(2*b_1*y)   if a_1 == a_2
    // m = r(x) * y
    let r = if a1 == a2 {
        // calculate tangent line
        let numerator = a1 * F::from(3u8) + DensePolynomial::from_coefficients_vec(vec![curve.a]);
        let denominator = div_rem(&(&x3_ax_b.naive_mul(b1) * F::from(2u8)), modulus);
        let (denominator_inv, _, gcd) = xgcd(&denominator, modulus);
        if gcd != one {
            return Err(UninvertablePolynomial(gcd));
        }
        &numerator * &denominator_inv
    } else {
        // calculate slope
        let numerator = b1 - b2;
        let denominator = a1 - a2;
        let (denominator_inv, _, gcd) = xgcd(&denominator, modulus);
        if gcd != one {
            return Err(UninvertablePolynomial(gcd));
        }
        &numerator * &denominator_inv
    };

    // note that y^2 = x^3 + A*x + B
    // a_3 = r^2 * y^2 - a_1 - a_2
    // b_3 = r * (a_1 - a_3) - b_1
    let rr = pow_mod(&r, 2u8.into(), modulus);
    let a3 = &div_rem(&rr.naive_mul(&x3_ax_b), modulus) - &(a1 - a2);
    let b3 = &div_rem(&r.naive_mul(&(a1 - &a3)), modulus) - b1;
    Ok((a3, b3))
}

/// Holds the GCD of the polynomial and the modulus
struct UninvertablePolynomial<F: PrimeField>(DensePolynomial<F>);

/// Returns a list of primes such that primes[0] * ... * primes[n-1] > 4*sqrt(p)
fn schoof_primes<F: PrimeField>() -> Vec<usize> {
    let hi = hasse_interval::<F>();
    let hasse_interval_len = hi.end - hi.start;
    let mut res = Vec::new();
    let mut product = BigUint::one();
    let mut small_primes = SMALL_PRIMES.into_iter();
    while product < hasse_interval_len {
        let prime = small_primes.next().unwrap();
        res.push(prime as usize);
        product *= prime;
    }
    res
}

/// Returns true if there is an even number of points on the curve.
fn has_even_order<F: PrimeField>(curve: &Curve<F>) -> bool {
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
    let Curve { a, b } = *curve;
    let p: BigUint = F::MODULUS.into();

    // Compute xp ≡ x^p (mod x^3 + A*x + B) by successive squaring
    let x3_ax_b = DensePolynomial::from_coefficients_vec(vec![b, a, F::zero(), F::one()]);
    let x = DensePolynomial::from_coefficients_vec(vec![F::zero(), F::one()]);
    let xp = pow_mod(&x, p, &x3_ax_b);

    // Compute gcd(xp - x, x^3 + A*x + B). If the gcd is 1, then there's no root and
    // the order is odd.
    gcd(&(&xp - &x), &x3_ax_b).degree() != 0
}

fn hasse_interval<F: PrimeField>() -> Range<BigUint> {
    let p: BigUint = F::MODULUS.into();
    (&p + 1u32 - p.sqrt() * 2u32)..(&p + 1u32 + p.sqrt() * 2u32)
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
            -curve.a.pow([3]) * F::from(2u8) - curve.b * F::from(16u8),
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
            let Curve { a, b } = *curve;
            let yy = DensePolynomial::from_coefficients_vec(vec![b, a, F::zero(), F::one()]);
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

/// Determines a points order `n` in `O(n)`
fn order_naive<F: PrimeField>(p: Point<F>) -> usize {
    let mut acc = Point::zero();
    let mut order = 0;
    while !acc.is_zero() {
        acc += p;
        order += 1;
    }
    order
}

// first 46 primes
const SMALL_PRIMES: [u32; 46] = [
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
