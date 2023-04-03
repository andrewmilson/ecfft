use crate::ec::Curve;
use crate::utils;
use ark_ff::One;
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use num_bigint::BigUint;
use num_integer::Integer;
use std::collections::BTreeMap;
use std::ops::Range;

/// Returns the cardinality of a curve using Schoofs Algorithm
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
        let pl = &p % l;
        let psi = division_polynomial(l, curve);
        todo!()
    }

    todo!()
}

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
    let xp = utils::pow_mod(&x, p, &x3_ax_b);

    // Compute gcd(xp - x, x^3 + A*x + B). If the gcd is 1, then there's no root and
    // the order is odd.
    utils::gcd(&(&xp - &x), &x3_ax_b).degree() != 0
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
