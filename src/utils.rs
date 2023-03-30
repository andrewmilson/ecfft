use ark_ff::PrimeField;
use ark_ff::Zero;
use ark_poly::univariate::DenseOrSparsePolynomial;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use num_bigint::BigUint;
use num_integer::Integer;
use rand::Rng;
use std::collections::BTreeMap;

type Degree = usize;

/// Finds all unique roots in the field F
pub fn find_roots<F: PrimeField>(poly: &DensePolynomial<F>) -> Vec<F> {
    let f = square_free_factors(poly);
    let ddf = distinct_degree_factors(&f);
    let degree_1_factors = equal_degree_factorization(ddf.get(&1).unwrap(), 1);
    let mut roots: Vec<F> = degree_1_factors
        .into_iter()
        .map(|factor| {
            // factor = x + c
            // root = -c
            assert_eq!(1, factor.degree());
            -factor[0]
        })
        .collect();
    roots.sort();
    roots
}

/// Returns a mapping of degrees d to the factors of the polynomial, f,
/// such that each factor is the product of all monic irreducible factors
/// of f of degree d.
///
/// The input must be a squarefree polynomial.
/// https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
fn distinct_degree_factors<F: PrimeField>(
    f: &DensePolynomial<F>,
) -> BTreeMap<Degree, DensePolynomial<F>> {
    let x = DensePolynomial::from_coefficients_slice(&[F::zero(), F::one()]);
    let mut res = BTreeMap::new();
    let mut f_star = f.clone();
    let mut i = 1;
    while f_star.degree() >= 2 * i {
        // TODO: only works for prime fields. Won't work for extension fields
        let p: BigUint = F::BasePrimeField::MODULUS.into();
        let xp = pow_mod(&x, p, &f_star);
        let xpi = pow_mod(&xp, i.into(), &f_star);
        let g = gcd(&f_star, &(&xpi - &x));
        if g.degree() != 0 {
            f_star = &f_star / &g;
            let prev = res.insert(i, g);
            assert!(prev.is_none());
        }
        i += 1;
    }
    if f_star.degree() != 0 {
        res.insert(f_star.degree(), f_star);
    } else if res.is_empty() {
        res.insert(1, f_star);
    }
    res
}

/// Takes as input a polynomial which is known the be the product of irreducible
/// polynomials of degree d. https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
fn equal_degree_factorization<F: PrimeField>(
    f: &DensePolynomial<F>,
    d: Degree,
) -> Vec<DensePolynomial<F>> {
    if d == 0 {
        return vec![];
    }
    let f = f.clone();
    let n = f.degree();
    let r = n / d;
    let one = DensePolynomial::from_coefficients_slice(&[F::one()]);
    let mut factors = vec![f.clone()];
    let mut rng = rand::thread_rng();
    while factors.len() < r {
        let h = rand_poly(n - 1, &mut rng);
        // TODO: only works for prime fields. Won't work for extension fields
        let p: BigUint = F::BasePrimeField::MODULUS.into().pow(d as u32);
        let g = &pow_mod(&h, (p - 1u32) / 2u32, &f) - &one;
        factors = factors
            .into_iter()
            .flat_map(|factor| {
                let gcd_res = gcd(&g, &factor);
                if gcd_res.degree() != 0 && gcd_res != factor {
                    vec![&factor / &gcd_res, gcd_res]
                } else {
                    vec![factor]
                }
            })
            .collect();
    }
    factors
}

/// Returns the polynomial's square free factorization.
/// https://mathrefresher.blogspot.com/2009/01/greatest-common-divisor-of-polynomial.html
/// https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
fn square_free_factors<F: PrimeField>(f: &DensePolynomial<F>) -> DensePolynomial<F> {
    let f_prime = derivative(f);
    if f_prime.is_zero() {
        // f is already square free
        f.clone()
    } else {
        // divide out gcd(f, f') to get square free factors
        f / &gcd(f, &f_prime)
    }
}

/// Returns the GCD of two polynomials.
/// The GCD is normalized to a monic polynomial.
/// https://math.stackexchange.com/a/3009888/393151
pub fn gcd<F: PrimeField>(a: &DensePolynomial<F>, b: &DensePolynomial<F>) -> DensePolynomial<F> {
    if a.is_zero() {
        DensePolynomial::zero()
    } else if b.is_zero() {
        let leading_coeff = a.last().unwrap();
        a * leading_coeff.inverse().unwrap()
    } else {
        gcd(b, &div_remainder(a, b))
    }
}

/// Returns numerator % denominator
fn div_remainder<F: PrimeField>(
    numerator: &DensePolynomial<F>,
    denominator: &DensePolynomial<F>,
) -> DensePolynomial<F> {
    let numerator = DenseOrSparsePolynomial::from(numerator);
    let denominator = DenseOrSparsePolynomial::from(denominator);
    numerator.divide_with_q_and_r(&denominator).unwrap().1
}

/// Calculates (a^exp) % modulus
fn pow_mod<F: PrimeField>(
    a: &DensePolynomial<F>,
    mut exp: BigUint,
    modulus: &DensePolynomial<F>,
) -> DensePolynomial<F> {
    let one = DensePolynomial::from_coefficients_slice(&[F::one()]);
    let mut res = one;
    let mut acc = a.clone();
    while !exp.is_zero() {
        if exp.is_odd() {
            res = div_remainder(&res.naive_mul(&acc), modulus);
        }
        acc = div_remainder(&acc.naive_mul(&acc), modulus);
        exp >>= 1;
    }
    res
}

/// Calculates and returns the derivative f' of f
fn derivative<F: PrimeField>(f: &DensePolynomial<F>) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(
        f.iter()
            .enumerate()
            .skip(1)
            .map(|(pow, coeff)| F::from(pow as u32) * coeff)
            .collect(),
    )
}

fn rand_poly<F: PrimeField, R: Rng>(d: Degree, rng: &mut R) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec((0..=d).map(|_| F::rand(rng)).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::One;
    use ark_ff_optimized::fp31::Fp;

    #[test]
    fn finds_roots_of_cubic() {
        // = x^3 + 16*x
        let f = DensePolynomial::from_coefficients_slice(&[
            Fp::zero(),
            -Fp::from(4u8),
            Fp::zero(),
            Fp::one(),
        ]);

        let actual = find_roots(&f);

        let expected = vec![Fp::zero(), Fp::from(2u32), Fp::from(2147483645u32)];
        assert_eq!(expected, actual);
    }
}
