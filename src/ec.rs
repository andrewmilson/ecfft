use crate::utils::find_roots;
use ark_ff::PrimeField;
use ark_ff::Zero;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use num_bigint::BigUint;
use num_integer::Integer;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Mul;
use std::ops::Neg;

/// Curve of the form y^2 = x^3 + a*x + b
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Curve<F> {
    a: F,
    b: F,
}

impl<F: PrimeField> Curve<F> {
    pub fn new(a: F, b: F) -> Self {
        Self { a, b }
    }

    pub fn two_isogenies(&self) -> Vec<Isogeny<F>> {
        // Find all the points of order two and then use Velu's formula to find all
        // the 2-isogenies. Velu's formula: https://math.mit.edu/~drew/CTNT2018.pdf
        self.two_torsion_points()
            .into_iter()
            .map(|point| {
                let x0 = point.x;
                let t = F::from(3u8) * x0 * x0 + self.a;

                // Our domain E(F) and codomain E'(F)
                let domain = *self;
                let codomain = Self::new(self.a - F::from(5u8) * t, self.b - F::from(7u8) * x0 * t);

                // ϕ: E -> E' := ((x^2 - x0*x + t)/(x - x0), ((x - x0)^2 - t)/(x - x0)^2 * y)
                let x_map_numerator = &[t, -x0, F::one()];
                let x_map_denominator = &[-x0, F::one()];
                let y_map_numerator = &[x0 * x0 - t, -(x0 + x0), F::one()];
                let y_map_denominator = &[x0 * x0, -(x0 + x0), F::one()];

                Isogeny::new(
                    domain,
                    codomain,
                    DensePolynomial::from_coefficients_slice(x_map_numerator),
                    DensePolynomial::from_coefficients_slice(x_map_denominator),
                    DensePolynomial::from_coefficients_slice(y_map_numerator),
                    DensePolynomial::from_coefficients_slice(y_map_denominator),
                )
            })
            .collect()
    }

    /// Returns all non-zero points on the curve that have order 2.
    fn two_torsion_points(&self) -> Vec<Point<F>> {
        // The two torsion points have a vertical tangent since 2*P = 0.
        // The points with vertical tangent are those with y = 0.
        // We can find the points if we find the values of x satisfy 0 = x^3 + a*x + b.
        let x3_ax_b = &[self.b, self.a, F::zero(), F::one()];
        let roots = find_roots(&DensePolynomial::from_coefficients_slice(x3_ax_b));
        roots
            .into_iter()
            .map(|root| Point::new(root, F::zero(), *self))
            .collect()
    }
}

// Defines an isogeny between curves
pub struct Isogeny<F: PrimeField> {
    pub domain: Curve<F>,
    pub codomain: Curve<F>,
    pub x_numerator_map: DensePolynomial<F>,
    pub x_denominator_map: DensePolynomial<F>,
    pub y_numerator_map: DensePolynomial<F>,
    pub y_denominator_map: DensePolynomial<F>,
}

impl<F: PrimeField> Isogeny<F> {
    pub fn new(
        domain: Curve<F>,
        codomain: Curve<F>,
        x_numerator_map: DensePolynomial<F>,
        x_denominator_map: DensePolynomial<F>,
        y_numerator_map: DensePolynomial<F>,
        y_denominator_map: DensePolynomial<F>,
    ) -> Self {
        Self {
            domain,
            codomain,
            x_numerator_map,
            x_denominator_map,
            y_numerator_map,
            y_denominator_map,
        }
    }

    pub fn map_x(&self, x: &F) -> Option<F> {
        Some(self.x_numerator_map.evaluate(x) * self.x_denominator_map.evaluate(x).inverse()?)
    }
}

/// Point on an elliptic curve
#[derive(Clone, Copy, Debug)]
pub struct Point<F> {
    pub x: F,
    pub y: F,
    pub curve: Option<Curve<F>>,
}

impl<F: PrimeField> Point<F> {
    pub fn new(x: F, y: F, curve: Curve<F>) -> Self {
        let curve = Some(curve);
        Self { x, y, curve }
    }
}

impl<F: PrimeField> Add for Point<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        // Addition law for elliptic curve groups
        // Source: Elliptic Curves book by LCW
        if self.is_zero() {
            rhs
        } else if rhs.is_zero() {
            self
        } else if self.curve != rhs.curve {
            panic!("points belong to different curves");
        } else if self.x == rhs.x && (self.y != rhs.y || rhs.y.is_zero()) {
            Self::zero()
        } else {
            let curve = self.curve.unwrap();
            let x1 = self.x;
            let y1 = self.y;
            let x2 = rhs.x;
            let y2 = rhs.y;
            if x1 == x2 {
                // use tangent line
                let x1x1 = x1 * x1;
                let m = ((x1x1 + x1x1 + x1x1) + curve.a) / y1.double();
                let x3 = m * m - (x1 + x1);
                let y3 = m * (x1 - x3) - y1;
                Self::new(x3, y3, curve)
            } else {
                // calculate slope through a and b
                let dx = x2 - x1;
                let dy = y2 - y1;
                let m = dy / dx;
                let x3 = m * m - x2 - x1;
                let y3 = m * (x1 - x3) - y1;
                Self::new(x3, y3, curve)
            }
        }
    }
}

impl<F: PrimeField> AddAssign for Point<F> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl<F: PrimeField> Mul<BigUint> for Point<F> {
    type Output = Self;

    fn mul(self, mut rhs: BigUint) -> Self {
        let mut res = Self::zero();
        let mut acc = self;
        while !rhs.is_zero() {
            if rhs.is_odd() {
                res = res + acc;
            }
            acc = acc + acc;
            rhs >>= 1;
        }
        res
    }
}

impl<F: PrimeField> Neg for Point<F> {
    type Output = Self;

    fn neg(self) -> Self {
        Self { y: -self.y, ..self }
    }
}

impl<F: PrimeField> PartialEq for Point<F> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_zero() && other.is_zero() {
            true
        } else {
            assert_eq!(self.curve, other.curve);
            self.x == other.x && self.y == other.y
        }
    }
}

impl<F: PrimeField> Zero for Point<F> {
    fn zero() -> Self {
        Self {
            x: F::zero(),
            y: F::zero(),
            curve: None,
        }
    }

    fn is_zero(&self) -> bool {
        self.curve.is_none()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::One;
    use ark_ff_optimized::fp31::Fp;

    #[test]
    fn two_torsion_points_are_degree_two() {
        let curve = Curve::new(Fp::one(), Fp::zero());

        let two_torsion_points = curve.two_torsion_points();

        for p in two_torsion_points {
            assert!(!p.is_zero());
            assert!((p + p).is_zero());
        }
    }

    #[test]
    fn two_isogenies_map_to_identity() {
        let curve = Curve::new(Fp::one(), Fp::zero());
        let two_torsion_points = curve.two_torsion_points();

        let two_isogenies = curve.two_isogenies();

        for p in two_torsion_points {
            for isogeny in &two_isogenies {
                assert!(isogeny.map_x(&p.x).is_none());
            }
        }
    }
}
