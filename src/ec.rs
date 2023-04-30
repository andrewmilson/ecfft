use crate::errors::TryFromGeneralToShortWeierstrassCurveError;
use crate::fftree::RationalMap;
use crate::utils::find_roots;
use crate::utils::two_adicity;
use crate::FFTree;
use ark_ff::vec;
use ark_ff::vec::Vec;
use ark_ff::Field;
use ark_ff::PrimeField;
use ark_ff::Zero;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use core::ops::Add;
use core::ops::AddAssign;
use core::ops::Mul;
use core::ops::Neg;
use num_bigint::BigUint;
use num_integer::Integer;

/// Short Weierstrass curve of the form y^2 = x^3 + a*x + b
#[derive(
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Debug,
    Default,
    CanonicalDeserialize,
    CanonicalSerialize,
)]
pub struct ShortWeierstrassCurve<F: Field> {
    pub a: F,
    pub b: F,
}

impl<F: Field> ShortWeierstrassCurve<F> {
    pub const fn new(a: F, b: F) -> Self {
        Self { a, b }
    }

    pub fn two_isogenies(&self) -> Vec<Isogeny<F>>
    where
        F: PrimeField,
    {
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

                let x_map = RationalMap {
                    numerator_map: DensePolynomial::from_coefficients_slice(x_map_numerator),
                    denominator_map: DensePolynomial::from_coefficients_slice(x_map_denominator),
                };

                let y_map = RationalMap {
                    numerator_map: DensePolynomial::from_coefficients_slice(y_map_numerator),
                    denominator_map: DensePolynomial::from_coefficients_slice(y_map_denominator),
                };

                Isogeny::new(domain.into(), codomain.into(), x_map, y_map)
            })
            .collect()
    }

    /// Returns all non-zero points on the curve that have order 2.
    fn two_torsion_points(&self) -> Vec<Point<F>>
    where
        F: PrimeField,
    {
        // The two torsion points have a vertical tangent since 2*P = 0.
        // The points with vertical tangent are those with y = 0.
        // We can find the points if we find the values of x satisfy 0 = x^3 + a*x + b.
        let roots = find_roots(&self.x3_ax_b());
        roots
            .into_iter()
            .map(|root| Point::new(root, F::zero(), (*self).into()))
            .collect()
    }

    // Returns the polynomial x^3 + A*x + B
    pub fn x3_ax_b(&self) -> DensePolynomial<F> {
        DensePolynomial::from_coefficients_vec(vec![self.b, self.a, F::zero(), F::one()])
    }
}

impl<F: Field> TryFrom<GeneralWeierstrassCurve<F>> for ShortWeierstrassCurve<F> {
    type Error = TryFromGeneralToShortWeierstrassCurveError;

    fn try_from(general_curve: GeneralWeierstrassCurve<F>) -> Result<Self, Self::Error> {
        let GeneralWeierstrassCurve { a1, a2, a3, a4, a6 } = general_curve;
        if !a1.is_zero() || !a2.is_zero() || !a3.is_zero() {
            return Err(TryFromGeneralToShortWeierstrassCurveError);
        }

        Ok(ShortWeierstrassCurve::new(a4, a6))
    }
}

/// General Weierstrass curve of the form:
/// `y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6`
#[derive(
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Debug,
    Default,
    CanonicalDeserialize,
    CanonicalSerialize,
)]
pub struct GeneralWeierstrassCurve<F: Field> {
    pub a1: F,
    pub a2: F,
    pub a3: F,
    pub a4: F,
    pub a6: F,
}

impl<F: Field> GeneralWeierstrassCurve<F> {
    pub const fn new(a1: F, a2: F, a3: F, a4: F, a6: F) -> Self {
        Self { a1, a2, a3, a4, a6 }
    }
}

impl<F: Field> From<ShortWeierstrassCurve<F>> for GeneralWeierstrassCurve<F> {
    fn from(short_curve: ShortWeierstrassCurve<F>) -> Self {
        let zero = F::zero();
        Self::new(zero, zero, zero, short_curve.a, short_curve.b)
    }
}

// Defines an isogeny between curves
#[derive(Clone, Debug, CanonicalDeserialize, CanonicalSerialize)]
pub struct Isogeny<F: Field> {
    pub domain: GeneralWeierstrassCurve<F>,
    pub codomain: GeneralWeierstrassCurve<F>,
    pub x_map: RationalMap<F>,
    pub y_map: RationalMap<F>,
}

impl<F: Field> Isogeny<F> {
    pub fn new(
        domain: GeneralWeierstrassCurve<F>,
        codomain: GeneralWeierstrassCurve<F>,
        x_map: RationalMap<F>,
        y_map: RationalMap<F>,
    ) -> Self {
        Self {
            domain,
            codomain,
            x_map,
            y_map,
        }
    }

    pub fn map(&self, p: &Point<F>) -> Point<F> {
        if p.is_zero() {
            // TODO: should this be handled differently?
            Point::zero()
        } else {
            assert_eq!(self.domain, p.curve.unwrap());
            let x_prime = self.x_map.map(&p.x);
            let y_prime = self.y_map.map(&p.x).map(|v| v * p.y);
            match (x_prime, y_prime) {
                (Some(x), Some(y)) => Point::new(x, y, self.codomain),
                _ => Point::zero(),
            }
        }
    }
}

/// Point on an elliptic curve
#[derive(Clone, Copy, Debug)]
pub struct Point<F: Field> {
    pub x: F,
    pub y: F,
    pub curve: Option<GeneralWeierstrassCurve<F>>,
}

impl<F: Field> Point<F> {
    pub const fn new(x: F, y: F, curve: GeneralWeierstrassCurve<F>) -> Self {
        let curve = Some(curve);
        Self { x, y, curve }
    }
}

impl<F: Field> Add for Point<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        // Addition law for elliptic curve groups
        // Source: "The arithmetic of elliptic curves, 2nd ed." Silverman, III.2.3
        if self.is_zero() {
            rhs
        } else if rhs.is_zero() {
            self
        } else if self.curve != rhs.curve {
            panic!("points belong to different curves");
        } else {
            let curve = self.curve.unwrap();
            let GeneralWeierstrassCurve { a1, a2, a3, a4, a6 } = curve;
            let x1 = self.x;
            let y1 = self.y;
            let x2 = rhs.x;
            let y2 = rhs.y;

            if x1 == x2 && (y1 + y2 + a1 * x2 + a3).is_zero() {
                Self::zero()
            } else {
                let lambda: F;
                let nu: F;
                if x1 == x2 {
                    // tangent line
                    let x1x1 = x1.square();
                    let a2x1 = a2 * x1;
                    let a1x1 = a1 * x1;
                    lambda =
                        (x1x1 + x1x1 + x1x1 + a2x1 + a2x1 + a4 - a1 * y1) / (y1 + y1 + a1x1 + a3);
                    nu = (-(x1x1 * x1) + a4 * x1 + a6 + a6 - a3 * y1) / (y1 + y1 + a1 * x1 + a3);
                } else {
                    // slope through the points
                    lambda = (y2 - y1) / (x2 - x1);
                    nu = (y1 * x2 - y2 * x1) / (x2 - x1);
                }
                let x3 = lambda.square() + a1 * lambda - a2 - x1 - x2;
                let y3 = -(lambda + a1) * x3 - nu - a3;
                Self::new(x3, y3, curve)
            }
        }
    }
}

impl<F: Field> AddAssign for Point<F> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl<F: Field> Mul<BigUint> for Point<F> {
    type Output = Self;

    fn mul(self, mut rhs: BigUint) -> Self {
        let mut res = Self::zero();
        let mut acc = self;
        while !rhs.is_zero() {
            if rhs.is_odd() {
                res += acc;
            }
            acc += acc;
            rhs >>= 1;
        }
        res
    }
}

impl<F: Field> Neg for Point<F> {
    type Output = Self;

    fn neg(self) -> Self {
        // Source: "The arithmetic of elliptic curves, 2nd ed." Silverman, III.2.3
        if self.is_zero() {
            return self;
        }

        let GeneralWeierstrassCurve { a1, a3, .. } = self.curve.unwrap();
        Self {
            y: -self.y - a1 * self.x - a3,
            ..self
        }
    }
}

impl<F: Field> PartialEq for Point<F> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_zero() && other.is_zero() {
            true
        } else {
            assert_eq!(self.curve, other.curve);
            self.x == other.x && self.y == other.y
        }
    }
}

impl<F: Field> Zero for Point<F> {
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

/// Builds an FFTree (from a elliptic curve point) that's capable of evaluating
/// degree n-1 polynomials. Returns None if no such FFTree exists.
///
/// `subgroup_generator` must generate a cyclic subgroup of order 2^k.
/// `subgroup_order` is the order of the subgroup generated by the generator.
/// `coset_offset` must be disjoint from the subgroup defined by the generator.
/// `n` is the power-of-2 size of the FFTree.
pub fn build_ec_fftree<F: PrimeField>(
    subgroup_generator: Point<F>,
    subgroup_order: usize,
    coset_offset: Point<F>,
    n: usize,
) -> Option<FFTree<F>> {
    assert_ne!(coset_offset, subgroup_generator);
    assert_eq!(coset_offset.curve, subgroup_generator.curve);
    assert!(n.is_power_of_two());
    assert!(subgroup_order.is_power_of_two());
    let subgroup_two_addicity = subgroup_order.ilog2();
    let log_n = n.ilog2();
    assert!(log_n < 32);

    // FFTree size is too large for our generator
    if log_n > subgroup_two_addicity {
        return None;
    }

    // get a generator of a subgroup with order `n`
    let mut generator = subgroup_generator;
    for _ in 0..subgroup_two_addicity - log_n {
        generator += generator;
    }

    // find our rational maps
    let mut rational_maps = Vec::new();
    let mut g = generator;
    for _ in 0..log_n {
        let short_weierstrass_curve = ShortWeierstrassCurve::try_from(g.curve.unwrap()).unwrap();
        let isogeny = short_weierstrass_curve
            .two_isogenies()
            .into_iter()
            .find_map(|isogeny| {
                let g_prime = isogeny.map(&g);
                if two_adicity(g_prime)? == two_adicity(g)? - 1 {
                    g = g_prime;
                    Some(isogeny)
                } else {
                    None
                }
            })
            .expect("cannot find a suitable isogeny");
        rational_maps.push(isogeny.x_map)
    }

    // generate the FFTree leaf nodes
    let mut leaves = vec![F::zero(); n];
    let mut acc = Point::zero();
    for x in &mut leaves {
        *x = (coset_offset + acc).x;
        acc += generator;
    }

    Some(FFTree::new(leaves, rational_maps))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::One;
    use ark_ff_optimized::fp31::Fp;

    #[test]
    fn two_torsion_points_have_order_two() {
        let curve = ShortWeierstrassCurve::new(Fp::one(), Fp::zero());

        let two_torsion_points = curve.two_torsion_points();

        for p in two_torsion_points {
            assert!(!p.is_zero());
            assert!((p + p).is_zero());
        }
    }

    #[test]
    fn two_isogenies_map_to_identity() {
        let curve = ShortWeierstrassCurve::new(Fp::one(), Fp::zero());
        let two_torsion_points = curve.two_torsion_points();

        let two_isogenies = curve.two_isogenies();

        for p in two_torsion_points {
            for isogeny in &two_isogenies {
                assert!(isogeny.x_map.map(&p.x).is_none());
            }
        }
    }
}
