use crate::fftree::RationalMap;
use crate::utils::find_roots;
use crate::utils::is_odd;
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
use ark_serialize::Valid;
use core::fmt::Debug;
use core::ops::Add;
use core::ops::AddAssign;
use core::ops::Mul;
use core::ops::Neg;
use num_bigint::BigUint;
use num_integer::Integer;

/// Good curve from ECFFT Part II
/// <https://www.math.toronto.edu/swastik/ECFFT2.pdf>
/// All Good Curves share the same 2-torsion point (0, 0)
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum GoodCurve<F: Field> {
    /// Even sized fields:
    /// E_B: y^2 + x*y = x^3 + B*x, B = b^2
    Even { b: F },
    /// Odd sized fields:
    /// E_{a, B}: y^2 = x^3 + a*x^2 + B*x, B = b^2
    Odd { a: F, b: F },
}

impl<F: Field> GoodCurve<F> {
    pub fn new_odd(a: F, bb: F) -> Self {
        assert!(is_odd::<F>());
        // requirement for the curve to be non-singluar
        assert!(!bb.is_zero() && !(a.square() - bb.double().double()).is_zero());
        let b = bb.sqrt().expect("must be a quadratic residue");
        assert!((a + b + b).sqrt().is_some());
        Self::Odd { a, b }
    }

    pub fn new_even(bb: F) -> Self {
        assert!(!is_odd::<F>());
        assert!(!bb.is_zero());
        let b = bb.sqrt().expect("must be a quadratic residue");
        Self::Even { b }
    }

    pub fn good_point(self) -> Point<Self> {
        match self {
            Self::Even { b } => Point::new(b, b, self),
            Self::Odd { a, b } => Point::new(a, b.square(), self),
        }
    }

    pub fn good_isogeny(self) -> Isogeny<Self> {
        match self {
            GoodCurve::Even { b } => {
                let domain = self;
                let codomain = Self::new_even(b);

                let bb = b.square();
                let one = F::one();
                let zero = F::zero();
                let r = RationalMap::new(&[bb, zero, one], &[zero, one]);
                let g = RationalMap::new(&[bb, b], &[zero, one]);
                let h = RationalMap::new(&[bb, zero, one], &[zero, zero, one]);
                Isogeny::new(domain, codomain, r, g, h)
            }
            GoodCurve::Odd { a, b } => {
                let domain = self;
                let bb = b.square();
                let a_prime = a + b.double().double() + b.double();
                let b_prime = (a * b).double().double() + bb.double().double().double();
                let codomain = Self::new_odd(a_prime, b_prime);

                let one = F::one();
                let zero = F::zero();
                let r = RationalMap::new(&[bb, -b.double(), one], &[zero, one]);
                let g = RationalMap::zero();
                let h = RationalMap::new(&[-bb, zero, one], &[zero, zero, one]);
                Isogeny::new(domain, codomain, r, g, h)
            }
        }
    }
}

impl<F: Field> CanonicalSerialize for GoodCurve<F> {
    fn serialize_with_mode<W: ark_serialize::Write>(
        &self,
        mut writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        match self {
            GoodCurve::Even { b } => b.serialize_with_mode(writer, compress),
            GoodCurve::Odd { a, b } => {
                a.serialize_with_mode(&mut writer, compress)?;
                b.serialize_with_mode(writer, compress)
            }
        }
    }

    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        match self {
            GoodCurve::Even { b } => b.serialized_size(compress),
            GoodCurve::Odd { a, b } => a.serialized_size(compress) + b.serialized_size(compress),
        }
    }
}

impl<F: Field> Valid for GoodCurve<F> {
    #[inline]
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        Ok(())
    }
}

impl<F: Field> CanonicalDeserialize for GoodCurve<F> {
    fn deserialize_with_mode<R: ark_serialize::Read>(
        mut reader: R,
        compress: ark_serialize::Compress,
        validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        Ok(if is_odd::<F>() {
            Self::Odd {
                a: F::deserialize_with_mode(&mut reader, compress, validate)?,
                b: F::deserialize_with_mode(reader, compress, validate)?,
            }
        } else {
            Self::Even {
                b: F::deserialize_with_mode(reader, compress, validate)?,
            }
        })
    }
}

impl<F: Field> WeierstrassCurve for GoodCurve<F> {
    type F = F;

    fn a1(&self) -> Self::F {
        match self {
            GoodCurve::Even { b: _ } => F::one(),
            GoodCurve::Odd { a: _, b: _ } => F::zero(),
        }
    }

    fn a2(&self) -> Self::F {
        match self {
            GoodCurve::Even { b: _ } => F::zero(),
            GoodCurve::Odd { a, b: _ } => *a,
        }
    }

    fn a3(&self) -> Self::F {
        F::zero()
    }

    fn a4(&self) -> Self::F {
        match self {
            GoodCurve::Even { b } => b.square(),
            GoodCurve::Odd { a: _, b } => b.square(),
        }
    }

    fn a6(&self) -> Self::F {
        F::zero()
    }
}

// Finds a chain of isogenies where all curves and all isogenies are good
// TODO: could implement this on Point<F>!?
pub fn find_isogeny_chain<F: Field>(generator: Point<GoodCurve<F>>) -> Vec<Isogeny<GoodCurve<F>>> {
    let k = two_adicity(generator).expect("not a point of order 2^k");
    let mut good_isogenies = Vec::new();
    let mut g = generator;
    for _ in 0..k {
        let good_isogeny = g.curve.unwrap().good_isogeny();
        let g_prime = good_isogeny.map(&g);
        assert_eq!(two_adicity(g).unwrap(), two_adicity(g_prime).unwrap() + 1);
        good_isogenies.push(good_isogeny);
        g = g_prime;
    }
    good_isogenies
}

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

    pub fn two_isogenies(&self) -> Vec<Isogeny<Self>>
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

                let r = RationalMap::new(x_map_numerator, x_map_denominator);
                let g = RationalMap::zero();
                let h = RationalMap::new(y_map_numerator, y_map_denominator);
                Isogeny::new(domain, codomain, r, g, h)
            })
            .collect()
    }

    /// Returns all non-zero points on the curve that have order 2.
    fn two_torsion_points(&self) -> Vec<Point<Self>>
    where
        F: PrimeField,
    {
        // The two torsion points have a vertical tangent since 2*P = 0.
        // The points with vertical tangent are those with y = 0.
        // We can find the points if we find the values of x satisfy 0 = x^3 + a*x + b.
        // TODO: can use Cantor–Zassenhaus algorithm for degree 3 polynomials
        // https://en.wikipedia.org/wiki/Cantor%E2%80%93Zassenhaus_algorithm
        let roots = find_roots(&self.x3_ax_b());
        roots
            .into_iter()
            .map(|root| Point::new(root, F::zero(), *self))
            .collect()
    }

    // Returns the polynomial x^3 + A*x + B
    pub fn x3_ax_b(&self) -> DensePolynomial<F> {
        DensePolynomial::from_coefficients_vec(vec![self.b, self.a, F::zero(), F::one()])
    }
}

impl<F: Field> WeierstrassCurve for ShortWeierstrassCurve<F> {
    type F = F;

    fn a1(&self) -> F {
        F::zero()
    }

    fn a2(&self) -> F {
        F::zero()
    }

    fn a3(&self) -> F {
        F::zero()
    }

    fn a4(&self) -> F {
        self.a
    }

    fn a6(&self) -> F {
        self.b
    }
}

/// General Weierstrass curve of the form:
/// `y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6`
pub trait WeierstrassCurve:
    Clone + Copy + Debug + CanonicalDeserialize + CanonicalSerialize + Eq + PartialEq + PartialOrd + Ord
{
    type F: Field;

    /// Returns Weierstrass equation coefficient `a1`
    fn a1(&self) -> Self::F;

    /// Returns Weierstrass equation coefficient `a2`
    fn a2(&self) -> Self::F;

    /// Returns Weierstrass equation coefficient `a3`
    fn a3(&self) -> Self::F;

    /// Returns Weierstrass equation coefficient `a4`
    fn a4(&self) -> Self::F;

    /// Returns Weierstrass equation coefficient `a6`
    fn a6(&self) -> Self::F;
}

/// Defines an isogeny between curves:
/// ϕ(x, y) = (r(x), g(x) + h(x) * y)
// TODO: r, g and h are pretty confusing symbols
#[derive(Clone, Debug, CanonicalDeserialize, CanonicalSerialize)]
pub struct Isogeny<C: WeierstrassCurve> {
    pub domain: C,
    pub codomain: C,
    pub r: RationalMap<C::F>,
    pub g: RationalMap<C::F>,
    pub h: RationalMap<C::F>,
}

impl<C: WeierstrassCurve> Isogeny<C> {
    /// Isogeny ϕ(x, y) = (r(x), g(x) + h(x) * y)
    pub fn new(
        domain: C,
        codomain: C,
        r: RationalMap<C::F>,
        g: RationalMap<C::F>,
        h: RationalMap<C::F>,
    ) -> Self {
        Self {
            domain,
            codomain,
            r,
            g,
            h,
        }
    }

    pub fn map(&self, p: &Point<C>) -> Point<C> {
        if p.is_zero() {
            // TODO: should this be handled differently?
            return Point::zero();
        }

        assert_eq!(self.domain, p.curve.unwrap());
        let rx = self.r.map(&p.x);
        let gx = self.g.map(&p.x);
        let hx = self.h.map(&p.x);
        match (rx, gx, hx) {
            (Some(rx), Some(gx), Some(hx)) => Point::new(rx, gx + hx * p.y, self.codomain),
            _ => Point::zero(),
        }
    }
}

/// Point on an elliptic curve
#[derive(Clone, Copy, Debug)]
pub struct Point<C: WeierstrassCurve> {
    pub x: C::F,
    pub y: C::F,
    pub curve: Option<C>,
}

impl<C: WeierstrassCurve> Point<C> {
    pub const fn new(x: C::F, y: C::F, curve: C) -> Self {
        let curve = Some(curve);
        Self { x, y, curve }
    }
}

impl<C: WeierstrassCurve> Add for Point<C> {
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
            let a1 = curve.a1();
            let a2 = curve.a2();
            let a3 = curve.a3();
            let a4 = curve.a4();
            let a6 = curve.a6();
            let x1 = self.x;
            let y1 = self.y;
            let x2 = rhs.x;
            let y2 = rhs.y;

            if x1 == x2 && (y1 + y2 + a1 * x2 + a3).is_zero() {
                Self::zero()
            } else {
                let lambda: C::F;
                let nu: C::F;
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

impl<C: WeierstrassCurve> AddAssign for Point<C> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl<C: WeierstrassCurve> Mul<BigUint> for Point<C> {
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

impl<C: WeierstrassCurve> Neg for Point<C> {
    type Output = Self;

    fn neg(self) -> Self {
        // Source: "The arithmetic of elliptic curves, 2nd ed." Silverman, III.2.3
        if self.is_zero() {
            return self;
        }

        let curve = self.curve.unwrap();
        Self {
            y: -self.y - curve.a1() * self.x - curve.a3(),
            ..self
        }
    }
}

impl<C: WeierstrassCurve> PartialEq for Point<C> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_zero() && other.is_zero() {
            true
        } else {
            assert_eq!(self.curve, other.curve);
            self.x == other.x && self.y == other.y
        }
    }
}

impl<C: WeierstrassCurve> Zero for Point<C> {
    fn zero() -> Self {
        Self {
            x: C::F::zero(),
            y: C::F::zero(),
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
    subgroup_generator: Point<ShortWeierstrassCurve<F>>,
    subgroup_order: usize,
    coset_offset: Point<ShortWeierstrassCurve<F>>,
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
        let isogeny = g
            .curve
            .unwrap()
            .two_isogenies()
            .into_iter()
            .find_map(|isogeny| {
                let g_prime = isogeny.map(&g);
                if two_adicity(g)? == two_adicity(g_prime)? + 1 {
                    g = g_prime;
                    Some(isogeny)
                } else {
                    None
                }
            })
            .expect("cannot find a suitable isogeny");
        rational_maps.push(isogeny.r)
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
                assert!(isogeny.r.map(&p.x).is_none());
            }
        }
    }
}
