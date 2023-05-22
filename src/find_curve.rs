use crate::ec::GoodCurve;
use crate::ec::Point;
use crate::utils::is_odd;
use ark_ff::Field;
use rand::Rng;

/// Returns the x-coordinate of Q if P is a half point of Q. Otherwise returns
/// None. y^2 = x(x^2 + a*x + B), B = b^2. `F` must be an odd order field.
/// <https://www.ams.org/journals/mcom/2005-74-249/S0025-5718-04-01640-0/S0025-5718-04-01640-0.pdf>
/// "If Q = 2P, we say that Q has a half point P".
fn double_point_x<F: Field>(px: F, a: F, bb: F) -> Option<F> {
    assert!(is_odd::<F>());
    let pxpx = px.square();
    let pypy = px * (pxpx + a * px + bb);
    if pypy.is_zero() {
        return None;
    }
    Some((pxpx - bb).square() / pypy.double().double())
}

/// Returns the x-coordinate of P if P is a half point of Q. Otherwise returns
/// None. y^2 = x(x^2 + a*x + B), B = b^2. `F` must be an odd order field.
/// <https://www.ams.org/journals/mcom/2005-74-249/S0025-5718-04-01640-0/S0025-5718-04-01640-0.pdf>
/// "If Q = 2P, we say that Q has a half point P".
fn half_point_x<F: Field>(qx: F, a: F, bb: F) -> Option<F> {
    assert!(is_odd::<F>());
    let roots = fi_roots(1, qx, a, bb).or_else(|| fi_roots(2, qx, a, bb))?;
    roots
        .into_iter()
        .find(|&x| (x * (x.square() + a * x + bb)).sqrt().is_some())
}

/// Finds the roots of a quadratic equation a*x^2 + b*x + c
/// F must be an odd order field
fn roots<F: Field>(a: F, b: F, c: F) -> Option<[F; 2]> {
    assert!(is_odd::<F>());
    let discriminant = b.square() - (a * c).double().double();
    let discriminant_sqrt = discriminant.sqrt()?;
    let two_inv = a.double().inverse()?;
    Some([
        (-b + discriminant_sqrt) * two_inv,
        (-b - discriminant_sqrt) * two_inv,
    ])
}

/// Finds the roots of f_{1, ξ} or f_{2, ξ} if they exists
/// Let Q be a point such that x(Q) = ξ
/// F must be an odd order field
fn fi_roots<F: Field>(i: usize, qx: F, a: F, bb: F) -> Option<[F; 2]> {
    assert!(i == 1 || i == 2);
    let delta = qx.square() + a * qx + bb;
    let delta_sqrt = delta.sqrt()?;
    let one = F::one();
    let x_coeff = -(qx.double() + (-one).pow([i as u64]) * delta_sqrt.double());
    roots(one, x_coeff, bb)
}

/// Determines the 2-sylow subgroup of curve y^2 = x(x^2 + a*x + B), B = b^2
/// Algorithm from: https://www.ams.org/journals/mcom/2005-74-249/S0025-5718-04-01640-0/S0025-5718-04-01640-0.pdf
/// Output is of the form (n, r, x1, x2) with n, r such that Sylow_2(E(Fq)) is
/// isomorphic to Z/2^nZ × Z/2^rZ and x-coordinates x1, x2 of points of order
/// 2^n and 2^r, respectively, generating this Sylow subgroup. Note that None
/// refers to infinity. Panics if `F` is not a field with odd order.
// TODO: finish this implementation
fn _find_two_sylow_subgroup<F: Field>(a: F, bb: F) -> (u32, u32, Option<F>, Option<F>) {
    assert!(is_odd::<F>());
    let discriminant = a.square() - bb.double().double();
    assert!(!discriminant.is_zero());
    let one = F::one();
    let zero = F::zero();

    if discriminant.sqrt().is_some() {
        // == non-cyclic case ==

        // the roots (r0, r1, r2) of our curve equation give us our 2-torsion points
        // x^3 + a*x^2 + B = x(x^2 + a*x + B), B = b^2
        let Some([r1, r2]) = roots(F::one(), a, bb) else { unreachable!() };
        let r0 = zero;

        // check if rational points of order 4 exist (lemma 3)
        let p4_cond1 = bb.sqrt().and_then(|b| (a - b.double()).sqrt()).is_some();
        let p4_cond2 = r1.sqrt().and_then(|_| (r1.double() + a).sqrt()).is_some();
        let p4_cond3 = r2.sqrt().and_then(|_| (r2.double() + a).sqrt()).is_some();

        if !p4_cond1 && !p4_cond2 && !p4_cond3 {
            // there is no point with order 4
            return (1, 1, Some(zero), Some(r1));
        }

        // find an x-coordinate of points of order 4 (lemma 3)
        // safe to unwrap since we've already checked if points of order 4 exist
        let p4_x0 = roots(one, r0, -bb).map(|[x, _]| x).unwrap();
        let p4_x1 = roots(one, -r1.double(), bb).map(|[x, _]| x).unwrap();
        let p4_x2 = roots(one, -r2.double(), bb).map(|[x, _]| x).unwrap();

        let has_half_point = |x: F| -> bool {
            // halving conditions (Proposition 2)
            if x.sqrt().is_none() {
                return false;
            }
            let delta = x.square() + a * x + bb;
            let delta_sqrt = match delta.sqrt() {
                Some(sqrt) => sqrt,
                None => return false,
            };
            (x.double() + a + delta_sqrt.double()).sqrt().is_some()
        };

        let mut k = 2;
        let mut h = 1;
        let mut acc0 = p4_x0;
        let mut acc1 = p4_x1;
        let mut acc2 = p4_x2;
        loop {
            let acc0_has_half_point = has_half_point(acc0);
            let acc1_has_half_point = has_half_point(acc1);
            let acc2_has_half_point = has_half_point(acc2);

            if !acc0_has_half_point && !acc1_has_half_point && !acc2_has_half_point {
                // no accumulators have a half point
                return (h, h, Some(acc0), Some(acc1));
            }

            if acc0_has_half_point && acc1_has_half_point && acc2_has_half_point {
                // all accumulators have a half point
                acc0 = half_point_x(acc0, a, bb).unwrap();
                acc1 = half_point_x(acc1, a, bb).unwrap();
                acc2 = half_point_x(acc2, a, bb).unwrap();
                k += 1;
                h = k;
                continue;
            }

            if h == 1 {
                // `x2` is the x-coord of a point of order 2 that does not have a half point
                // `acc` is TODO
                let _x2 = match (
                    acc0_has_half_point,
                    acc1_has_half_point,
                    acc2_has_half_point,
                ) {
                    (false, _, _) => r0,
                    (_, false, _) => r1,
                    (_, _, false) => r2,
                    _ => unreachable!(),
                };

                todo!()
            }

            todo!()
        }

        // TODO
    } else if let Some(b) = bb.sqrt() {
        // == cyclic case ==
        let p4_x = match ((a + b.double()).sqrt(), (a - b.double()).sqrt()) {
            (Some(_), _) => b,
            (_, Some(_)) => -b,
            _ => unreachable!(),
        };

        if double_point_x(p4_x, a, bb).is_none() {
            return (1, 0, Some(zero), None);
        }

        let mut k = 2;
        let mut acc = p4_x;
        while let Some(x) = half_point_x(acc, a, bb) {
            k += 1;
            acc = x;
        }

        (k, 0, Some(acc), None)
    } else {
        (0, 0, None, None)
    }
}

/// Determines the 2-sylow subgroup of curve y^2 = x(x^2 + a*x + B), B = b^2
/// Output is of the form (n, x) with n such that Sylow_2(E(Fq)) is
/// isomorphic to Z/2^nZ and x-coordinates x, the x-coordinate of the point
/// of order 2^n generating this Sylow subgroup. Note that None refers to
/// infinity. Panics if `F` is not a field with odd order.
///
/// Algorithm from: https://www.ams.org/journals/mcom/2005-74-249/S0025-5718-04-01640-0/S0025-5718-04-01640-0.pdf
/// note that the algorithm in the paper finds the 2-slyow subgoup of non-cyclic
/// curves as well. This has been ommited since we are only interested in cyclic
/// subgroups for ECFFT.
fn cyclic_two_sylow_subgroup<F: Field>(a: F, bb: F) -> (u32, Option<F>) {
    assert!(is_odd::<F>());
    let discriminant = a.square() - bb.double().double();
    assert!(!discriminant.is_zero());

    if let (Some(b), None) = (bb.sqrt(), discriminant.sqrt()) {
        // 2-sylow subgroup is cyclic
        let p4_x = match ((a + b.double()).sqrt(), (a - b.double()).sqrt()) {
            (Some(_), _) => b,
            (_, Some(_)) => -b,
            _ => unreachable!(),
        };

        if double_point_x(p4_x, a, bb).is_none() {
            return (1, Some(F::zero()));
        }

        let mut k = 2;
        let mut acc = p4_x;
        while let Some(x) = half_point_x(acc, a, bb) {
            k += 1;
            acc = x;
        }

        (k, Some(acc))
    } else {
        (0, None)
    }
}

/// Returns a point with order 2^n on a Good Curve such that n >= k.
/// Output is of the form (n, subgroup_generator_point)
/// Based on `find_curve` algorithm from "ECFFT part II":
/// <https://www.math.toronto.edu/swastik/ECFFT2.pdf>
pub fn find_curve<F: Field>(mut rng: impl Rng, k: u32) -> (u32, Point<GoodCurve<F>>) {
    let k = k.max(2);
    if is_odd::<F>() {
        loop {
            // curve: y^2 = x*(x^2 + a*x + b)
            let a = F::rand(&mut rng);
            let bb = F::rand(&mut rng);
            // TODO: may be faster this way but would be nice to
            // finish the generic two sylow algorithm.
            let (n, x) = cyclic_two_sylow_subgroup(a, bb);
            if n >= k {
                let x = x.unwrap();
                let yy = x * (x.square() + a * x + bb);
                let y = yy.sqrt().unwrap();
                let good_curve = GoodCurve::new_odd(a, bb);
                let p = Point::new(x, y, good_curve);
                return (n, p);
            }
        }
    } else {
        todo!()
    }
}
