#![feature(array_windows, array_chunks)]

pub mod ec;
pub mod ecfft;
mod utils;

#[cfg(test)]
mod tests {
    use crate::ec::Curve;
    use crate::ec::Point;
    use crate::ecfft::FFTree;
    use ark_ff::One;
    use ark_ff::Zero;
    use ark_ff_optimized::fp31::Fp;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use ark_poly::Polynomial;
    use num_bigint::BigUint;

    #[test]
    fn it_works() {
        let curve = Curve::new(Fp::one(), Fp::zero());
        let n = 1u32 << 3;
        let cardinality = BigUint::from(1u32 << 31);
        let primitive_point = Point::new(Fp::from(1048755163u32), Fp::from(279503108u32), curve);
        let generator = primitive_point * (cardinality / n);
        let fftree = FFTree::new(primitive_point, generator);

        let one = Fp::one();
        let two = one + one;
        let poly = DensePolynomial::from_coefficients_vec(vec![one, two, one, two]);
        let mut ys = Vec::new();
        for x in fftree.f.leaves() {
            ys.push(poly.evaluate(x))
        }

        println!("YS: {:?}", ys);
        let original_ys = ys.iter().copied().step_by(2).collect::<Vec<Fp>>();
        println!("original: {:?}", original_ys);
        let expected_ys = ys.iter().copied().skip(1).step_by(2).collect::<Vec<Fp>>();
        let actual_ys = fftree.extend(&original_ys);
        println!("YOOYOYO: {:?}", actual_ys);

        println!("ENTER: {:?}", fftree.enter(&poly));

        // // let fftree = FFForrest::new(primitive_point, generator);

        // println!("FFTree: {:?}", fftree);

        // let one = Fp::one();
        // let two = one + one;
        // let poly = DensePolynomial::from_coefficients_vec(vec![one, two, one,
        // two]); // two, two
        // let layer = ffforrest.get_layer(0);
        // // let mut ys = Vec::new();
        // // for [x, _] in layer.l.array_chunks() {
        // //     ys.push(poly.evaluate(x));
        // // }

        // for [x, _] in
        // crate::utils::get_layers(&ffforrest.vertices)[1].array_chunks() {
        //     println!("OUTPUT: {}", poly.evaluate(x));
        // }
        // for x in crate::utils::get_layers(&ffforrest.vertices)[3] {
        //     println!("BOT: {}", poly.evaluate(x));
        // }
        // for x in crate::utils::get_layers(&ffforrest.vertices)[2] {
        //     println!("SHOT: {}", poly.evaluate(x));
        // }

        // let curve = Curve::new(Fp::one(), Fp::zero());
        // let n = 1u32 << 10;
        // let cardinality = BigUint::from(1u32 << 31);
        // let primitive_point = Point::new(Fp::from(1048755163u32),
        // Fp::from(279503108u32), curve); let generator =
        // primitive_point * (cardinality / n); let ffforrest =
        // FFForrest::new(primitive_point, generator);
        // let one = Fp::one();
        // let two = one + one;
        // let poly = DensePolynomial::from_coefficients_vec(vec![one, two, one,
        // two]); // two, two let layer = ffforrest.get_layer(0);
        // // let mut ys = Vec::new();
        // // for [x, _] in layer.l.array_chunks() {
        // //     ys.push(poly.evaluate(x));
        // // }

        // for [x, _] in
        // crate::utils::get_layers(&ffforrest.vertices)[1].array_chunks() {
        //     println!("OUTPUT: {}", poly.evaluate(x));
        // }
        // for x in crate::utils::get_layers(&ffforrest.vertices)[3] {
        //     println!("BOT: {}", poly.evaluate(x));
        // }
        // for x in crate::utils::get_layers(&ffforrest.vertices)[2] {
        //     println!("SHOT: {}", poly.evaluate(x));
        // }

        // let actual = ffforrest.extend(&ys);

        // // let mut expected = Vec::new();
        // // for [_, x] in layer.l.array_chunks() {
        // //     expected.push(poly.evaluate(x));
        // // }
        // // println!("ys: {:?}", ys);
        // // println!("actual: {:?}", actual);
        // // println!("expected: {:?}", expected);
        // // assert_eq!(expected, actual);

        // // let enter_res = ffforrest.enter(&[
        // //     one, two, one,
        // //     two, /* two,
        // //          * two,
        // //          * Fp::zero(),
        // //          * Fp::zero(),
        // //          * Fp::zero(),
        // //          * Fp::zero() */
        // // ]);
        // // println!("entered: {:?}", enter_res);
        // // println!("EXT: {:?}", ffforrest.extend(&enter_res));
        // // let enter_res2 = ffforrest.enter(&[
        // //     one,
        // //     two,
        // //     // two,
        // //     // two,
        // //     // Fp::zero(),
        // //     // Fp::zero(),
        // //     Fp::zero(),
        // //     Fp::zero(),
        // // ]);
        // // println!("entered: {:?}", enter_res2);
    }
}
