use ark_ff::Field;
use ark_ff::One;
use ark_ff::Zero;
use ark_poly::univariate::DenseOrSparsePolynomial;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use ecfft::m31::Fp;
use ecfft::utils::lagrange_interpolate;
use ecfft::EcFftField;
use std::env;

fn main() {
    let fftree = Fp::build_fftree(1 << 5).unwrap();
    let fftree_8 = fftree.subtree_with_size(8);
    let s_8 = fftree_8.f.get_layer(0).to_vec();
    let (s0_8, s1_8): (Vec<Fp>, Vec<Fp>) = s_8.chunks_exact(2).map(|s| (s[0], s[1])).unzip();
    let z0_8_evals = fftree_8
        .z0_s1
        .iter()
        .flat_map(|v| [Fp::zero(), *v])
        .collect::<Vec<_>>();
    let z0_8_coeffs = lagrange_interpolate(&s_8, &z0_8_evals);
    println!("Eval is: {}", z0_8_coeffs.evaluate(&s1_8[0]));
    println!("GAS: {:?}", z0_8_coeffs);
    let z0z0 = z0_8_coeffs.naive_mul(&z0_8_coeffs);
    let a_4 = DensePolynomial::from_coefficients_vec(vec![
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::one(),
    ]);
    println!("z0z0: {:?}", z0z0);
    println!("z0z0 rem a: {:?}", z0z0);
    println!("GRRR: {:?}", a_4);
    println!(
        "ACTUAL: {:?}",
        DenseOrSparsePolynomial::divide_with_q_and_r(&z0z0.into(), &a_4.into())
            .unwrap()
            .1,
    );
    println!(
        "EXPECTED: {:?}",
        lagrange_interpolate(&s_8, &fftree_8.z0z0_rem_xnn_s)
    );

    let z0_sub_xnn = z0_8_evals
        .iter()
        .zip(&fftree_8.xnn_s)
        .map(|(&z0, xnn)| (z0 - xnn).pow([4]))
        .collect::<Vec<Fp>>();
    let reduced_vals = fftree.redc_z0(&z0_sub_xnn, &fftree_8.xnn_s);
    println!("redu: {:?}", reduced_vals);
    println!("redu: {:?}", lagrange_interpolate(&s_8, &reduced_vals));

    // todo!();

    let fftree_16 = fftree.subtree_with_size(16);
    let s_16 = fftree_8.f.get_layer(0).to_vec();
    let z0_8_ext = fftree_16.extend(
        &z0_8_evals
            .iter()
            .zip(&fftree_8.xnn_s)
            .map(|(&z0, xnn)| z0)
            .collect::<Vec<Fp>>(),
        ecfft::ecfft::Moiety::S1,
    );
    let z0_8_16 = z0_8_evals
        .iter()
        .zip(&z0_8_ext)
        .zip(&fftree_8.xnn_s)
        .flat_map(|((&v0, &v1), &xnn)| [v0, v1])
        .collect::<Vec<_>>();
    let z0_8_16_pow_4 = z0_8_16
        .iter()
        .zip(&s_16)
        .map(|(y, s)| (y - &s.square().square()).square().square())
        .collect::<Vec<_>>();
    println!("{:?}", lagrange_interpolate(&s_16, &z0_8_16));
    println!("{:?}", lagrange_interpolate(&s_16, &z0_8_16_pow_4));
    println!("{:?}", z0_8_16_pow_4);
    // println!("{:?}", )

    // println!(
    //     "evals: {:?}",
    //     lagrange_interpolate(
    //         &s0_8,
    //         &fftree_8
    //             .z0_s1
    //             .iter()
    //             .enumerate()
    //             .map(|(i, v)| v.square() / fftree_8.xnn_s[i * 2])
    //             .collect::<Vec<_>>()
    //     )
    // );
}
