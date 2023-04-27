use ark_ff::FftField;
use ark_ff::Field;
use ark_ff::One;
use ark_ff::Zero;
use ark_ff_optimized::fp64::Fp;
use ark_poly::univariate::DenseOrSparsePolynomial;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::Polynomial;
use ark_poly::Radix2EvaluationDomain;

fn main() {
    let f = DensePolynomial::from_coefficients_vec(vec![
        Fp::from(6u32),
        Fp::from(2u32),
        Fp::from(9u32),
        Fp::from(7u32),
        Fp::from(7u32),
        Fp::from(7u32),
        Fp::from(7u32),
        Fp::from(7u32),
    ]);

    let res = transpose(f.clone().coeffs, 2);
    let g = res[0].clone();
    let h = res[1].clone();

    let half_domain = Radix2EvaluationDomain::<Fp>::new(f.len() / 2).unwrap();
    let full_domain = Radix2EvaluationDomain::<Fp>::new(f.len()).unwrap();
    println!("g_evals: {:?}", half_domain.fft(&g));
    println!(
        "h_evals: {:?}",
        half_domain
            .fft(&h)
            .into_iter()
            .zip(full_domain.elements())
            .map(|(e, x)| e * x)
            .collect::<Vec<_>>()
    );

    let domain = Radix2EvaluationDomain::new(f.len()).unwrap();
    let mut evals = f.evaluate_over_domain(domain).evals;
    println!("evals: {:?}", evals);

    let mut transposed_evals = transpose(evals, 4);
    let domain2 = Radix2EvaluationDomain::<Fp>::new(2).unwrap();
    let t_res = transposed_evals
        .into_iter()
        .map(|c| domain2.ifft(&c))
        .collect::<Vec<_>>()
        .concat();
    println!("YO: {:?}", transpose(t_res, 2));

    let f = DensePolynomial::from_coefficients_vec(vec![
        Fp::from(5u32),
        Fp::from(1u32),
        Fp::from(5u32),
        Fp::from(1u32),
    ]);
    let g = DensePolynomial::from_coefficients_vec(vec![Fp::from(2u32), Fp::from(3u32)]);
    let div_res = DenseOrSparsePolynomial::divide_with_q_and_r(&(&f).into(), &(&g).into()).unwrap();
    println!("Quotient: {:?}", div_res.0);
    println!("Remainder: {:?}", div_res.1);

    let f0 = DensePolynomial::from_coefficients_vec(vec![f[0], f[2]]);
    let f1 = DensePolynomial::from_coefficients_vec(vec![f[1], f[3]]);
    let g0 = DensePolynomial::from_coefficients_vec(vec![g[0], g[2]]);
    let g1 = DensePolynomial::from_coefficients_vec(vec![g[1], g[3]]);

    // evals.chunks_mut(2).for_each(|c| ifft_in_place(domain2, c));
    // println!("evals: {:?}", evals);

    // let g = DensePolynomial::from_coefficients_vec(vec![f[0], f[2]]);
    // let h = DensePolynomial::from_coefficients_vec(vec![f[1], f[3]]);

    // let domain = Radix2EvaluationDomain::new(f.len()).unwrap();
    // let evals = f.evaluate_over_domain(domain).evals;
    // println!("evals: {:?}", evals);

    // let half_domain = Radix2EvaluationDomain::<Fp>::new(2).unwrap();
    // let coeffs_a = half_domain.ifft(&vec![evals[0], evals[2]]);
    // let coeffs_b = half_domain.ifft(&vec![evals[1], evals[3]]);
    // println!("CA: {:?}", coeffs_a);
    // println!("CB: {:?}", coeffs_b);

    // let g_evals = g.evaluate_over_domain(half_domain).evals;
    // let h_evals = h.evaluate_over_domain(half_domain).evals;
    // println!("g_evals: {:?}", g_evals);
    // println!("h_evals: {:?}", h_evals);
}

fn ifft_in_place<F: FftField>(domain: Radix2EvaluationDomain<F>, evals: &mut [F]) {
    let coeffs = domain.ifft(&evals);
    evals.iter_mut().zip(coeffs).for_each(|(mut e, c)| *e = c);
}

fn transpose<F: Field>(coefficients: Vec<F>, num_columns: usize) -> Vec<Vec<F>> {
    let column_len = coefficients.len() / num_columns;

    let mut result = (0..num_columns)
        .map(|_| vec![F::zero(); column_len])
        .collect::<Vec<_>>();

    // TODO: implement multi-threaded version
    for (i, coeff) in coefficients.into_iter().enumerate() {
        let row_idx = i / num_columns;
        let col_idx = i % num_columns;
        result[col_idx][row_idx] = coeff;
    }

    result
}

// pub fn transpose_slice<T: Copy + Send + Sync + Default, const N: usize>(
//     source: &[T],
// ) -> Vec<[T; N]> {
//     let row_count = source.len() / N;
//     assert_eq!(
//         row_count * N,
//         source.len(),
//         "source length must be divisible by {}, but was {}",
//         N,
//         source.len()
//     );

//     let mut result = group_vector_elements(vec![T::default(); row_count *
// N]);     result.iter().enumerate().for_each(|(i, element)| {
//         for j in 0..N {
//             element[j] = source[i + j * row_count]
//         }
//     });
//     result
// }

// pub fn group_vector_elements<T, const N: usize>(source: Vec<T>) -> Vec<[T;
// N]> {     assert_eq!(
//         source.len() % N,
//         0,
//         "source length must be divisible by {}, but was {}",
//         N,
//         source.len()
//     );
//     let mut v = core::mem::ManuallyDrop::new(source);
//     let p = v.as_mut_ptr();
//     let len = v.len() / N;
//     let cap = v.capacity() / N;
//     unsafe { Vec::from_raw_parts(p as *mut [T; N], len, cap) }
// }

// f(x) = 1x^3 + 5x^2 + 1x^1 + 5x^0
// g(x) = 0x^3 + 0x^2 + 3x^1 + 2x^0
// f / g =  quotient: 12297829379609722881x^2 + 16397105839479630509x^1 +
// 7515340176428163982          remainder: 3416063716558256362
//
// f_0(x) = 5x^1 + 5x^0
// f_1(x) = 1x^1 + 1x^0
// f(x) = f_0(x^2) + x * f_1(x^2)
//
// g_0(x) = 0x^1 + 2x^0
// g_1(x) = 0x^1 + 3x^0
// g(x) = g_0(x^2) + x * g_1(x^2)
//
// f mod g = (f_0 mod g_0) + x * (f_1 mod g_1)

// FFTree size 2
// =============
// s = [a, b]
// s0 = [a]
// s1 = [b]
// x^1 = [a, b]
// z0 = (x - a) = [0, b - a]
// z0 mod x = z0 - x = [-a, -a]
// (z0 mod x)^2 = [-a^2, -a^2]
// z1^2 mod x = [-b^2, -b^2]
//
// FFTree size 4
// =============
// s = [a, c, b, d]
// s0 = [a, b]
// s1 = [c, d]
// x^2 = [a^2, c^2, b^2, d^2]
// z0 = (x - a)(x - b) = z0_2 * z1_2 = x^2 + -(a + b)*x + ab
// z0^2 = x^4 + (a + b)*(a + b)x^2 + ab^2
// aa + 2ab + bb
//
// z0^2 mod x^2 = (z0_2^2 mod x) * (z1_2^2 mod x)
// z1^2 mod x^2 = MOD( z0^2 mod x^2, )
