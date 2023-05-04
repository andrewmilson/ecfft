use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use ecfft::secp256k1::Fp;
use ecfft::FftreeField;
use std::ops::Deref;
use std::time::Instant;

fn main() {
    let n = 1 << 14;
    let mut rng = rand::thread_rng();

    let now = Instant::now();
    let fftree = Fp::build_fftree(n).unwrap();
    println!("FFTree generation time: {:?}", now.elapsed());

    let poly = DensePolynomial::rand(n - 1, &mut rng);

    let now = Instant::now();
    let ecfft_evals = fftree.enter(&poly);
    println!("evaluation time (fft): {:?}", now.elapsed());

    let now = Instant::now();
    let xs = fftree.subtree_with_size(n).f.get_layer(0);
    let ys = xs.iter().map(|x| poly.evaluate(x)).collect::<Vec<_>>();
    assert_eq!(ecfft_evals, ys);
    println!("naive O(n^2) eval: {:?}", now.elapsed());

    let now = Instant::now();
    let ecfft_coeffs = fftree.exit(&ecfft_evals);
    println!("interpolation time (ifft): {:?}", now.elapsed());

    assert_eq!(poly.deref(), ecfft_coeffs);
}
