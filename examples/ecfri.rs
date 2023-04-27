use ark_ff::Field;
use ark_ff::One;
use ark_ff::Zero;
use ark_poly::univariate::DenseOrSparsePolynomial;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use ecfft::m31::Fp;
use ecfft::EcFftField;
use std::env;

fn main() {
    let fftree = Fp::build_fftree(1 << 5).unwrap();
    let fftree_8 = fftree.subtree_with_size(8);
}
