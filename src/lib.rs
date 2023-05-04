#![feature(array_chunks, let_chains)]

extern crate alloc;

pub mod ec;
pub mod fftree;
pub mod find_curve;
pub mod utils;

use ark_ff::PrimeField;
use ec::Point;
pub use fftree::FFTree;
pub use fftree::Moiety;

/// The interface for fields to build FFTrees.
pub trait FftreeField: PrimeField {
    fn build_fftree(n: usize) -> Option<FFTree<Self>>;
}

pub mod secp256k1 {
    use super::FFTree;
    use super::FftreeField;
    use super::Point;
    use crate::ec::find_isogeny_chain;
    use crate::ec::GoodCurve;
    use ark_ff::Fp256;
    use ark_ff::MontBackend;
    use ark_ff::MontConfig;
    use ark_ff::MontFp as F;
    use ark_ff::Zero;

    /// Secp256k1 field
    #[derive(MontConfig)]
    #[modulus = "115792089237316195423570985008687907853269984665640564039457584007908834671663"]
    #[generator = "3"]
    #[small_subgroup_base = "3"]
    #[small_subgroup_power = "1"]
    pub struct FqConfig;
    pub type Fp = Fp256<MontBackend<FqConfig, 4>>;

    impl FftreeField for Fp {
        fn build_fftree(n: usize) -> Option<FFTree<Self>> {
            assert!(n.is_power_of_two());
            let log_n = n.ilog2();

            // Curve with 2^36 | #E
            let curve = GoodCurve::new_odd(
                F!("31172306031375832341232376275243462303334845584808513005362718476441963632613"),
                F!("45508371059383884471556188660911097844526467659576498497548207627741160623272"),
            );
            let coset_offset = Point::new(
                F!("105623886150579165427389078198493427091405550492761682382732004625374789850161"),
                F!("7709812624542158994629670452026922591039826164720902911013234773380889499231"),
                curve,
            );
            let subgroup_generator = Point::new(
                F!("41293412487153066667050767300223451435019201659857889215769525847559135483332"),
                F!("73754924733368840065089190002333366411120578552679996887076912271884749237510"),
                curve,
            );
            let subgroup_two_addicity = 36;

            // FFTree size is too large for our generator
            if log_n > subgroup_two_addicity {
                return None;
            }

            // get a generator of a subgroup with order `n`
            let mut generator = subgroup_generator;
            for _ in 0..subgroup_two_addicity - log_n {
                generator += generator;
            }

            // generate the FFTree leaf nodes
            let mut leaves = vec![Self::zero(); n];
            let mut acc = Point::zero();
            for x in &mut leaves {
                *x = (coset_offset + acc).x;
                acc += generator;
            }

            let isogenies = find_isogeny_chain(generator);
            let rational_maps = isogenies.into_iter().map(|isogeny| isogeny.r).collect();

            Some(FFTree::new(leaves, rational_maps))
        }
    }

    #[cfg(test)]
    mod tests {
        use super::FftreeField;
        use super::Fp;
        use crate::FFTree;
        use crate::Moiety;
        use ark_poly::univariate::DensePolynomial;
        use ark_poly::DenseUVPolynomial;
        use ark_poly::Polynomial;
        use ark_serialize::CanonicalDeserialize;
        use ark_serialize::CanonicalSerialize;
        use rand::rngs::StdRng;
        use rand::SeedableRng;
        use std::sync::OnceLock;

        static FFTREE: OnceLock<FFTree<Fp>> = OnceLock::new();

        fn get_fftree() -> &'static FFTree<Fp> {
            FFTREE.get_or_init(|| Fp::build_fftree(64).unwrap())
        }

        #[test]
        fn evaluates_polynomial() {
            let n = 64;
            let fftree = get_fftree();
            let mut rng = StdRng::from_seed([1; 32]);
            let poly = DensePolynomial::rand(n - 1, &mut rng);
            let eval_domain = fftree.subtree_with_size(n).eval_domain();

            let ecfft_evals = fftree.enter(&poly);

            let expected_evals: Vec<Fp> = eval_domain.iter().map(|x| poly.evaluate(x)).collect();
            assert_eq!(expected_evals, ecfft_evals);
        }

        #[test]
        fn extends_evaluations_from_s0_to_s1() {
            let n = 64;
            let fftree = get_fftree();
            let eval_domain = fftree.subtree_with_size(n).eval_domain();
            let mut rng = StdRng::from_seed([1; 32]);
            let poly = DensePolynomial::rand(n / 2 - 1, &mut rng);
            let (s0, s1): (Vec<Fp>, Vec<Fp>) = eval_domain.chunks(2).map(|s| (s[0], s[1])).unzip();
            let s0_evals: Vec<Fp> = s0.iter().map(|x| poly.evaluate(x)).collect();

            let s1_evals_actual = fftree.extend(&s0_evals, Moiety::S1);

            let s1_evals_expected: Vec<Fp> = s1.iter().map(|x| poly.evaluate(x)).collect();
            assert_eq!(s1_evals_expected, s1_evals_actual)
        }

        #[test]
        fn extends_evaluations_from_s1_to_s0() {
            let n = 64;
            let fftree = get_fftree();
            let eval_domain = fftree.subtree_with_size(n).eval_domain();
            let mut rng = StdRng::from_seed([1; 32]);
            let poly = DensePolynomial::rand(n / 2 - 1, &mut rng);
            let (s0, s1): (Vec<Fp>, Vec<Fp>) = eval_domain.chunks(2).map(|c| (c[0], c[1])).unzip();
            let s1_evals: Vec<Fp> = s1.iter().map(|x| poly.evaluate(x)).collect();

            let s0_evals_actual = fftree.extend(&s1_evals, Moiety::S0);

            let s0_evals_expected: Vec<Fp> = s0.iter().map(|x| poly.evaluate(x)).collect();
            assert_eq!(s0_evals_expected, s0_evals_actual)
        }

        #[test]
        fn deserialized_uncompressed_tree_works() {
            let n = 64;
            let fftree = get_fftree();
            let mut rng = StdRng::from_seed([1; 32]);
            let poly = DensePolynomial::rand(n - 1, &mut rng);
            let eval_domain = fftree.subtree_with_size(n).eval_domain();

            let mut fftree_bytes = Vec::new();
            fftree.serialize_uncompressed(&mut fftree_bytes).unwrap();
            let fftree = FFTree::deserialize_uncompressed(&*fftree_bytes).unwrap();
            let ecfft_evals = fftree.enter(&poly);

            let expected_evals: Vec<Fp> = eval_domain.iter().map(|x| poly.evaluate(x)).collect();
            assert_eq!(expected_evals, ecfft_evals);
        }

        #[test]
        fn deserialized_compressed_tree_works() {
            let n = 64;
            let fftree = get_fftree();
            let mut rng = StdRng::from_seed([1; 32]);
            let poly = DensePolynomial::rand(n - 1, &mut rng);
            let eval_domain = fftree.subtree_with_size(n).eval_domain();

            let mut fftree_bytes = Vec::new();
            fftree.serialize_compressed(&mut fftree_bytes).unwrap();
            let fftree = FFTree::deserialize_compressed(&*fftree_bytes).unwrap();
            let ecfft_evals = fftree.enter(&poly);

            let expected_evals: Vec<Fp> = eval_domain.iter().map(|x| poly.evaluate(x)).collect();
            assert_eq!(expected_evals, ecfft_evals);
        }
    }
}

pub mod m31 {
    use super::FFTree;
    use super::FftreeField;
    use super::Point;
    use crate::ec::build_ec_fftree;
    use crate::ec::ShortWeierstrassCurve;
    pub use ark_ff_optimized::fp31::Fp;

    impl FftreeField for Fp {
        fn build_fftree(n: usize) -> Option<FFTree<Fp>> {
            /// Supersingular curve with 2^31 | #E
            const CURVE: ShortWeierstrassCurve<Fp> = ShortWeierstrassCurve::new(Fp(1), Fp(0));
            const COSET_OFFSET: Point<ShortWeierstrassCurve<Fp>> =
                Point::new(Fp(1048755163), Fp(279503108), CURVE);
            const SUBGROUP_GENERATOR: Point<ShortWeierstrassCurve<Fp>> =
                Point::new(Fp(1273083559), Fp(804329170), CURVE);
            const SUBGORUP_TWO_ADDICITY: u32 = 28;

            build_ec_fftree(
                SUBGROUP_GENERATOR,
                1 << SUBGORUP_TWO_ADDICITY,
                COSET_OFFSET,
                n,
            )
        }
    }

    // TODO: there's a lot of repetition between field tests. Should implement macro
    // or loop at solutions to remove duplication of test logic.
    #[cfg(test)]
    mod tests {
        use super::Fp;
        use crate::fftree::FFTree;
        use crate::FftreeField;
        use ark_ff::One;
        use ark_ff::Zero;
        use ark_poly::univariate::DensePolynomial;
        use ark_poly::DenseUVPolynomial;
        use ark_poly::Polynomial;
        use rand::rngs::StdRng;
        use rand::SeedableRng;
        use std::sync::OnceLock;

        static FFTREE: OnceLock<FFTree<Fp>> = OnceLock::new();

        fn get_fftree() -> &'static FFTree<Fp> {
            FFTREE.get_or_init(|| Fp::build_fftree(64).unwrap())
        }

        #[test]
        fn evaluates_polynomial() {
            let n = 64;
            let fftree = get_fftree();
            let mut rng = StdRng::from_seed([1; 32]);
            let poly = DensePolynomial::rand(n - 1, &mut rng);
            let eval_domain = fftree.subtree_with_size(n).eval_domain();

            let ecfft_evals = fftree.enter(&poly);

            let expected_evals: Vec<Fp> = eval_domain.iter().map(|x| poly.evaluate(x)).collect();
            assert_eq!(expected_evals, ecfft_evals);
        }

        #[test]
        fn interpolates_evaluations() {
            let fftree = get_fftree();
            let one = Fp::one();
            let zero = Fp::zero();
            let coeffs: &[Fp] = &[one, one, Fp::from(5u8), zero, zero, one, zero, zero];
            let evals = fftree.enter(coeffs);

            let exit_coeffs = fftree.exit(&evals);

            assert_eq!(coeffs, &exit_coeffs);
        }

        #[test]
        fn determines_degree() {
            let fftree = get_fftree();
            let one = Fp::one();
            let zero = Fp::zero();
            let coeffs = &[one, one, one, zero, zero, one, zero, zero];
            let evals = fftree.enter(coeffs);

            let degree = fftree.degree(&evals);

            let poly = DensePolynomial::from_coefficients_slice(coeffs);
            assert_eq!(poly.degree(), degree);
        }
    }
}
