#![feature(array_chunks, let_chains)]

extern crate alloc;

pub mod ec;
pub mod errors;
pub mod fftree;
pub mod find_curve;
pub mod utils;

use ark_ff::PrimeField;
use ec::Point;
pub use fftree::FFTree;
pub use fftree::Moiety;

/// The interface for fields that are able to be used in ECFFTs.
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

            // Curve with 2^31 | #E
            let curve = GoodCurve::new_odd(
                F!("68841697511640215656821373706374234598217962241974318227639791148600733690420"),
                F!("114313710132054377775997502842410245686288156214174113961380608383562274249190"),
            );
            let coset_offset = Point::new(
                F!("77457053793169897889304038652609345026746378478577985075189525651767568940193"),
                F!("35008124385965218187120632765572629524990942960751417275802166882882266399874"),
                curve,
            );
            let subgroup_generator = Point::new(
                F!("47338370574030999257492065562148397416925214763814892326587854834443644679284"),
                F!("26548029037302126578093317363856407980656459126501265688960367566322417565573"),
                curve,
            );
            let subgroup_two_addicity = 31;

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

            // println!("EVAL DOMAIN: {:?}", fftree.subtree_with_size(8).eval_domain());
            // // println!("EVAL DOMAIN: {:?}", fftree.z0_s1_inv);
            // let coeffs2 = fftree.exit(&evals);
            // // println!("TOOOO: {:?}", evals);
            // // println!("TOOOO: {:?}", coeffs2);

            let degree = fftree.degree(&evals);

            let poly = DensePolynomial::from_coefficients_slice(coeffs);
            assert_eq!(poly.degree(), degree);
        }
    }
}
