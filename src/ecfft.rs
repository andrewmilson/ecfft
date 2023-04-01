use crate::ec::Isogeny;
use crate::ec::Point;
use crate::ffforrest::Mat2x2;
use crate::utils::BinaryTree;
// use crate::utils::BinaryTree;
// use crate::utils::Mat2x2;
use ark_ff::PrimeField;
use ark_ff::Zero;
use ark_poly::Polynomial;
use std::iter::zip;

#[derive(Clone, Debug)]
pub struct FFTree<F: PrimeField> {
    f: BinaryTree<F>,
    recomp_matrices: BinaryTree<Mat2x2<F>>,
    decomp_matrices: BinaryTree<Mat2x2<F>>,
    isogenies: Vec<Isogeny<F>>,
    subtree: Option<Box<Self>>,
}

impl<F: PrimeField> FFTree<F> {
    pub fn new(coset_offset: Point<F>, generator: Point<F>) -> Self {
        assert_ne!(coset_offset, generator);
        assert_eq!(coset_offset.curve, generator.curve);
        let log_n = two_adicity(generator).expect("generator must have order 2^k");
        assert!(log_n < 32, "generator has unsupported two adicity");
        let n = 1 << log_n;

        let mut l = vec![F::zero(); n];
        let mut acc = Point::zero();
        for x in &mut l {
            *x = (coset_offset + acc).x;
            acc += generator;
        }

        let mut isogenies = Vec::new();
        let mut g = generator;
        for _ in 0..log_n {
            let isogeny = g
                .curve
                .unwrap()
                .two_isogenies()
                .into_iter()
                .find_map(|isogeny| {
                    let g_prime = isogeny.map(&g);
                    println!("YO: {:?} {:?}", two_adicity(g_prime), two_adicity(g));
                    if two_adicity(g_prime)? == two_adicity(g)? - 1 {
                        g = g_prime;
                        Some(isogeny)
                    } else {
                        None
                    }
                })
                .expect("cannot find a suitable isogeny");
            isogenies.push(isogeny);
        }

        // copy leaf nodes
        let mut f = BinaryTree::from(vec![F::zero(); 2 * n]);
        f[n..].copy_from_slice(&l);

        // generate internal nodes
        // TODO: use array_windows_mut
        let mut f_layers = f.get_layers_mut();
        for (i, isogeny) in isogenies.iter().enumerate() {
            let (prev_layer, layer) = {
                let (prev_layers, layers) = f_layers.split_at_mut(i + 1);
                (prev_layers.last_mut().unwrap(), layers.first_mut().unwrap())
            };

            let layer_size = layer.len();
            for (i, s) in layer.iter_mut().enumerate() {
                *s = isogeny.map_x(&prev_layer[i]).unwrap();
                debug_assert_eq!(*s, isogeny.map_x(&prev_layer[i + layer_size]).unwrap());
            }
        }

        Self::from_tree(f, isogenies)
    }

    fn from_tree(f: BinaryTree<F>, isogenies: Vec<Isogeny<F>>) -> Self {
        let mut f_layers = f.get_layers();
        let n = f.leaves().len();

        // Generate polynomial decomposition matrices
        // Lemma 3.2 (M_t) https://arxiv.org/abs/2107.08473
        // TODO: change notation
        let mut recomp_matrices = BinaryTree::from(vec![Mat2x2::identity(); n / 2]);
        let mut decomp_matrices = BinaryTree::from(vec![Mat2x2::identity(); n / 2]);
        let recomp_layers = recomp_matrices.get_layers_mut();
        let decomp_layers = decomp_matrices.get_layers_mut();
        for ((recomp_layer, decomp_layer), (l, iso)) in
            zip(zip(recomp_layers, decomp_layers), zip(f_layers, &isogenies))
        {
            let d = l.len() / 2;
            let v = &iso.x_denominator_map;
            for (i, (rmat, dmat)) in zip(recomp_layer, decomp_layer).enumerate() {
                let s0 = l[i * 2];
                let s1 = l[i * 2 + d];
                let v0 = v.evaluate(&s0).pow([(d / 2 - 1) as u64]);
                let v1 = v.evaluate(&s1).pow([(d / 2 - 1) as u64]);
                *rmat = Mat2x2([[v0, s0 * v0], [v1, s1 * v1]]).inverse().unwrap();

                let s0_prime = l[i * 2 + 1];
                let s1_prime = l[i * 2 + 1 + d];
                let v0_prime = v.evaluate(&s0_prime).pow([(d / 2 - 1) as u64]);
                let v1_prime = v.evaluate(&s1_prime).pow([(d / 2 - 1) as u64]);
                *dmat = Mat2x2([
                    [v0_prime, s0_prime * v0_prime],
                    [v1_prime, s1_prime * v1_prime],
                ]);
            }
        }

        let mut tree = Self {
            f,
            recomp_matrices,
            decomp_matrices,
            isogenies,
            subtree: None,
        };
        tree.subtree = tree.derive_subtree().map(Box::new);
        tree
    }

    fn derive_subtree(&self) -> Option<Self> {
        let Self { isogenies, f, .. } = self;
        let n = f.leaves().len() / 2;
        if n == 0 {
            return None;
        }

        let mut f_prime = BinaryTree::from(vec![F::zero(); n * 2]);
        let f_prime_layers = f_prime.get_layers_mut();
        let f_layers = f.get_layers();
        for (l_prime, l) in zip(f_prime_layers, f_layers) {
            // TODO: confusing notation compared to paper.
            for (s_prime, s) in l_prime.iter_mut().zip(l.iter().step_by(2)) {
                *s_prime = *s;
            }
        }

        Some(Self::from_tree(f_prime, isogenies.clone()))
    }

    // pub fn extend(&self, evals: &[F]) -> Vec<F> {
    //     // if self.f.leaves().len()

    //     if evals.len() == 1 {
    //         return evals.to_vec();
    //     }

    //     let n = evals.len();
    //     let layer_id = self.num_layers() - 1 - n.ilog2() as usize;
    //     let layer = self.get_layer(layer_id);

    //     // π_0 and π_1
    //     let (evals0, evals1): (Vec<F>, Vec<F>) = layer
    //         .recomp_matrices
    //         .iter()
    //         .enumerate()
    //         .map(|(i, m)| {
    //             let res = m * &[evals[i], evals[i + n / 2]];
    //             (res[0], res[1])
    //         })
    //         .unzip();

    //     let evals0_prime = self.extend(&evals0);
    //     let evals1_prime = self.extend(&evals1);

    //     // π_0' and π_1'
    //     let (lhs, rhs): (Vec<F>, Vec<F>) = layer
    //         .decomp_matrices
    //         .iter()
    //         .enumerate()
    //         .map(|(i, m)| {
    //             let res = m * &[evals0_prime[i], evals1_prime[i]];
    //             (res[0], res[1])
    //         })
    //         .unzip();

    //     [lhs, rhs].concat()
    // }
}

/// Returns the two adicity of a point i.e. returns `n` such that 2^n * p = 0.
/// Returns `None` if `p` isn't a point of order 2^n.
fn two_adicity<F: PrimeField>(p: Point<F>) -> Option<u32> {
    let mut acc = p;
    for i in 0..4096 {
        if acc.is_zero() {
            return Some(i);
        }
        acc += acc;
    }
    None
}
