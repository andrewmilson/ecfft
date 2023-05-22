use crate::utils::BinaryTree;
use crate::utils::Mat2x2;
use alloc::boxed::Box;
use ark_ff::batch_inversion;
use ark_ff::vec;
use ark_ff::vec::Vec;
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use ark_serialize::Compress;
use ark_serialize::Valid;
use core::cmp::Ordering;
use core::iter::zip;

#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct RationalMap<F: Field> {
    pub numerator: DensePolynomial<F>,
    pub denominator: DensePolynomial<F>,
}

impl<F: Field> RationalMap<F> {
    pub fn new(numerator: &[F], denominator: &[F]) -> Self {
        let numerator = DensePolynomial::from_coefficients_slice(numerator);
        let denominator = DensePolynomial::from_coefficients_slice(denominator);
        Self {
            numerator,
            denominator,
        }
    }

    pub fn map(&self, x: &F) -> Option<F> {
        Some(self.numerator.evaluate(x) * self.denominator.evaluate(x).inverse()?)
    }

    pub fn zero() -> Self {
        Self::new(&[], &[F::one()])
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Moiety {
    S0,
    S1,
}

#[derive(Clone, Debug)]
pub struct FFTree<F: Field> {
    pub f: BinaryTree<F>,
    pub recombine_matrices: BinaryTree<Mat2x2<F>>,
    pub decompose_matrices: BinaryTree<Mat2x2<F>>,
    pub rational_maps: Vec<RationalMap<F>>,
    pub subtree: Option<Box<Self>>,
    pub xnn_s: Vec<F>,          // = <X^(n/2) ≀ S>
    pub xnn_s_inv: Vec<F>,      // = <1/X^(n/2) ≀ S>
    pub z0_s1: Vec<F>,          // = <Z_0 ≀ S_1>
    pub z1_s0: Vec<F>,          // = <Z_1 ≀ S_0>
    pub z0_inv_s1: Vec<F>,      // = <1/Z_0 ≀ S_1>
    pub z1_inv_s0: Vec<F>,      // = <1/Z_1 ≀ S_0>
    pub z0z0_rem_xnn_s: Vec<F>, // = <Z_0^2 mod X^(n/2) ≀ S>
    pub z1z1_rem_xnn_s: Vec<F>, // = <Z_0^2 mod X^(n/2) ≀ S>
}

// TODO: errors
impl<F: Field> FFTree<F> {
    pub fn new(leaves: Vec<F>, rational_maps: Vec<RationalMap<F>>) -> Self {
        let n = leaves.len();
        assert!(n.is_power_of_two());
        let log_n = n.ilog2();
        assert_eq!(log_n, rational_maps.len() as u32);
        let n = 1 << log_n;

        // copy leaf nodes
        let mut f = BinaryTree::from(vec![F::zero(); 2 * n]);
        f[n..].copy_from_slice(&leaves);

        // generate internal nodes
        // TODO: would be cool to use array_windows_mut
        let mut f_layers = f.get_layers_mut();
        for (i, rational_map) in rational_maps.iter().enumerate() {
            let (prev_layer, layer) = {
                let (prev_layers, layers) = f_layers.split_at_mut(i + 1);
                (prev_layers.last_mut().unwrap(), layers.first_mut().unwrap())
            };

            let layer_size = layer.len();
            for (i, s) in layer.iter_mut().enumerate() {
                *s = rational_map.map(&prev_layer[i]).unwrap();
                debug_assert_eq!(*s, rational_map.map(&prev_layer[i + layer_size]).unwrap());
            }
        }

        Self::from_tree(f, rational_maps)
    }

    fn extend_impl(&self, evals: &[F], moiety: Moiety) -> Vec<F> {
        let n = evals.len();
        if n == 1 {
            return evals.to_vec();
        }

        let layer = (self.f.num_layers() - 2 - n.ilog2()) as usize;

        // π_0 and π_1
        let mut evals0 = vec![F::zero(); n / 2];
        let mut evals1 = vec![F::zero(); n / 2];
        for (i, m) in self
            .decompose_matrices
            .get_layer(layer)
            .iter()
            .skip(match moiety {
                Moiety::S0 => 1,
                Moiety::S1 => 0,
            })
            .step_by(2)
            .enumerate()
        {
            let [v0, v1] = m * &[evals[i], evals[i + n / 2]];
            evals0[i] = v0;
            evals1[i] = v1;
        }

        // π_0' and π_1'
        let evals0_prime = self.extend_impl(&evals0, moiety);
        let evals1_prime = self.extend_impl(&evals1, moiety);

        let mut res = vec![F::zero(); n];
        for (i, m) in self
            .recombine_matrices
            .get_layer(layer)
            .iter()
            .skip(match moiety {
                Moiety::S0 => 0,
                Moiety::S1 => 1,
            })
            .step_by(2)
            .enumerate()
        {
            let [v0, v1] = m * &[evals0_prime[i], evals1_prime[i]];
            res[i] = v0;
            res[i + n / 2] = v1;
        }
        res
    }

    /// Extends evals on the chosen moiety
    pub fn extend(&self, evals: &[F], moiety: Moiety) -> Vec<F> {
        let tree = self.subtree_with_size(evals.len() * 2);
        tree.extend_impl(evals, moiety)
    }

    fn mextend_impl(&self, evals: &[F], moiety: Moiety) -> Vec<F> {
        let e = self.extend_impl(evals, moiety);
        let z = match moiety {
            Moiety::S1 => &self.z0_s1,
            Moiety::S0 => &self.z1_s0,
        };
        zip(e, z).map(|(e, z)| e + z).collect()
    }

    /// Extends special monic polynomials
    pub fn mextend(&self, evals: &[F], moiety: Moiety) -> Vec<F> {
        let tree = self.subtree_with_size(evals.len() * 2);
        tree.mextend_impl(evals, moiety)
    }

    fn enter_impl(&self, coeffs: &[F]) -> Vec<F> {
        let n = coeffs.len();
        if n == 1 {
            return coeffs.to_vec();
        }

        let subtree = self.subtree().unwrap();
        let u0 = subtree.enter(&coeffs[0..n / 2]);
        let v0 = subtree.enter(&coeffs[n / 2..]);
        let u1 = self.extend(&u0, Moiety::S1);
        let v1 = self.extend(&v0, Moiety::S1);

        let mut res = Vec::new();
        for i in 0..n / 2 {
            res.push(u0[i] + v0[i] * self.xnn_s[i * 2]);
            res.push(u1[i] + v1[i] * self.xnn_s[i * 2 + 1]);
        }
        res
    }

    /// Converts from coefficient to evaluation representation of a polynomial
    pub fn enter(&self, coeffs: &[F]) -> Vec<F> {
        let tree = self.subtree_with_size(coeffs.len());
        tree.enter_impl(coeffs)
    }

    fn degree_impl(&self, evals: &[F]) -> usize {
        let n = evals.len();
        if n == 1 {
            return 0;
        }

        let subtree = self.subtree().unwrap();
        let (e0, e1): (Vec<F>, Vec<F>) = evals.chunks(2).map(|e| (e[0], e[1])).unzip();

        // The intuition is that if `degree < n/2` then all coefficients on the RHS
        // will be 0. Therefore if `degree < n/2` then e1 == extend(e0)
        let g1 = self.extend_impl(&e0, Moiety::S1);
        if g1 == e1 {
            return subtree.degree_impl(&e0);
        }

        // compute `<(π-g)/Z_0 ≀ S1>`
        // isolate the evaluations of the coefficients on the RHS
        let t1: Vec<F> = zip(zip(e1, g1), &self.z0_inv_s1)
            .map(|((e1, g1), z0_inv)| (e1 - g1) * z0_inv)
            .collect();
        let t0 = self.extend_impl(&t1, Moiety::S0);
        n / 2 + subtree.degree_impl(&t0)
    }

    /// Evaluates the degree of an evaluation table in `O(n log n)`
    pub fn degree(&self, evals: &[F]) -> usize {
        let tree = self.subtree_with_size(evals.len());
        tree.degree_impl(evals)
    }

    pub fn exit_impl(&self, evals: &[F]) -> Vec<F> {
        let n = evals.len();
        if n == 1 {
            return evals.to_vec();
        }

        let u0: Vec<F> = self
            .modular_reduce_impl(evals, &self.xnn_s, &self.z0z0_rem_xnn_s)
            .into_iter()
            .step_by(2)
            .collect();

        let subtree = self.subtree().unwrap();
        let mut a = subtree.exit_impl(&u0);

        let xnn0_inv = self.xnn_s_inv.iter().step_by(2);
        let e0 = evals.iter().step_by(2);
        let v0: Vec<F> = zip(zip(e0, u0), xnn0_inv)
            .map(|((&e0, u0), xnn0_inv)| (e0 - u0) * xnn0_inv)
            .collect();
        let mut b = subtree.exit_impl(&v0);

        a.append(&mut b);
        a
    }

    /// Converts from evaluation to coefficient representation of a polynomial
    pub fn exit(&self, evals: &[F]) -> Vec<F> {
        let tree = self.subtree_with_size(evals.len());
        tree.exit_impl(evals)
    }

    fn redc_impl(&self, evals: &[F], a: &[F], moiety: Moiety) -> Vec<F> {
        let (e0, e1): (Vec<F>, Vec<F>) = evals.chunks(2).map(|e| (e[0], e[1])).unzip();
        let (mut a0_inv, a1): (Vec<F>, Vec<F>) = a.chunks(2).map(|a| (a[0], a[1])).unzip();
        batch_inversion(&mut a0_inv);

        // compute <π/a ≀ S>
        let t0: Vec<F> = zip(e0, a0_inv).map(|(e0, a0_inv)| e0 * a0_inv).collect();
        let g1 = self.extend_impl(
            &t0,
            match moiety {
                Moiety::S1 => Moiety::S0,
                Moiety::S0 => Moiety::S1,
            },
        );

        let z_inv = match moiety {
            Moiety::S0 => &self.z0_inv_s1,
            Moiety::S1 => &self.z1_inv_s0,
        };

        // compute `<(π - a*g)/Z ≀ S'>`
        let h1: Vec<F> = zip(zip(e1, g1), zip(a1, z_inv))
            .map(|((e1, g1), (a1, z0_inv))| (e1 - g1 * a1) * z0_inv)
            .collect();
        let h0 = self.extend_impl(&h1, moiety);

        zip(h0, h1).flat_map(|h| [h.0, h.1]).collect()
    }

    /// Computes <P(X)*Z_0(x)^(-1) mod a ≀ S>
    /// Z_0 is the vanishing polynomial of S_0
    /// `a` must be a polynomial of max degree `n/2` having no zeroes in `S_0`
    pub fn redc_z0(&self, evals: &[F], a: &[F]) -> Vec<F> {
        let tree = self.subtree_with_size(evals.len());
        tree.redc_impl(evals, a, Moiety::S0)
    }

    /// Computes `<P(X)*Z_1(x)^(-1) mod A ≀ S>`
    /// `Z_1` is the vanishing polynomial of `S_1`
    /// `A` must be a polynomial of max degree `n/2` having no zeroes in `S_1`
    pub fn redc_z1(&self, evals: &[F], a: &[F]) -> Vec<F> {
        let tree = self.subtree_with_size(evals.len());
        tree.redc_impl(evals, a, Moiety::S1)
    }

    fn modular_reduce_impl(&self, evals: &[F], a: &[F], c: &[F]) -> Vec<F> {
        let h = self.redc_impl(evals, a, Moiety::S0);
        let hc: Vec<F> = zip(h, c).map(|(h, c)| h * c).collect();
        self.redc_impl(&hc, a, Moiety::S0)
    }

    /// Computes MOD algorithm
    /// `a` must be a polynomial of max degree `n/2` having no zeroes in `S_0`
    /// `c` must be the evaluation table `<Z_0^2 mod a ≀ S>`
    pub fn modular_reduce(&self, evals: &[F], a: &[F], c: &[F]) -> Vec<F> {
        let tree = self.subtree_with_size(evals.len());
        tree.modular_reduce_impl(evals, a, c)
    }

    fn vanish_impl(&self, vanish_domain: &[F]) -> Vec<F> {
        let n = vanish_domain.len();
        if n == 1 {
            let l = self.f.leaves();
            assert_eq!(2, l.len());
            let alpha = vanish_domain[0];
            return vec![alpha - l[0], alpha - l[1]];
        }

        let subtree = self.subtree().unwrap();
        let qp = subtree.vanish_impl(&vanish_domain[0..n / 2]);
        let qpp = subtree.vanish_impl(&vanish_domain[n / 2..n]);
        let q_s0: Vec<F> = zip(qp, qpp).map(|(qp, qpp)| qp * qpp).collect();
        let q_s1 = self.mextend(&q_s0, Moiety::S1);
        zip(q_s0, q_s1)
            .flat_map(|(q_s0, q_s1)| [q_s0, q_s1])
            .collect()
    }

    /// Returns an evaluation of the vanishing polynomial `Z(x) = ∏ (x - a_i)`
    /// Runtime `O(n log^2 n)`. `vanishi_domain = [a_0, a_1, ..., a_(n - 1)]`
    /// Section 7.1 https://arxiv.org/pdf/2107.08473.pdf
    pub fn vanish(&self, vanish_domain: &[F]) -> Vec<F> {
        let tree = self.subtree_with_size(vanish_domain.len() * 2);
        tree.vanish_impl(vanish_domain)
    }

    fn from_tree(f: BinaryTree<F>, rational_maps: Vec<RationalMap<F>>) -> Self {
        let subtree = Self::derive_subtree(&f, &rational_maps).map(Box::new);
        let f_layers = f.get_layers();
        let n = f.leaves().len();
        let nn = n as u64 / 2;
        let nnnn = n as u64 / 4;
        let s = f_layers[0];

        // Precompute eval table <X^(n/2) ≀ S>
        // TODO: compute xnn from xnnnn for n != 2
        let xnnnn_s: Vec<F> = f_layers[0].iter().map(|x| x.pow([nnnn])).collect();
        let mut xnnnn_s_inv = xnnnn_s.clone();
        batch_inversion(&mut xnnnn_s_inv);
        let xnn_s: Vec<F> = f_layers[0].iter().map(|x| x.pow([nn])).collect();
        let mut xnn_s_inv = xnn_s.clone();
        batch_inversion(&mut xnn_s_inv);

        // Split S into its two moieties S0 and S1
        let (s0, s1): (Vec<F>, Vec<F>) = s.chunks_exact(2).map(|s| (s[0], s[1])).unzip();

        // Generate polynomial decomposition matrices
        // Lemma 3.2 (M_t) https://arxiv.org/abs/2107.08473
        // TODO: change notation
        let mut recombine_matrices = BinaryTree::from(vec![Mat2x2::identity(); n]);
        let mut decompose_matrices = BinaryTree::from(vec![Mat2x2::identity(); n]);
        let recombine_layers = recombine_matrices.get_layers_mut();
        let decompose_layers = decompose_matrices.get_layers_mut();
        for ((recombine_layer, decompose_layer), (l, map)) in zip(
            zip(recombine_layers, decompose_layers),
            zip(f_layers, &rational_maps),
        ) {
            let d = l.len() / 2;
            if d == 1 {
                continue;
            }

            let v = &map.denominator;
            for (i, (rmat, dmat)) in zip(recombine_layer, decompose_layer).enumerate() {
                let s0 = l[i];
                let s1 = l[i + d];
                let v0 = v.evaluate(&s0).pow([(d / 2 - 1) as u64]);
                let v1 = v.evaluate(&s1).pow([(d / 2 - 1) as u64]);
                *rmat = Mat2x2([[v0, s0 * v0], [v1, s1 * v1]]);
                *dmat = rmat.inverse().unwrap();
            }
        }

        let mut tree = Self {
            f,
            recombine_matrices,
            decompose_matrices,
            rational_maps,
            subtree,
            xnn_s,
            xnn_s_inv,
            z0_s1: Vec::new(),
            z1_s0: Vec::new(),
            z1_inv_s0: Vec::new(),
            z0_inv_s1: Vec::new(),
            z0z0_rem_xnn_s: Vec::new(),
            z1z1_rem_xnn_s: Vec::new(),
        };

        // Precompute eval tables <Z_0 ≀ S_1> and <Z_1 ≀ S_0> using our partial FFTree
        // Z_0 is the vanishing polynomial of S_0 i.e. Z_0(x) = Π(x - s0_i)
        // TODO: this code is a little brittle, might be nice to find a cleaner solution
        match n.cmp(&2) {
            Ordering::Greater => {
                // compute z0_s1 in O(n log n) using the subtree's vanishing polynomials
                let zero = F::zero();
                let st = tree.subtree.as_ref().unwrap();
                let st_z0_s0: Vec<F> = st.z0_s1.iter().flat_map(|&y| [zero, y]).collect();
                let st_z1_s0: Vec<F> = st.z1_s0.iter().flat_map(|&y| [y, zero]).collect();
                let st_z0_s1 = tree.extend(&st_z0_s0, Moiety::S1);
                let st_z1_s1 = tree.extend(&st_z1_s0, Moiety::S1);
                tree.z0_s1 = zip(st_z0_s1, st_z1_s1).map(|(z0, z1)| z0 * z1).collect();

                // compute z1_s in O(n log^2 n) - .vanish() uses z0_s1
                let z1_s = tree.vanish(&s1);
                tree.z1_s0 = z1_s.into_iter().step_by(2).collect();
            }
            Ordering::Equal => {
                // base cases
                tree.z0_s1 = vec![s1[0] - s0[0]];
                tree.z1_s0 = vec![s0[0] - s1[0]];
            }
            Ordering::Less => {}
        }

        tree.z0_inv_s1 = tree.z0_s1.clone();
        tree.z1_inv_s0 = tree.z1_s0.clone();
        batch_inversion(&mut tree.z0_inv_s1);
        batch_inversion(&mut tree.z1_inv_s0);

        // Precompute evaluation tables <Z_0^2 rem X^(n/2) ≀ S> and
        // <Z_1^2 rem X^(n/2) ≀ S> using our partial FFTree.
        // TODO: alternative approach - https://www.math.toronto.edu/swastik/ECFFT2.pdf
        // section 5.1 equation 11 gives an algorithm for the vanishing polynomial.
        // Might be nice for a O(log n) verifier vanishing polynomial evaluation.
        match n.cmp(&2) {
            Ordering::Greater => {
                // compute z0z0_rem_xnn_s in O(n log n)
                let st = tree.subtree.as_ref().unwrap();
                let z0_rem_xnnnn_sq_s0 = zip(&st.z0z0_rem_xnn_s, &st.z1z1_rem_xnn_s)
                    .map(|(y0, y1)| *y0 * y1)
                    .collect::<Vec<F>>();
                let z0z0_rem_xnnnn_s0 =
                    st.modular_reduce(&z0_rem_xnnnn_sq_s0, &st.xnn_s, &st.z0z0_rem_xnn_s);
                let z0z0_rem_xnnnn_s1 = tree.extend(&z0z0_rem_xnnnn_s0, Moiety::S1);
                let z0z0_rem_xnnnn_s = zip(z0z0_rem_xnnnn_s0, z0z0_rem_xnnnn_s1)
                    .flat_map(|(y0, y1)| [y0, y1])
                    .collect::<Vec<F>>();
                let z0_s = tree.z0_s1.iter().flat_map(|&y1| [F::zero(), y1]);
                let z0_rem_xnn_s = zip(z0_s, &tree.xnn_s).map(|(z0, xnn)| z0 - xnn);
                let z0_rem_xnn_sq_s = z0_rem_xnn_s.map(|y| y.square()).collect::<Vec<F>>();
                let z0_rem_xnn_sq_div_xnnnn_s =
                    zip(&z0_rem_xnn_sq_s, zip(&z0z0_rem_xnnnn_s, &xnnnn_s_inv))
                        .map(|(z0_rem_xnn_sq, (z0z0_rem_xnnnn, xnnnn_inv))| {
                            (*z0_rem_xnn_sq - z0z0_rem_xnnnn) * xnnnn_inv
                        })
                        .collect::<Vec<F>>();
                let z0z0_div_xnnnn_rem_xnnnn_s =
                    tree.modular_reduce(&z0_rem_xnn_sq_div_xnnnn_s, &xnnnn_s, &z0z0_rem_xnnnn_s);
                tree.z0z0_rem_xnn_s =
                    zip(z0z0_rem_xnnnn_s, zip(z0z0_div_xnnnn_rem_xnnnn_s, xnnnn_s))
                        .map(|(z0z0_rem_xnnnn, (z0z0_div_xnnnn_rem_xnnnn, xnnnn))| {
                            z0z0_rem_xnnnn + xnnnn * z0z0_div_xnnnn_rem_xnnnn
                        })
                        .collect();

                // compute z1z1_rem_xnn_s in O(n log n)
                let z1_s = tree.z1_s0.iter().flat_map(|&y0| [y0, F::zero()]);
                let z1_rem_xnn_s = zip(z1_s, &tree.xnn_s).map(|(z1, xnn)| z1 - xnn);
                let z1z1 = z1_rem_xnn_s.map(|y| y.square()).collect::<Vec<F>>();
                tree.z1z1_rem_xnn_s = tree.modular_reduce(&z1z1, &tree.xnn_s, &tree.z0z0_rem_xnn_s);
            }
            Ordering::Equal => {
                // base cases
                tree.z0z0_rem_xnn_s = vec![s0[0].square(); 2];
                tree.z1z1_rem_xnn_s = vec![s1[0].square(); 2];
            }
            Ordering::Less => {}
        }

        tree
    }

    fn derive_subtree(f: &BinaryTree<F>, rational_maps: &[RationalMap<F>]) -> Option<Self> {
        let n = f.leaves().len() / 2;
        if n == 0 {
            return None;
        }

        let mut f_prime = BinaryTree::from(vec![F::zero(); n * 2]);
        let f_prime_layers = f_prime.get_layers_mut();
        let f_layers = f.get_layers();
        for (l_prime, l) in zip(f_prime_layers, f_layers) {
            for (s_prime, s) in l_prime.iter_mut().zip(l.iter().step_by(2)) {
                *s_prime = *s;
            }
        }

        let rational_maps = rational_maps.split_last().map_or(vec![], |i| i.1.to_vec());
        Some(Self::from_tree(f_prime, rational_maps))
    }

    pub(crate) fn subtree(&self) -> Option<&Self> {
        Some(self.subtree.as_ref()?)
    }

    /// Returns a FFTree with `n` leaves
    pub fn subtree_with_size(&self, n: usize) -> &Self {
        assert!(n.is_power_of_two());
        match usize::cmp(&n, &self.f.leaves().len()) {
            Ordering::Less => self.subtree().unwrap().subtree_with_size(n),
            Ordering::Equal => self,
            Ordering::Greater => panic!("FFTree is too small"),
        }
    }
}

#[cfg(test)]
impl<F: Field> FFTree<F> {
    // TODO: maybe remove
    pub(crate) fn eval_domain(&self) -> &[F] {
        self.f.leaves()
    }
}

// TODO: implement CanonicalSerialize for Box in arkworks
// TODO: Derive bug "error[E0275]" with recursive field subtree
// TODO: add compress version (with flags perhaps)
impl<F: Field> CanonicalSerialize for FFTree<F> {
    fn serialize_with_mode<W: ark_serialize::Write>(
        &self,
        mut writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        #[deny(unused_variables)]
        let Self {
            f,
            recombine_matrices,
            decompose_matrices,
            rational_maps,
            xnn_s,
            xnn_s_inv,
            z0_s1,
            z1_s0,
            z0_inv_s1,
            z1_inv_s0,
            z0z0_rem_xnn_s,
            z1z1_rem_xnn_s,
            subtree,
        } = self;
        f.serialize_with_mode(&mut writer, compress)?;
        recombine_matrices.serialize_with_mode(&mut writer, compress)?;
        decompose_matrices.serialize_with_mode(&mut writer, compress)?;
        rational_maps.serialize_with_mode(&mut writer, compress)?;
        xnn_s.serialize_with_mode(&mut writer, compress)?;
        z0_s1.serialize_with_mode(&mut writer, compress)?;
        z1_s0.serialize_with_mode(&mut writer, compress)?;
        if compress == ark_serialize::Compress::No {
            // Inverses are the cheapest to regenerate (1 mul per element)
            xnn_s_inv.serialize_with_mode(&mut writer, compress)?;
            z0_inv_s1.serialize_with_mode(&mut writer, compress)?;
            z1_inv_s0.serialize_with_mode(&mut writer, compress)?;
        }
        z0z0_rem_xnn_s.serialize_with_mode(&mut writer, compress)?;
        z1z1_rem_xnn_s.serialize_with_mode(&mut writer, compress)?;
        // TODO: get "error[E0275]: overflow evaluating the requirement" for:
        // subtree.as_ref().map(Box::as_ref).serialize_with_mode(...)
        (subtree.is_some()).serialize_with_mode(&mut writer, compress)?;
        if let Some(subtree) = subtree.as_ref().map(Box::as_ref) {
            subtree.serialize_with_mode(writer, compress)?;
        }
        Ok(())
    }

    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        #[deny(unused_variables)]
        let Self {
            f,
            recombine_matrices,
            decompose_matrices,
            rational_maps,
            xnn_s,
            xnn_s_inv,
            z0_s1,
            z1_s0,
            z0_inv_s1,
            z1_inv_s0,
            z0z0_rem_xnn_s,
            z1z1_rem_xnn_s,
            subtree,
        } = self;
        let mut size = f.serialized_size(compress)
            + recombine_matrices.serialized_size(compress)
            + decompose_matrices.serialized_size(compress)
            + rational_maps.serialized_size(compress)
            + xnn_s.serialized_size(compress)
            + z0_s1.serialized_size(compress)
            + z1_s0.serialized_size(compress)
            + z0z0_rem_xnn_s.serialized_size(compress)
            + z1z1_rem_xnn_s.serialized_size(compress)
            // subtree: 1 (for Option state) + subtree size
            + 1 + subtree.as_ref().map_or(0, |v| v.as_ref().serialized_size(compress));
        if compress == Compress::No {
            size += xnn_s_inv.serialized_size(compress)
                + z0_inv_s1.serialized_size(compress)
                + z1_inv_s0.serialized_size(compress);
        }
        size
    }
}

impl<F: Field> Valid for FFTree<F> {
    #[inline]
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        Ok(())
    }
}

// TODO: implement CanonicalDeserialize for Box in arkworks
// TODO: Derive bug "error[E0275]" with recursive field subtree
impl<F: Field> CanonicalDeserialize for FFTree<F> {
    fn deserialize_with_mode<R: ark_serialize::Read>(
        mut reader: R,
        compress: ark_serialize::Compress,
        validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        let f = BinaryTree::deserialize_with_mode(&mut reader, compress, validate)?;
        let recombine_matrices =
            BinaryTree::deserialize_with_mode(&mut reader, compress, validate)?;
        let decompose_matrices =
            BinaryTree::deserialize_with_mode(&mut reader, compress, validate)?;
        let rational_maps = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
        let xnn_s = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
        let z0_s1 = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
        let z1_s0 = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
        let mut xnn_s_inv: Vec<F>;
        let mut z0_inv_s1: Vec<F>;
        let mut z1_inv_s0: Vec<F>;
        match compress {
            Compress::Yes => {
                xnn_s_inv = xnn_s.clone();
                z0_inv_s1 = z0_s1.clone();
                z1_inv_s0 = z1_s0.clone();
                batch_inversion(&mut xnn_s_inv);
                batch_inversion(&mut z0_inv_s1);
                batch_inversion(&mut z1_inv_s0);
            }
            Compress::No => {
                xnn_s_inv = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
                z0_inv_s1 = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
                z1_inv_s0 = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
            }
        }
        let z0z0_rem_xnn_s = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
        let z1z1_rem_xnn_s = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
        let subtree = if bool::deserialize_with_mode(&mut reader, compress, validate)? {
            Some(Box::new(Self::deserialize_with_mode(
                reader, compress, validate,
            )?))
        } else {
            None
        };
        Ok(Self {
            f,
            recombine_matrices,
            decompose_matrices,
            rational_maps,
            subtree,
            xnn_s,
            xnn_s_inv,
            z0_s1,
            z1_s0,
            z0_inv_s1,
            z1_inv_s0,
            z0z0_rem_xnn_s,
            z1z1_rem_xnn_s,
        })
    }
}
