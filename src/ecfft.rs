use crate::ec::Isogeny;
use crate::ec::Point;
use crate::utils::BinaryTree;
use crate::utils::Mat2x2;
use ark_ff::batch_inversion;
use ark_ff::PrimeField;
use ark_ff::Zero;
use ark_poly::Polynomial;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use ark_serialize::Compress;
use ark_serialize::Valid;
use std::cmp::Ordering;
use std::iter::zip;

#[derive(Clone, Copy, Debug)]
pub enum Moiety {
    S0 = 0,
    S1 = 1,
}

struct FFTreeLayer<'a, F: PrimeField> {
    pub _l: &'a [F],
    pub _isogeny: &'a Isogeny<F>,
    pub decomp_matrices: &'a [Mat2x2<F>],
    pub recomp_matrices: &'a [Mat2x2<F>],
}

// TODO: consider:
// pub enum Capability {
//     ExtendS0,
//     ExtendS1,
//     Enter,
//     Exit,
//     Degree,
// }

#[derive(Clone, Debug)]
pub struct FFTree<F: PrimeField> {
    pub f: BinaryTree<F>,
    recomp_matrices: BinaryTree<Mat2x2<F>>,
    decomp_matrices: BinaryTree<Mat2x2<F>>,
    isogenies: Vec<Isogeny<F>>,
    subtree: Option<Box<Self>>,
    xnn_s: Vec<F>,         // = <X^(n/2) ≀ S>
    xnn_s_inv: Vec<F>,     // = <1/X^(n/2) ≀ S>
    pub z0_s1: Vec<F>,     // = <Z_0 ≀ S_1>
    pub z1_s0: Vec<F>,     // = <Z_1 ≀ S_0>
    pub z0_s1_inv: Vec<F>, // = <1/Z_0 ≀ S_1>
    pub z1_s0_inv: Vec<F>, // = <1/Z_1 ≀ S_0>
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

    fn extend_impl(&self, evals: &[F], moiety: Moiety) -> Vec<F> {
        let n = evals.len();
        if n == 1 {
            return evals.to_vec();
        }

        let i = self.f.num_layers() - 2 - n.ilog2();
        let layer = self.get_layer(i as usize);

        // π_0 and π_1
        let mut evals0 = vec![F::zero(); n / 2];
        let mut evals1 = vec![F::zero(); n / 2];
        for (i, m) in layer
            .recomp_matrices
            .iter()
            .skip(1 - moiety as usize)
            .step_by(2)
            .enumerate()
        {
            let v = m * &[evals[i], evals[i + n / 2]];
            evals0[i] = v[0];
            evals1[i] = v[1];
        }

        // π_0' and π_1'
        let evals0_prime = self.extend_impl(&evals0, moiety);
        let evals1_prime = self.extend_impl(&evals1, moiety);

        let mut res = vec![F::zero(); n];
        for (i, m) in layer
            .decomp_matrices
            .iter()
            .skip(moiety as usize)
            .step_by(2)
            .enumerate()
        {
            let v = m * &[evals0_prime[i], evals1_prime[i]];
            res[i] = v[0];
            res[i + n / 2] = v[1];
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

    /// Extends evals on the chosen moiety
    /// TODO: docs
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

        self.xnn_s
            .array_chunks()
            .enumerate()
            .flat_map(|(i, [s0, s1])| [u0[i] + v0[i] * s0, u1[i] + v1[i] * s1])
            .collect()
    }

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

        // The intuition is that if degree < n / 2 then all coefficients on the RHS
        // will be 0. Therefore if degree < n / 2 then e1 == extend(e0)
        let g1 = self.extend_impl(&e0, Moiety::S1);
        if g1 == e1 {
            return subtree.degree_impl(&e0);
        }

        // compute <(π-g)/Z_0 ≀ S1>
        // isolate the evaluations of the coefficients on the RHS
        let t1: Vec<F> = zip(zip(e1, g1), &self.z0_s1_inv)
            .map(|((e1, g1), z0_inv)| (e1 - g1) * z0_inv)
            .collect();
        let t0 = self.extend_impl(&t1, Moiety::S0);
        n / 2 + subtree.degree_impl(&t0)
    }

    /// Evaluates the degree of an evaluation table in O(n log n)
    pub fn degree(&self, evals: &[F]) -> usize {
        let tree = self.subtree_with_size(evals.len());
        tree.degree_impl(evals)
    }

    pub fn exit_impl(&self, evals: &[F]) -> Vec<F> {
        let n = evals.len();
        if n == 1 {
            return evals.to_vec();
        }

        println!("d{}", evals.len());
        let u0: Vec<F> = self
            .modular_reduce_impl(evals, &self.xnn_s)
            .array_chunks()
            .map(|[u0, _]| *u0)
            .collect();
        println!("u{}", u0.len());

        let subtree = self.subtree().unwrap();
        let mut a = subtree.exit_impl(&u0);

        let xnn0_inv = self.xnn_s_inv.array_chunks().map(|[xnn0, _]| xnn0);
        let e0 = evals.array_chunks().map(|[e0, _]| e0);
        let v0: Vec<F> = zip(zip(e0, u0), xnn0_inv)
            .map(|((&e0, u0), xnn0_inv)| (e0 - u0) * xnn0_inv)
            .collect();
        let mut b = subtree.exit_impl(&v0);

        println!("{} {}", a.len(), b.len());
        a.append(&mut b);
        a
    }

    pub fn exit(&self, evals: &[F]) -> Vec<F> {
        let tree = self.subtree_with_size(evals.len());
        tree.exit_impl(evals)
    }

    fn redc_z0_impl(&self, evals: &[F], a: &[F]) -> Vec<F> {
        let (e0, e1): (Vec<F>, Vec<F>) = evals.chunks(2).map(|e| (e[0], e[1])).unzip();
        let (mut a0_inv, a1): (Vec<F>, Vec<F>) = a.chunks(2).map(|a| (a[0], a[1])).unzip();
        batch_inversion(&mut a0_inv);

        // compute <π/a ≀ S0>
        let t0: Vec<F> = zip(e0, a0_inv).map(|(e0, a0_inv)| e0 * a0_inv).collect();
        let g1 = self.extend_impl(&t0, Moiety::S1);
        // compute <(π - a*g)/Z_0 ≀ S1>
        let h1: Vec<F> = zip(zip(e1, g1), zip(a1, &self.z0_s1_inv))
            .map(|((e1, g1), (a1, z0_inv))| (e1 - g1 * a1) * z0_inv)
            .collect();
        let h0 = self.extend_impl(&h1, Moiety::S0);

        zip(h0, h1).flat_map(|h| [h.0, h.1]).collect()
    }

    /// Computes <P(X)*Z_0(x)^(-1) mod a ≀ S>
    /// Z_0 is the vanishing polynomial of S_0
    /// `a` must be a polynomial of degree at most n/2 having no zeroes in S_0
    pub fn redc_z0(&self, evals: &[F], a: &[F]) -> Vec<F> {
        let tree = self.subtree_with_size(evals.len());
        tree.redc_z0_impl(evals, a)
    }

    fn redc_z1_impl(&self, evals: &[F], a: &[F]) -> Vec<F> {
        let (e0, e1): (Vec<F>, Vec<F>) = evals.chunks(2).map(|e| (e[0], e[1])).unzip();
        let (a0, mut a1_inv): (Vec<F>, Vec<F>) = a.chunks(2).map(|a| (a[0], a[1])).unzip();
        batch_inversion(&mut a1_inv);

        // compute <π/a ≀ S1>
        let t1: Vec<F> = zip(e1, a1_inv).map(|(e1, a1_inv)| e1 * a1_inv).collect();
        let g0 = self.extend_impl(&t1, Moiety::S0);
        // compute <(π - a*g)/Z_1 ≀ S1>
        let h0: Vec<F> = zip(zip(e0, g0), zip(a0, &self.z1_s0_inv))
            .map(|((e0, g0), (a0, z1_inv))| (e0 - g0 * a0) * z1_inv)
            .collect();
        let h1 = self.extend_impl(&h0, Moiety::S1);

        zip(h0, h1).flat_map(|h| [h.0, h.1]).collect()
    }

    /// Computes <P(X)*Z_1(x)^(-1) mod A ≀ S>
    /// Z_1 is the vanishing polynomial of S_1
    /// `A` must be a polynomial of degree at most n/2 having no zeroes in S_1
    pub fn redc_z1(&self, evals: &[F], a: &[F]) -> Vec<F> {
        let tree = self.subtree_with_size(evals.len());
        tree.redc_z1_impl(evals, a)
    }

    fn modular_reduce_impl(&self, evals: &[F], a: &[F]) -> Vec<F> {
        // TODO: precompute. surely a better way to calculate C(X) = Z0(X)^2 rem A(X)

        // let z0z0: Vec<F> = self.z0_s1.iter().flat_map(|&y| [F::zero(), y]).collect();
        let z0: Vec<F> = self.z0_s1.iter().flat_map(|y| [F::zero(), *y]).collect();
        let z0_inv_z1_mod_a = self.redc_z1_impl(&z0, a);
        // println!("z000000: {:?}", z0_inv_z1_mod_a);

        let z1: Vec<F> = self.z1_s0.iter().flat_map(|y| [*y, F::zero()]).collect();
        let z1_inv_z0_mod_a = self.redc_z0_impl(&z1, a);

        println!(
            "{:?}",
            zip(z0_inv_z1_mod_a, z1_inv_z0_mod_a)
                .map(|(a, b)| a + b)
                .collect::<Vec<F>>(),
        );

        // println!("z111111: {:?}", z1_inv_z0_mod_a);
        // println!(
        //     "yo: {:?}",
        //     zip(z0_inv_z1_mod_a, z1_inv_z0_mod_a)
        //         .map(|(a, b)| a * b)
        //         .collect::<Vec<F>>()
        // );

        // println!("tm[: {:?}", tmp);
        // // let tmp: Vec<F> = zip(tmp, &self.z1_s0).map(|(v, z1)| v * z1).collect();
        // println!("tm[: {:?}", tmp);
        todo!()

        // println!("REDC");
        // // let zero = F::zero();
        // // let z0_s_inv: Vec<F> = self
        // //     .z0_s1_inv
        // //     .iter()
        // //     .flat_map(|&y| [zero, y.square().inverse().unwrap()])
        // //     .collect();
        // // let mut c = self.redc_impl(&z0_s_inv, a);
        // // // batch_inversion(&mut c);
        // // println!(
        // //     "g: {:?}",
        // //     self.z0_s1_inv
        // //         .iter()
        // //         .map(|v| v.inverse().unwrap())
        // // //         .collect::<Vec<F>>()
        // // // );
        // // // println!("c: {:?}", c);
        // // // println!("DEGREE: {}", self.degree_impl(&c));
        // // let zero = F::zero();
        // // let z0_s_inv: Vec<F> = self.z0_s1_inv.iter().flat_map(|&y| [zero,
        // // y]).collect(); let mut c = self.redc_impl(&z0_s_inv, a);
        // // // let mut c2 = self.redc_impl(&c, a);
        // // // batch_inversion(&mut c);
        // // println!("z0 original: {:?}", self.z0_s1_inv);
        // // println!("g: {:?}", z0_s_inv);
        // // println!("c: {:?}", c);
        // // // println!("HERE: {:?}", )
        // // // println!("c2: {:?}", c2);
        // // println!("DEGREE: {}", self.degree_impl(&c));
        // let z0z0_s: Vec<F> = self
        //     .z0_s1_inv
        //     .iter()
        //     // .flat_map(|&y| [F::zero(), y.inverse().unwrap().pow([2])])
        //     // .flat_map(|&y| [F::zero(), y.inverse().unwrap()])
        //     .flat_map(|&y| [F::one(), F::one()])
        //     .collect();

        // // let c: Vec<F> = zip(&z0z0_s, a_inv).map(|(z, a)| *z *
        // a).collect(); // println!("C: {:?}", c);
        // // println!("C: {:?}", self.degree(&c));
        // // println!("z0z0_s: {:?}", z0z0_s);
        // // let s = self.f.get_layer(0);
        // let c = self.redc_impl(&z0z0_s, a);
        // let c = self.redc_impl(&c, a);
        // let c = self.redc_impl(&c, a);
        // println!("c: {:?}", c);
        // // let c = self.redc_impl(&c, a);
        // // println!("c: {:?}", c);
        // // let (mut c0, _): (Vec<F>, Vec<F>) = c.chunks(2).map(|c| (c[0],
        // // c[1])).unzip(); batch_inversion(&mut c0);
        // // let c1 = self.extend_impl(&c0, Moiety::S1);
        // // println!("a: {:?}, {:?}", c0, c1);

        // // let h = self.redc_impl(evals, a);
        // // println!("h len {}", h.len());
        // // let hc: Vec<F> = zip(h, c).map(|(h, c)| h * c).collect();
        // // self.redc_impl(&hc, a)
        // println!("Up: {:?}", self.z0_s1_inv);

        // todo!()
    }

    pub fn modular_reduce(&self, evals: &[F], a: &[F]) -> Vec<F> {
        let tree = self.subtree_with_size(evals.len());
        tree.modular_reduce_impl(evals, a)
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

    /// Returns an evaluation of the vanishing polynomial Z(x) = ∏ (x - S_i)
    /// Runtime O(n log^2 n)
    /// Section 7.1 https://arxiv.org/pdf/2107.08473.pdf
    pub fn vanish(&self, vanish_domain: &[F]) -> Vec<F> {
        let tree = self.subtree_with_size(vanish_domain.len() * 2);
        tree.vanish_impl(vanish_domain)
    }

    fn from_tree(f: BinaryTree<F>, isogenies: Vec<Isogeny<F>>) -> Self {
        let subtree = Self::derive_subtree(&f, &isogenies).map(Box::new);
        let f_layers = f.get_layers();
        let n = f.leaves().len();
        let nn = n as u64 / 2;
        let s = f_layers[0];

        // Precompute eval table <X^(n/2) ≀ S>
        let xnn_s = f_layers[0].iter().map(|x| x.pow([nn])).collect::<Vec<F>>();
        let mut xnn_s_inv = xnn_s.clone();
        batch_inversion(&mut xnn_s_inv);

        // Split S into its two moieties S0 and S1
        let (s0, s1): (Vec<F>, Vec<F>) = s.chunks_exact(2).map(|s| (s[0], s[1])).unzip();

        // Generate polynomial decomposition matrices
        // Lemma 3.2 (M_t) https://arxiv.org/abs/2107.08473
        // TODO: change notation
        let mut recomp_matrices = BinaryTree::from(vec![Mat2x2::identity(); n]);
        let mut decomp_matrices = BinaryTree::from(vec![Mat2x2::identity(); n]);
        let recomp_layers = recomp_matrices.get_layers_mut();
        let decomp_layers = decomp_matrices.get_layers_mut();
        for ((recomp_layer, decomp_layer), (l, iso)) in
            zip(zip(recomp_layers, decomp_layers), zip(f_layers, &isogenies))
        {
            let d = l.len() / 2;
            if d == 1 {
                continue;
            }

            let v = &iso.x_denominator_map;
            for (i, (rmat, dmat)) in zip(recomp_layer, decomp_layer).enumerate() {
                let s0 = l[i];
                let s1 = l[i + d];
                let v0 = v.evaluate(&s0).pow([(d / 2 - 1) as u64]);
                let v1 = v.evaluate(&s1).pow([(d / 2 - 1) as u64]);
                *dmat = Mat2x2([[v0, s0 * v0], [v1, s1 * v1]]);
                *rmat = dmat.inverse().unwrap();
            }
        }

        let mut tree = Self {
            f,
            recomp_matrices,
            decomp_matrices,
            isogenies,
            subtree,
            xnn_s,
            xnn_s_inv,
            z0_s1: Vec::new(),
            z1_s0: Vec::new(),
            z1_s0_inv: Vec::new(),
            z0_s1_inv: Vec::new(),
        };

        // Precompute eval tables <Z_0 ≀ S_1> and <Z_1 ≀ S_0> using our partial FFTree
        // Z_0 is the vanishing polynomial of S_0 i.e. Z_0(x) = Π(x - s0_i)
        // TODO: this is a little brittle, might be nice to find a cleaner solution
        match n.cmp(&2) {
            Ordering::Greater => {
                // Compute z0_s1 in O(n log n) using the subtree's vanishing polynomials
                let zero = F::zero();
                let st = tree.subtree.as_ref().unwrap();
                let st_z0_s0: Vec<F> = st.z0_s1.iter().flat_map(|&y| [zero, y]).collect();
                let st_z1_s0: Vec<F> = st.z1_s0.iter().flat_map(|&y| [y, zero]).collect();
                let st_z0_s1 = tree.extend(&st_z0_s0, Moiety::S1);
                let st_z1_s1 = tree.extend(&st_z1_s0, Moiety::S1);
                tree.z0_s1 = zip(st_z0_s1, st_z1_s1).map(|(z0, z1)| z0 * z1).collect();

                // Compute z1_s in O(n log^2 n) - .vanish() uses z0_s1.
                let z1_s = tree.vanish(&s1);
                tree.z1_s0 = z1_s.array_chunks().map(|[z1_s0, _]| *z1_s0).collect();
            }
            Ordering::Equal => {
                tree.z0_s1 = vec![s1[0] - s0[0]];
                tree.z1_s0 = vec![s0[0] - s1[0]];
            }
            Ordering::Less => {}
        }

        tree.z0_s1_inv = tree.z0_s1.clone();
        tree.z1_s0_inv = tree.z1_s0.clone();
        batch_inversion(&mut tree.z0_s1_inv);
        batch_inversion(&mut tree.z1_s0_inv);

        tree
    }

    fn derive_subtree(f: &BinaryTree<F>, isogenies: &[Isogeny<F>]) -> Option<Self> {
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

        Some(Self::from_tree(f_prime, isogenies.to_vec()))
    }

    fn get_layer(&self, i: usize) -> FFTreeLayer<'_, F> {
        FFTreeLayer {
            _l: self.f.get_layer(i),
            _isogeny: &self.isogenies[i],
            decomp_matrices: self.decomp_matrices.get_layer(i),
            recomp_matrices: self.recomp_matrices.get_layer(i),
        }
    }

    pub(crate) fn subtree(&self) -> Option<&Self> {
        Some(self.subtree.as_ref()?)
    }

    /// Returns a FFTree with `n` leaves
    pub(crate) fn subtree_with_size(&self, n: usize) -> &Self {
        assert!(n.is_power_of_two());
        match usize::cmp(&n, &self.f.leaves().len()) {
            Ordering::Less => self.subtree().unwrap().subtree_with_size(n),
            Ordering::Equal => self,
            Ordering::Greater => panic!("FFTree is too small"),
        }
    }
}

#[cfg(test)]
impl<F: PrimeField> FFTree<F> {
    // TODO: maybe convert to method to get subtree of size n?
    // this could be used to simplify public algorithm interfaces on FFTree as well
    pub(crate) fn eval_domain(&self) -> &[F] {
        self.f.leaves()
    }
}

// TODO: implement CanonicalSerialize for Box in arkworks
// TODO: Derive bug "error[E0275]" with recursive field subtree
// TODO: add compress version (with flags perhaps)
impl<F: PrimeField> CanonicalSerialize for FFTree<F> {
    fn serialize_with_mode<W: ark_serialize::Write>(
        &self,
        mut writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        #[deny(unused_variables)]
        let Self {
            f,
            recomp_matrices,
            decomp_matrices,
            isogenies,
            xnn_s,
            xnn_s_inv,
            z0_s1,
            z1_s0,
            z0_s1_inv,
            z1_s0_inv,
            subtree,
        } = self;
        f.serialize_with_mode(&mut writer, compress)?;
        recomp_matrices.serialize_with_mode(&mut writer, compress)?;
        decomp_matrices.serialize_with_mode(&mut writer, compress)?;
        isogenies.serialize_with_mode(&mut writer, compress)?;
        xnn_s.serialize_with_mode(&mut writer, compress)?;
        z0_s1.serialize_with_mode(&mut writer, compress)?;
        z1_s0.serialize_with_mode(&mut writer, compress)?;
        if compress == ark_serialize::Compress::No {
            // Inverses are the cheapest to regenerate (1 mul per element)
            xnn_s_inv.serialize_with_mode(&mut writer, compress)?;
            z0_s1_inv.serialize_with_mode(&mut writer, compress)?;
            z1_s0_inv.serialize_with_mode(&mut writer, compress)?;
        }
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
            recomp_matrices,
            decomp_matrices,
            isogenies,
            xnn_s,
            xnn_s_inv,
            z0_s1,
            z1_s0,
            z0_s1_inv,
            z1_s0_inv,
            subtree,
        } = self;
        let mut size = f.serialized_size(compress)
            + recomp_matrices.serialized_size(compress)
            + decomp_matrices.serialized_size(compress)
            + isogenies.serialized_size(compress)
            + xnn_s.serialized_size(compress)
            + z0_s1.serialized_size(compress)
            + z1_s0.serialized_size(compress)
            // subtree: 1 (for Option state) + subtree size
            + 1 + subtree.as_ref().map_or(0, |v| v.as_ref().serialized_size(compress));
        if compress == Compress::No {
            size += xnn_s_inv.serialized_size(compress)
                + z0_s1_inv.serialized_size(compress)
                + z1_s0_inv.serialized_size(compress);
        }
        size
    }
}

impl<F: PrimeField> Valid for FFTree<F> {
    #[inline]
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        Ok(())
    }
}

// TODO: implement CanonicalDeserialize for Box in arkworks
// TODO: Derive bug "error[E0275]" with recursive field subtree
impl<F: PrimeField> CanonicalDeserialize for FFTree<F> {
    fn deserialize_with_mode<R: ark_serialize::Read>(
        mut reader: R,
        compress: ark_serialize::Compress,
        validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        let f = BinaryTree::deserialize_with_mode(&mut reader, compress, validate)?;
        let recomp_matrices = BinaryTree::deserialize_with_mode(&mut reader, compress, validate)?;
        let decomp_matrices = BinaryTree::deserialize_with_mode(&mut reader, compress, validate)?;
        let isogenies = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
        let xnn_s = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
        let z0_s1 = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
        let z1_s0 = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
        let mut xnn_s_inv: Vec<F>;
        let mut z0_s1_inv: Vec<F>;
        let mut z1_s0_inv: Vec<F>;
        match compress {
            Compress::Yes => {
                xnn_s_inv = xnn_s.clone();
                z0_s1_inv = z0_s1.clone();
                z1_s0_inv = z1_s0.clone();
                batch_inversion(&mut xnn_s_inv);
                batch_inversion(&mut z0_s1_inv);
                batch_inversion(&mut z1_s0_inv);
            }
            Compress::No => {
                xnn_s_inv = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
                z0_s1_inv = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
                z1_s0_inv = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
            }
        }
        let subtree = if bool::deserialize_with_mode(&mut reader, compress, validate)? {
            Some(Box::new(Self::deserialize_with_mode(
                reader, compress, validate,
            )?))
        } else {
            None
        };
        Ok(Self {
            f,
            recomp_matrices,
            decomp_matrices,
            isogenies,
            subtree,
            xnn_s,
            xnn_s_inv,
            z0_s1,
            z1_s0,
            z0_s1_inv,
            z1_s0_inv,
        })
    }
}

/// Returns the two adicity of a point i.e. returns `k` such that 2^k * p = 0.
/// Returns `None` if `p` isn't a point of order 2^k.
pub fn two_adicity<F: PrimeField>(p: Point<F>) -> Option<u32> {
    let mut acc = p;
    for i in 0..2048 {
        if acc.is_zero() {
            return Some(i);
        }
        acc += acc;
    }
    None
}
