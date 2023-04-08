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
use ark_serialize::Valid;
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
// pub struct FFTree<F: PrimeField> {
//     pub f: BinaryTree<F>,
//     pub recomp: BinaryTree<Mat2x2<F>>,
//     pub decomp: BinaryTree<Mat2x2<F>>,
//     pub s0: Option<Box<Self>>,
//     pub s1: Option<Box<Self>>,
//     pub z: Option<Vec<F>>,
// }

#[derive(Clone, Debug)]
pub struct FFTree<F: PrimeField> {
    pub f: BinaryTree<F>,
    recomp_matrices: BinaryTree<Mat2x2<F>>,
    decomp_matrices: BinaryTree<Mat2x2<F>>,
    isogenies: Vec<Isogeny<F>>,
    subtree: Option<Box<Self>>,
    xnn_s: Vec<F>,     // = <X^(n/2) ≀ S>
    z0_s1_inv: Vec<F>, // = <1/Z_0 ≀ S_1>
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

    pub fn extend(&self, evals: &[F]) -> Vec<F> {
        self.extend_moiety(evals, Moiety::S1)
    }

    /// Extends evals on the chosen moiety
    pub fn extend_moiety(&self, evals: &[F], moiety: Moiety) -> Vec<F> {
        assert!(evals.len().is_power_of_two());
        match usize::cmp(&(evals.len() * 2), &self.f.leaves().len()) {
            std::cmp::Ordering::Less => self.subtree().unwrap().extend(evals),
            std::cmp::Ordering::Equal => self.extend_impl(evals, moiety),
            std::cmp::Ordering::Greater => panic!("FFTree is too small"),
        }
    }

    pub fn enter_impl(&self, coeffs: &[F]) -> Vec<F> {
        let n = coeffs.len();
        if n == 1 {
            return coeffs.to_vec();
        }

        let subtree = self.subtree().unwrap();
        let u0 = subtree.enter(&coeffs[0..n / 2]);
        let v0 = subtree.enter(&coeffs[n / 2..]);
        let u1 = self.extend(&u0);
        let v1 = self.extend(&v0);

        self.xnn_s
            .array_chunks()
            .enumerate()
            .flat_map(|(i, [s0, s1])| [u0[i] + v0[i] * s0, u1[i] + v1[i] * s1])
            .collect()
    }

    pub fn enter(&self, coeffs: &[F]) -> Vec<F> {
        assert!(coeffs.len().is_power_of_two());
        match usize::cmp(&coeffs.len(), &self.f.leaves().len()) {
            std::cmp::Ordering::Less => self.subtree().unwrap().enter(coeffs),
            std::cmp::Ordering::Equal => self.enter_impl(coeffs),
            std::cmp::Ordering::Greater => panic!("FFTree is too small"),
        }
    }

    pub fn degree_impl(&self, evals: &[F]) -> isize {
        let n = evals.len() as isize;
        if n == 1 {
            return 0;
        }

        let subtree = self.subtree().unwrap();
        let (e0, e1): (Vec<F>, Vec<F>) = evals.chunks(2).map(|e| (e[0], e[1])).unzip();

        // The intuition is that if degree < n / 2 then all coefficients on the right
        // side will be 0. Therefore if degree < n / 2 then e1 == extend(e0)
        let g1 = self.extend_impl(&e0, Moiety::S1);
        if g1 == e1 {
            return subtree.degree_impl(&e0);
        }

        // compute <(π-g)/Z_0 ≀ S1>
        let t1: Vec<F> = zip(zip(e1, g1), &self.z0_s1_inv)
            .map(|((e1, g1), z0_inv)| (e1 - g1) * z0_inv)
            .collect();
        let t0 = self.extend_impl(&t1, Moiety::S0);
        n / 2 + subtree.degree_impl(&t0)
    }

    pub fn degree(&self, evals: &[F]) -> usize {
        assert!(evals.len().is_power_of_two());
        match usize::cmp(&evals.len(), &self.f.leaves().len()) {
            std::cmp::Ordering::Less => self.subtree().unwrap().degree(evals),
            std::cmp::Ordering::Equal => self.degree_impl(evals).max(0).unsigned_abs(),
            std::cmp::Ordering::Greater => panic!("FFTree is too small"),
        }
    }

    pub fn exit_impl(&self, evals: &[F]) -> Vec<F> {
        let n = evals.len();
        if n == 1 {
            return evals.to_vec();
        }

        todo!()
    }

    pub fn exit(&self, evals: &[F]) -> Vec<F> {
        assert!(evals.len().is_power_of_two());
        match usize::cmp(&evals.len(), &self.f.leaves().len()) {
            std::cmp::Ordering::Less => self.subtree().unwrap().exit(evals),
            std::cmp::Ordering::Equal => self.exit_impl(evals),
            std::cmp::Ordering::Greater => panic!("FFTree is too small"),
        }
    }

    fn from_tree(f: BinaryTree<F>, isogenies: Vec<Isogeny<F>>) -> Self {
        let f_layers = f.get_layers();
        let n = f.leaves().len();
        let nn = n as u64 / 2;
        let s = f_layers[0];

        // Precompute eval table <X^(n/2) ≀ S>
        let xnn_s = f_layers[0].iter().map(|x| x.pow([nn])).collect();

        let (s0, s1): (Vec<F>, Vec<F>) = s.chunks_exact(2).map(|s| (s[0], s[1])).unzip();

        // Precompute eval table <Z_0 ≀ S_1>
        // i.e. evaluate S_1 over Z_0.
        // Z_0 is the vanishing polynomial of S_0 i.e. Z_0(x) = Π(x - s0_i)
        let mut z0_s1_inv: Vec<F> = s1
            .iter()
            .map(|&s1| s0.iter().fold(F::one(), |acc, &s0| acc * (s1 - s0)))
            .collect();
        batch_inversion(&mut z0_s1_inv);

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
            subtree: None,
            xnn_s,
            z0_s1_inv,
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
            for (s_prime, s) in l_prime.iter_mut().zip(l.iter().step_by(2)) {
                *s_prime = *s;
            }
        }

        Some(Self::from_tree(f_prime, isogenies.clone()))
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
}

// TODO: implement CanonicalSerialize for Box in arkworks
// TODO: Derive bug "error[E0275]" with recursive field subtree
impl<F: PrimeField> CanonicalSerialize for FFTree<F> {
    fn serialize_with_mode<W: ark_serialize::Write>(
        &self,
        mut writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        self.f.serialize_with_mode(&mut writer, compress)?;
        self.recomp_matrices
            .serialize_with_mode(&mut writer, compress)?;
        self.decomp_matrices
            .serialize_with_mode(&mut writer, compress)?;
        self.isogenies.serialize_with_mode(&mut writer, compress)?;
        self.xnn_s.serialize_with_mode(&mut writer, compress)?;
        self.z0_s1_inv.serialize_with_mode(&mut writer, compress)?;
        // TODO: get "error[E0275]: overflow evaluating the requirement" for:
        // self.subtree.as_ref().map(Box::as_ref).serialize_with_mode(...)
        self.subtree
            .is_some()
            .serialize_with_mode(&mut writer, compress)?;
        if let Some(subtree) = self.subtree.as_ref().map(Box::as_ref) {
            subtree.serialize_with_mode(writer, compress)?;
        }
        Ok(())
    }

    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        self.f.serialized_size(compress)
            + self.recomp_matrices.serialized_size(compress)
            + self.recomp_matrices.serialized_size(compress)
            // self.subtree: 1 (for Option state) + subtree size
            + 1
            + self
                .subtree
                .as_ref()
                .map_or(0, |v| v.as_ref().serialized_size(compress))
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
        let z0_s1_inv = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
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
            z0_s1_inv,
        })
    }
}

#[cfg(test)]
impl<F: PrimeField> FFTree<F> {
    pub(crate) fn eval_domain(&self, n: usize) -> &[F] {
        assert!(n.is_power_of_two());
        let eval_domain = self.f.leaves();
        match usize::cmp(&n, &eval_domain.len()) {
            std::cmp::Ordering::Less => self.subtree().unwrap().eval_domain(n),
            std::cmp::Ordering::Equal => eval_domain,
            std::cmp::Ordering::Greater => panic!("FFTree is too small"),
        }
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
