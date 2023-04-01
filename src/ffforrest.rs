use crate::ec::Isogeny;
use crate::ec::Point;
use crate::utils::get_layers;
use crate::utils::get_layers_mut;
use ark_ff::PrimeField;
use ark_ff::Zero;
use ark_poly::Polynomial;
use std::iter::zip;
use std::ops::Mul;

pub struct FFTreeLayer<'a, F: PrimeField> {
    pub l: &'a [F],
    pub decomp_matrices: &'a [Mat2x2<F>],
    pub recomp_matrices: &'a [Mat2x2<F>],
    pub isogeny: &'a Isogeny<F>,
}

pub struct FFForrest<F: PrimeField> {
    coset_offset: Point<F>,
    generator: Point<F>,
    pub vertices: Vec<F>,
    isogenies: Vec<Isogeny<F>>,
    decomp_matrices: Vec<Mat2x2<F>>,
    recomp_matrices: Vec<Mat2x2<F>>,
}

impl<F: PrimeField> FFForrest<F> {
    pub fn new(coset_offset: Point<F>, generator: Point<F>) -> Self {
        assert_ne!(coset_offset, generator);
        assert_eq!(coset_offset.curve, generator.curve);
        let log_n = two_adicity(generator).expect("generator must have order 2^k");
        assert!(log_n < 32, "generator has unsupported two adicity");
        let n = 1 << log_n;

        // generate leaf vertices
        let mut vertices = vec![F::zero(); 2 * n];
        let mut acc = Point::zero();
        for vertex in &mut vertices[n..] {
            *vertex = (coset_offset + acc).x;
            acc += generator;
        }

        // generate remaining vertices and collect isogenies.
        let mut vertex_layers = get_layers_mut(&mut vertices);
        let mut prev_layer_domain = generator.curve.unwrap();
        let mut isogenies = Vec::new();
        let mut g = generator;
        // TODO: use array_windows_mut
        for i in 1..vertex_layers.len() {
            println!("KKLLKLKLKLKLKLKLK{i}");
            let (prev_layers, layers) = vertex_layers.split_at_mut(i);
            let prev_layer = prev_layers.last_mut().unwrap();
            let layer = layers.first_mut().unwrap();
            let layer_size = layer.len();
            let isogeny = prev_layer_domain
                .two_isogenies()
                .into_iter()
                .find_map(|isogeny| {
                    // let new_two = two_adicity(isogeny.map(&g));
                    // println!("NT: {:?} {}", new_two, layer.len().ilog2());
                    // if new_two.unwrap() == layer.len().ilog2() {
                    //     println!("TWO ADDICITY: {:?}", new_two);
                    //     g = isogeny.map(&g);
                    // }

                    for (i, dst) in layer.iter_mut().enumerate() {
                        let prev_lhs = prev_layer[i];
                        let prev_rhs = prev_layer[i + layer_size];
                        let lhs_mapped = isogeny.map_x(&prev_lhs)?;
                        let rhs_mapped = isogeny.map_x(&prev_rhs)?;
                        if lhs_mapped != rhs_mapped {
                            return None;
                        }
                        *dst = lhs_mapped;
                    }

                    Some(isogeny)
                })
                .expect("cannot find a suitable isogeny");
            prev_layer_domain = isogeny.codomain;
            isogenies.push(isogeny);
        }

        // Generate polynomial decomposition matrices
        // Lemma 3.2 (M_t) https://arxiv.org/abs/2107.08473
        let mut decomp_matrices = vec![Mat2x2::identity(); n / 2];
        let mut recomp_matrices = vec![Mat2x2::identity(); n / 2];
        let decomp_layers = get_layers_mut(&mut decomp_matrices);
        let recomp_layers = get_layers_mut(&mut recomp_matrices);
        for ((decomp_layer, recomp_layer), (xs, iso)) in zip(
            zip(decomp_layers, recomp_layers),
            zip(vertex_layers, &isogenies),
        ) {
            let d = xs.len();
            let v = &iso.x_denominator_map;
            for (i, (dmat, rmat)) in zip(decomp_layer, recomp_layer).enumerate() {
                let s0 = xs[i * 2];
                let s1 = xs[i * 2 + d / 2];
                let v0 = v.evaluate(&s0).pow([(d / 4 - 1) as u64]);
                let v1 = v.evaluate(&s1).pow([(d / 4 - 1) as u64]);
                *rmat = Mat2x2([[v0, s0 * v0], [v1, s1 * v1]]).inverse().unwrap();

                let s0_prime = xs[i * 2 + 1];
                let s1_prime = xs[i * 2 + 1 + d / 2];
                let v0_prime = v.evaluate(&s0_prime).pow([(d / 4 - 1) as u64]);
                let v1_prime = v.evaluate(&s1_prime).pow([(d / 4 - 1) as u64]);
                *dmat = Mat2x2([
                    [v0_prime, s0_prime * v0_prime],
                    [v1_prime, s1_prime * v1_prime],
                ]);
                println!("LL: {:?}", dmat);
            }
        }

        println!("{:?}", recomp_matrices);

        Self {
            coset_offset,
            generator,
            vertices,
            isogenies,
            decomp_matrices,
            recomp_matrices,
        }
    }

    fn num_layers(&self) -> usize {
        self.isogenies.len()
    }

    pub fn get_layer(&self, i: usize) -> FFTreeLayer<'_, F> {
        let isogeny = &self.isogenies[i];
        let decomp_matrices = get_layers(&self.decomp_matrices)[i];
        let recomp_matrices = get_layers(&self.recomp_matrices)[i];
        let l = get_layers(&self.vertices)[i];
        FFTreeLayer {
            l,
            decomp_matrices,
            recomp_matrices,
            isogeny,
        }
    }

    // EXTEND
    pub fn extend(&self, evals: &[F]) -> Vec<F> {
        if evals.len() == 1 {
            evals.to_vec()
        } else {
            let n = evals.len();
            let layer_id = self.num_layers() - 1 - n.ilog2() as usize;
            let layer = self.get_layer(layer_id);

            // π_0 and π_1
            let (evals0, evals1): (Vec<F>, Vec<F>) = layer
                .recomp_matrices
                .iter()
                .enumerate()
                .map(|(i, m)| {
                    let res = m * &[evals[i], evals[i + n / 2]];
                    (res[0], res[1])
                })
                .unzip();

            let evals0_prime = self.extend(&evals0);
            let evals1_prime = self.extend(&evals1);

            // π_0' and π_1'
            let (lhs, rhs): (Vec<F>, Vec<F>) = layer
                .decomp_matrices
                .iter()
                .enumerate()
                .map(|(i, m)| {
                    let res = m * &[evals0_prime[i], evals1_prime[i]];
                    (res[0], res[1])
                })
                .unzip();

            [lhs, rhs].concat()
        }
    }

    // ENTER
    pub fn enter(&self, coeffs: &[F]) -> Vec<F> {
        if coeffs.len() == 1 {
            coeffs.to_vec()
        } else {
            let n = coeffs.len();

            let u0 = self.enter(&coeffs[0..n / 2]);
            let u1 = self.extend(&u0);
            let v0 = self.enter(&coeffs[n / 2..]);
            let v1 = self.extend(&v0);

            let layer_id = self.num_layers() - 1 - n.ilog2() as usize;
            let l = get_layers(&self.vertices)[layer_id];

            l.array_chunks()
                .enumerate()
                // .take(n / 2)
                .flat_map(|(i, [s0, s1])| {
                    [
                        u0[i] + s0.pow([(n / 2) as u64]) * v0[i],
                        u1[i] + s1.pow([(n / 2) as u64]) * v1[i],
                    ]
                })
                .collect()

            // [lhs, rhs].concat()
        }
    }
}

/// Returns the two adicity of a point i.e. returns `n` such that 2^n * p = 0.
/// Returns `None` if `p` isn't a point of order 2^n.
fn two_adicity<F: PrimeField>(p: Point<F>) -> Option<u32> {
    if p.is_zero() {
        return None;
    }
    let mut acc = p;
    for i in 0..4096 {
        if acc.is_zero() {
            return Some(i);
        }
        acc += acc;
    }
    None
}

#[derive(Clone, Copy, Default, Debug)]
pub struct Mat2x2<F>(pub [[F; 2]; 2]);

impl<F: PrimeField> Mat2x2<F> {
    pub fn identity() -> Mat2x2<F> {
        Self([[F::one(), F::zero()], [F::zero(), F::one()]])
    }

    pub fn inverse(&self) -> Option<Self> {
        let det_inv = self.determinant().inverse();
        if det_inv.is_none() {
            println!("{:?}", self);
            return None;
        }
        let det_inv = det_inv.unwrap();
        Some(Self([
            [self.0[1][1] * det_inv, -self.0[0][1] * det_inv],
            [-self.0[1][0] * det_inv, self.0[0][0] * det_inv],
        ]))
    }

    pub fn determinant(&self) -> F {
        self.0[0][0] * self.0[1][1] - self.0[0][1] * self.0[1][0]
    }
}

impl<F: PrimeField> Mul<&[F; 2]> for &Mat2x2<F> {
    type Output = [F; 2];

    fn mul(self, rhs: &[F; 2]) -> Self::Output {
        [
            self.0[0][0] * rhs[0] + self.0[0][1] * rhs[1],
            self.0[1][0] * rhs[0] + self.0[1][1] * rhs[1],
        ]
    }
}
