use crate::ec::Isogeny;
use crate::ec::Point;
use ark_ff::PrimeField;
use ark_ff::Zero;

pub struct FFForrest<F: PrimeField> {
    coset_offset: Point<F>,
    generator: Point<F>,
    vertices: Vec<F>,
    isogenies: Vec<Isogeny<F>>,
}

impl<F: PrimeField> FFForrest<F> {
    pub fn new(coset_offset: Point<F>, generator: Point<F>) -> Self {
        assert_ne!(coset_offset, generator);
        assert_eq!(coset_offset.curve, generator.curve);
        let two_addicity = two_addicity(generator).expect("invalid generator");
        let n = 1 << two_addicity;

        // generate leaf vertices
        let mut vertices = vec![F::zero(); 2 * n];
        let mut acc = Point::zero();
        for vertex in &mut vertices[n..] {
            *vertex = (coset_offset + acc).x;
            acc += acc;
        }

        // generate remaining vertices and collect isogenies.
        let mut layers = get_tree_layers(&mut vertices);
        let mut prev_layer_domain = generator.curve.unwrap();
        let mut isogenies = Vec::new();
        for [prev_layer, layer] in layers.array_chunks_mut() {
            let layer_size = layer.len();
            let isogeny = prev_layer_domain
                .two_isogenies()
                .into_iter()
                .find_map(|isogeny| {
                    for (i, dst) in layer.iter_mut().enumerate() {
                        let prev_lhs = prev_layer[i];
                        let prev_rhs = prev_layer[i + layer_size];
                        let prev_lhs_mapped = isogeny.map_x(&prev_lhs)?;
                        let prev_rhs_mapped = isogeny.map_x(&prev_rhs)?;
                        if prev_lhs_mapped != prev_rhs_mapped {
                            return None;
                        }
                        *dst = prev_lhs_mapped;
                    }
                    Some(isogeny)
                })
                .expect("cannot find a suitable isogeny");
            prev_layer_domain = isogeny.codomain;
            isogenies.push(isogeny);
        }

        Self {
            coset_offset,
            generator,
            vertices,
            isogenies,
        }
    }
}

/// Returns all layers of a binary tree.
/// ```text
///          a         <- layers[d-1] = [a]
///         / \
///        b   c       <- layers[d-2] = [b, c]
///       / \ / \
///    ...  ...  ...
///    / \       / \
///   w   x ... y   z  <- layers[0] = [w, x, ..., y, z]
/// ```
fn get_tree_layers<T>(vertecies: &mut [T]) -> Vec<&mut [T]> {
    let n = vertecies.len();
    assert!(n.is_power_of_two());
    let depth = n.ilog2();
    let mut res = Vec::new();
    (0..depth).rev().fold(vertecies, |rem, i| {
        let (lhs, rhs) = rem.split_at_mut(1 << i);
        res.push(rhs);
        lhs
    });
    res
}

/// Returns the two addicity of a point i.e. returns `n` such that 2^n * p = 0.
/// Returns `None` if `p` isn't a point of order 2^n or if `n` exceeds 31.
fn two_addicity<F: PrimeField>(p: Point<F>) -> Option<u32> {
    if p.is_zero() {
        return None;
    }
    let mut acc = p;
    for i in 0..=31 {
        if acc.is_zero() {
            return Some(i);
        }
        acc += acc;
    }
    None
}
