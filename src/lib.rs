#![feature(array_chunks)]

pub mod ec;
pub mod ffforrest;
mod utils;

struct FFTree<F> {
    verticies: Vec<F>,
    // isogenies: Vec<F>,
}

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ec::Curve;
    use crate::ec::Point;
    use crate::ffforrest::FFForrest;
    use ark_ff::Field;
    use ark_ff::One;
    use ark_ff::Zero;
    use ark_ff_optimized::fp31::Fp;
    use num_bigint::BigUint;
    use std::ops::Mul;

    #[test]
    fn it_works() {
        let curve = Curve::new(Fp::one(), Fp::zero());
        let n = 1u32 << 23;
        let cardinality = BigUint::from(1u32 << 31);
        let primitive_point = Point::new(Fp::from(1048755163u32), Fp::from(279503108u32), curve);
        let generator = primitive_point * (cardinality / n);

        let ffforrest = FFForrest::new(primitive_point, generator);
        println!("YO");
    }
}
