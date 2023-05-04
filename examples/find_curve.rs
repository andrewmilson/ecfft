#![feature(let_chains)]

use ark_ff::Field;
use ecfft::ec::GoodCurve;
use ecfft::find_curve::find_curve;
// use ecfft::m31::Fp;
use ecfft::secp256k1::Fp;
use std::sync::atomic::AtomicU32;
use std::sync::atomic::Ordering;

fn main() {
    let max_n = AtomicU32::new(10);

    rayon::scope(|s| {
        for _ in 0..10 {
            s.spawn(|_| {
                let mut rng = rand::thread_rng();
                loop {
                    let (n, g) = find_curve::<Fp>(&mut rng, max_n.load(Ordering::Relaxed) + 1);
                    if n > max_n.load(Ordering::Relaxed) {
                        max_n.store(n, Ordering::Relaxed);
                        let curve = g.curve.unwrap();
                        match curve {
                            GoodCurve::Odd { a, b } => {
                                println!("n={n}, a={a}, bb={}, p=({}, {})", b.square(), g.x, g.y)
                            }
                            GoodCurve::Even { b: _ } => todo!(),
                        }
                    }
                }
            });
        }
    });
}
