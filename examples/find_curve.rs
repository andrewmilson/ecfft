#![feature(let_chains)]

use ark_ff::Field;
use ark_ff::Fp256;
use ark_ff::MontBackend;
use ark_ff::MontConfig;
use ecfft::ec::GoodCurve;
use ecfft::ec::Point;
use ecfft::find_curve::find_curve;
use std::sync::atomic::AtomicU32;
use std::sync::atomic::Ordering;

#[derive(MontConfig)]
#[modulus = "3618502788666131213697322783095070105623107215331596699973092056135872020481"]
#[generator = "3"]
pub struct FpMontConfig;

/// The 252-bit prime field used by StarkWare for Cairo
/// Field has modulus `2^251 + 17 * 2^192 + 1`
pub type Fp = Fp256<MontBackend<FpMontConfig, 4>>;

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
                                let bb = b.square();
                                let Point { x, y, .. } = g;
                                println!("n={n}, a = {a}, bb = {bb}, p=({x}, {y})",)
                            }
                            GoodCurve::Even { b: _ } => todo!(),
                        }
                    }
                }
            });
        }
    });
}
