use ark_ff::Fp256;
use ark_ff::MontBackend;
use ark_ff::MontConfig;
use ark_ff::UniformRand;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use criterion::criterion_group;
use criterion::criterion_main;
use criterion::BenchmarkId;
use criterion::Criterion;
use ecfft::secp256k1::Fp as EcfftField;
use ecfft::FftreeField;
use rand::rngs::StdRng;
use rand::SeedableRng;

const BENCHMARK_INPUT_SIZES: [usize; 1] = [8192];

#[derive(MontConfig)]
#[modulus = "3618502788666131213697322783095070105623107215331596699973092056135872020481"]
#[generator = "3"]
pub struct FpMontConfig;

pub type FftField = Fp256<MontBackend<FpMontConfig, 4>>;

// TODO: naive evaluation/interpolation
fn ecfft_vs_fft_benches(c: &mut Criterion) {
    let mut rng = StdRng::from_seed([1; 32]);
    let mut group = c.benchmark_group("ECFFT vs FFT");
    group.sample_size(10);

    for n in BENCHMARK_INPUT_SIZES {
        let ecfft_vals: Vec<EcfftField> = (0..n).map(|_| EcfftField::rand(&mut rng)).collect();
        let fft_vals: Vec<FftField> = (0..n).map(|_| FftField::rand(&mut rng)).collect();
        let domain = Radix2EvaluationDomain::<FftField>::new(n).unwrap();
        let fftree = EcfftField::build_fftree(n * 2).unwrap();

        group.bench_with_input(BenchmarkId::new("evaluate/ECFFT", n), &n, |b, _| {
            b.iter(|| fftree.enter(&ecfft_vals))
        });

        group.bench_with_input(BenchmarkId::new("interpolate/ECFFT", n), &n, |b, _| {
            b.iter(|| fftree.exit(&ecfft_vals))
        });

        group.bench_with_input(BenchmarkId::new("evaluate/FFT", n), &n, |b, _| {
            b.iter(|| domain.fft(&fft_vals))
        });

        group.bench_with_input(BenchmarkId::new("interpolate/FFT", n), &n, |b, _| {
            b.iter(|| domain.ifft(&fft_vals))
        });
    }

    group.finish();
}

criterion_group!(benches, ecfft_vs_fft_benches);
criterion_main!(benches);
