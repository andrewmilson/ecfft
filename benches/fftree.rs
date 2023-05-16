use criterion::criterion_group;
use criterion::criterion_main;
use criterion::BenchmarkId;
use criterion::Criterion;
use ecfft::m31;
use ecfft::secp256k1;
use ecfft::FftreeField;
use rand::rngs::StdRng;
use rand::SeedableRng;

const BENCHMARK_INPUT_SIZES: [usize; 1] = [2048];

const FIELD_DESCRIPTION_M31: &str = "31 bit Mersenne prime field";
const FIELD_DESCRIPTION_SECP256K1: &str = "secp256k1's prime field";

fn bench_ecfft_algorithms<F: FftreeField>(c: &mut Criterion, field_description: &str) {
    let mut rng = StdRng::from_seed([1; 32]);
    let mut group = c.benchmark_group(format!("ECFFT algorithms ({field_description})"));
    group.sample_size(10);

    for n in BENCHMARK_INPUT_SIZES {
        let vals: Vec<F> = (0..n).map(|_| F::rand(&mut rng)).collect();
        let fftree = F::build_fftree(n * 2).unwrap();

        group.bench_with_input(BenchmarkId::new("ENTER", n), &n, |b, _| {
            b.iter(|| fftree.enter(&vals))
        });

        group.bench_with_input(BenchmarkId::new("EXIT", n), &n, |b, _| {
            b.iter(|| fftree.exit(&vals))
        });

        group.bench_with_input(BenchmarkId::new("DEGREE", n), &n, |b, _| {
            b.iter(|| fftree.degree(&vals))
        });

        group.bench_with_input(BenchmarkId::new("EXTEND", n), &n, |b, _| {
            b.iter(|| fftree.extend(&vals, ecfft::Moiety::S1))
        });

        group.bench_with_input(BenchmarkId::new("MEXTEND", n), &n, |b, _| {
            b.iter(|| fftree.mextend(&vals, ecfft::Moiety::S1))
        });

        group.bench_with_input(BenchmarkId::new("MOD", n), &n, |b, _| {
            b.iter(|| fftree.modular_reduce(&vals, &fftree.xnn_s, &fftree.z0z0_rem_xnn_s))
        });

        group.bench_with_input(BenchmarkId::new("REDC", n), &n, |b, _| {
            b.iter(|| fftree.redc_z0(&vals, &fftree.xnn_s))
        });

        group.bench_with_input(BenchmarkId::new("VANISH", n), &n, |b, _| {
            b.iter(|| fftree.vanish(&vals))
        });
    }

    group.finish();
}

fn bench_fftree<F: FftreeField>(c: &mut Criterion, field_description: &str) {
    let mut group = c.benchmark_group(format!("FFTree ({field_description})"));
    group.sample_size(10);

    for n in BENCHMARK_INPUT_SIZES {
        group.bench_with_input(BenchmarkId::new("generate", n), &n, |b, _| {
            b.iter(|| F::build_fftree(n).unwrap())
        });
    }

    group.finish();
}

fn ecfft_algorithm_benches(c: &mut Criterion) {
    bench_ecfft_algorithms::<m31::Fp>(c, FIELD_DESCRIPTION_M31);
    bench_ecfft_algorithms::<secp256k1::Fp>(c, FIELD_DESCRIPTION_SECP256K1);
}

fn fftree_benches(c: &mut Criterion) {
    bench_fftree::<m31::Fp>(c, FIELD_DESCRIPTION_M31);
    bench_fftree::<secp256k1::Fp>(c, FIELD_DESCRIPTION_SECP256K1);
}

criterion_group!(fftree_group, ecfft_algorithm_benches, fftree_benches);
criterion_main!(fftree_group);
