use ecfft::m31::Fp;
use ecfft::EcFftField;

fn main() {
    let fftree = Fp::build_fftree(1 << 5).unwrap();
    let _ = fftree.subtree_with_size(8);
}
