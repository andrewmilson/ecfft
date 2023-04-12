
# ECFFT




This library enables fast polynomial arithmetic over any prime field by implementing all the algorithms outlined in [Elliptic Curve Fast Fourier Transform (ECFFT) Part I](https://arxiv.org/pdf/2107.08473.pdf). The implemented algorithms are:

|Algorithm|Description|Runtime|
|:-|:-|:-|
|ENTER|Coefficients to evaluations (fft analogue)|$\mathcal{O}(n\log^2{n})$|
|EXIT|Evaluations to coefficients (ifft analogue)|$\mathcal{O}(n\log^2{n})$|
|DEGREE|Computes a polynomial's degree|$\mathcal{O}(n\log{n})$|
|EXTEND|Extends evaluations from one set to another|$\mathcal{O}(n\log{n})$|
|MEXTEND|EXTEND for special monic polynomials|$\mathcal{O}(n\log{n})$|
|MOD|Calculates the remainder of polynomial division|$\mathcal{O}(n\log{n})$|
|REDC|Computes polynomial analogue of Montgomery's REDC|$\mathcal{O}(n\log{n})$|

## Build FFTrees at compile time

FFTrees are the core datastructure that the ECFFT algorithms are built apon. FFTrees can be generated and serialized at compile time and then be deserialised and used at runtime. This can be preferable since generating FFTrees involves a significant amount of computation. While this approach improves runtime it will significantly blow up a binary's size. TODO: mention about how much runtime and space of generated FFTree is.

```rust
// build.rs

use ecfft::EcFftField;
use ecfft::secp256k1::Fp;

fn main() -> Result<()> {
    let mut f = File::open(concat!(env!("OUT_DIR"), "/fftree.bin"))?;
    let fftree = Fp::build_fftree(1 << 12);
    fftree.serialize_compressed(&mut f)?;
    println!("cargo:rerun-if-changed=build.rs");
    Ok(())
}
```

```rust
// src/main.rs

use ecfft::FFTree;
use ecfft::secp256k1::Fp;
use std::sync::OnceLock;

static FFTREE: OnceLock<FFTree<Fp>> = OnceLock::new();

fn get_fftree() -> &'static FFTree<Fp> {
    FFTREE.get_or_init(|| {
        const PATH = concat!(env!("OUT_DIR"), "/fftree.bin"));
        const BYTES = include_bytes!(path);
        FFTree::deserialize_compressed(BYTES).unwrap()
    })
}

fn main() {
    let fftree = get_fftree();
    let one = Fp::one();
    let zero = Fp::zero();
    // = x^3 + x^2 + 1
    let poly = &[one, zero, one, one];
    let evals = fftree.enter(poly);
    let coeffs = fftree.exit(&evals);
    assert_eq!(poly, coeffs);
}
```

# Performance

TODO