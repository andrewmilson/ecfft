# ECFFT

This library enables fast polynomial arithmetic over any prime field by implementing all the algorithms outlined in [Elliptic Curve Fast Fourier Transform (ECFFT) Part I](https://arxiv.org/pdf/2107.08473.pdf). The implemented algorithms are:
* **_ENTER_** - coefficients to evaluations (fft analogue) in $\mathcal{O}(n\log^2{n})$
* **_EXIT_** - evaluations to coefficients (ifft analogue) in $\mathcal{O}(n\log^2{n})$
* **_DEGREE_** - computes a polynomial's degree in $\mathcal{O}(n\log{n})$
* **_EXTEND_** - extends evaluations from one set to another in $\mathcal{O}(n\log{n})$
* **_MEXTEND_** - _EXTEND_  for special monic polynomials in $\mathcal{O}(n\log{n})$
* **_MOD_** - calculates the remainder of polynomial division in $\mathcal{O}(n\log{n})$
* **_REDC_** - computes polynomial analogue of Montgomery's REDC in $\mathcal{O}(n\log{n})$

## Build FFTrees at compile time

FFTrees can be generated and serialized at compile time and can then be deserialised and used at runtime. This can be preferable since generating FFTrees involves a significant amount of computation. This approach improves runtime but can significantly blow up the binary in size.

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

TODO: Given the runtime of generating an FFTree can be O(n^2) the library provides the ability to specify what algorithms will be required at runtime. This can dramatically speed up the time it takes to generate an FFTree since only certain precomputations can be omitted since they are only used for specific algorithms.

TODO: Another issue is the sheer size of the FFTree. The precomputation FFTree space complexity of O(n) but this can quickly get to many MB/GB if the


# Performance
