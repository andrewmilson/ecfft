> **Warning**
> this code is a work in progress 👷‍♂️🚧. All the algorithms work so please try it out and provide feedback. If it addresses a use-case that is important to you please open an issue to discuss it or get in touch [andrew.j.milson@gmail.com](mailto:andrew.j.milson@gmail.com)

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
|VANISH|Generates a vanishing polynomial ([from section 7.1](https://arxiv.org/pdf/2107.08473.pdf))|$\mathcal{O}(n\log^2{n})$|

There are also some relevant algorithms implemented from [ECFFT Part II](https://www.math.toronto.edu/swastik/ECFFT2.pdf). The algorithms are:

|Algorithm|Description|Runtime|
|:-|:-|:-|
|FIND_CURVE|Finds a curve over $\mathbb{F}_q$ with a cyclic subgroup of order $2^k$ |$\mathcal{O}(2^k\log{q})$|



## Build FFTrees at compile time

FFTrees are the core datastructure that the ECFFT algorithms are built apon. FFTrees can be generated and serialized at compile time and then be deserialised and used at runtime. This can be preferable since generating FFTrees involves a significant amount of computation. While this approach improves runtime it will significantly blow up a binary's size. Generating a FFTree capable of evaluating/interpolating degree $n$ polynomials takes $\mathcal{O}(n\log^3{n})$ - the space complexity of this FFTree is $\mathcal{O}(n)$.

```rust
// build.rs

use ark_serialize::CanonicalSerialize;
use ecfft::{secp256k1::Fp, FftreeField};
use std::{env, fs::File, io, path::Path};

fn main() -> io::Result<()> {
    let fftree = Fp::build_fftree(1 << 16).unwrap();
    let out_dir = env::var_os("OUT_DIR").unwrap();
    let path = Path::new(&out_dir).join("fftree");
    fftree.serialize_compressed(File::create(path)?).unwrap();
    println!("cargo:rerun-if-changed=build.rs");
    Ok(())
}
```

```rust
// src/main.rs

use ark_ff::One;
use ark_serialize::CanonicalDeserialize;
use ecfft::{ecfft::FFTree, secp256k1::Fp};
use std::sync::OnceLock;

static FFTREE: OnceLock<FFTree<Fp>> = OnceLock::new();

fn get_fftree() -> &'static FFTree<Fp> {
    const BYTES: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/fftree"));
    FFTREE.get_or_init(|| FFTree::deserialize_compressed(BYTES).unwrap())
}

fn main() {
    let fftree = get_fftree();
    // = x^65535 + x^65534 + ... + x + 1
    let poly = vec![Fp::one(); 1 << 16];
    let evals = fftree.enter(&poly);
    let coeffs = fftree.exit(&evals);
    assert_eq!(poly, coeffs);
}
```

# Performance

TODO