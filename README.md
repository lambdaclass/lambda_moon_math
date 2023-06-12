# LambdaWorks

<div>

[![Telegram Chat][tg-badge]][tg-url]
[![codecov](https://img.shields.io/codecov/c/github/lambdaclass/lambdaworks)](https://codecov.io/gh/lambdaclass/lambdaworks)

[tg-badge]: https://img.shields.io/static/v1?color=green&logo=telegram&label=chat&style=flat&message=join
[tg-url]: https://t.me/+98Whlzql7Hs0MDZh

</div>

From the heights of these towers of fields, forty centuries of mathematics look down on us. The library for kids who wanna learn how to do SNARKs and learn other cryptographic stuff too.

![FqJmqwtWABozHi_](https://user-images.githubusercontent.com/569014/222294400-4a9e2a46-2823-4a98-b0f4-2293d372cccb.jpeg)

If you use ```Lambdaworks``` libraries in your research projects, please cite them using the following template:

```
@software{Lambdaworks,
  author={Lambdaworks contributors},
  title={Lambdaworks},
  url={https://github.com/lambdaclass/lambdaworks},
  year={2023}
}
```
## Main features

- [Finite Field Algebra](https://github.com/lambdaclass/lambdaworks/tree/main/math/src/field)
- [Polynomial operations](https://github.com/lambdaclass/lambdaworks/blob/main/math/src/polynomial.rs)
- [Fast Fourier Transform](https://github.com/lambdaclass/lambdaworks/tree/main/fft)
- [Elliptic curves](https://github.com/lambdaclass/lambdaworks/tree/main/math/src/elliptic_curve)
- [Multiscalar multiplication](https://github.com/lambdaclass/lambdaworks/tree/main/math/src/msm)
- [ABI compatible KZG commitment scheme - EIP4844](https://github.com/lambdaclass/lambdaworks_kzg)
- [CAIRO prover](https://github.com/lambdaclass/lambdaworks_cairo_prover/tree/main)
- [PLONK prover](https://github.com/lambdaclass/lambdaworks_plonk_prover)

## Building blocks

- Finite Field Algebra
- Elliptic curve models
- Elliptic curve operations
- Arithmetization schemes
- Polynomial commitment schemes
- PIOP
- Cryptographic tools
- Advanced tools: aggregation, recursion, accumulation
- Protocols
- Gadgets

## Blocks
### Finite Field Algebra
- Big integer representation
- Basic algebra: addition, multiplication, subtraction, inversion, square root (Tonelli–Shanks)
- Field extensions
- Number theoretic transform
- Polynomial operations
- Fast Fourier Transform
- Montgomery and Barrett

### Elliptic curve models
- BLS12-381 (H)
- BLS12-377 (H)
- secp256k1 (H)
- Ed25519 (H)
- Jubjub (M)
- BN254 (M)
- Pasta: Pallas and Vesta (L)
- Forms:
1. Affine (H)
2. Projective (H)
3. Montgomery (M)
4. Twisted Edwards (H)
5. Jacobi (L)

### Elliptic curve operations
- Add, double, scalar multiplication.
- Multiscalar multiplication (Pippenger)
- Weyl, Tate and Ate pairings.

### Arithmetization
- R1CS - gadgets (H)
- AIR (M)
- Plonkish (H)
- ACIR (L)

### Polynomial commitment schemes
- KZG and variants
- Hashing
- Inner product arguments
- Dory (L)

### PIOP/PCS
- Groth16
- Plonk
- Marlin
- FRI

### [Crypto primitives](https://github.com/RustCrypto)
- Pseudorandom generator
- Hashes
- Blake2
- Keccak
- Poseidon
- Pedersen
- Encryption schemes
- AES
- ChaCha20
- Rescue
- ElGamal

### Protocol
- Fiat-Shamir

### Gadgets

## Documentation

To serve the documentation locally, first install both [mdbook](https://rust-lang.github.io/mdBook/guide/installation.html) and the [Katex preprocessor](https://github.com/lzanini/mdbook-katex#getting-started) to render LaTeX, then run

```
make docs
```

## 📊 Benchmarks

Benchmark results are hosted [here](https://lambdaclass.github.io/lambdaworks/bench).

To run them locally, you will need `cargo-criterion` and `cargo-flamegraph`. Install it with:

```bash
cargo install cargo-criterion
```

Run the complete benchmark suite with:

```bash
make benchmarks
```

Run a specific benchmark suite with `cargo`, for example to run the one for `field`:

```bash
make benchmark BENCH=field
```

You can check the generated HTML report in `target/criterion/reports/index.html`
