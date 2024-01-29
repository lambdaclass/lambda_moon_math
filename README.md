# lambdaworks
> From the heights of these towers of fields, forty centuries of mathematics look down on us. 

This library provides efficient implementation of cryptographic primitives used to build proving systems. Along with it, many backends for proving systems are shipped, and compatibility with different frontends is supported.

- [Transforming the Future with Zero-Knowledge Proofs, Fully Homomorphic Encryption and new Distributed Systems algorithms](https://blog.lambdaclass.com/transforming-the-future-with-zero-knowledge-proofs-fully-homomorphic-encryption-and-new-distributed-systems-algorithms/)
- [Lambda Crypto Doctrine](https://blog.lambdaclass.com/lambda-crypto-doctrine/)

[![Telegram Chat][tg-badge]][tg-url]
[![codecov](https://img.shields.io/codecov/c/github/lambdaclass/lambdaworks)](https://codecov.io/gh/lambdaclass/lambdaworks)

[tg-badge]: https://img.shields.io/endpoint?url=https%3A%2F%2Ftg.sumanjay.workers.dev%2Flambdaworks%2F&logo=telegram&label=chat&color=neon
[tg-url]: https://t.me/lambdaworks

</div>

## Why we built lambdaworks

Zero-Knowledge and Validity Proofs have gained a lot of attention over the last few years. We strongly believe in this potential and that is why we decided to start working in this challenging ecosystem, where math, cryptography and distributed systems meet. The main barrier in the beginning was not the cryptography or math but the lack of good libraries which are performant and developer friendly. There are some exceptions, though, like gnark or halo2. Some have nice APIs and are easy to work with, but they are not written in Rust, and some are written in Rust but have poor programming and engineering practices. Most of them don't have support for CUDA, Metal and WebGPU or distributed FFT calculation using schedulers like Dask.

So, we decided to build our library, focusing on performance, with clear documentation and developer-focused. Our core team is a group of passionate people from different backgrounds and different strengths; we think that the whole is greater than just the addition of the parts. We don't want to be a compilation of every research result in the ZK space. We want this to be a library that can be used in production, not just in academic research. We want to offer developers the main building blocks and proof systems so that they can build their applications on top of this library.

## [Documentation](https://lambdaclass.github.io/lambdaworks)

## Main crates

- [Math](https://github.com/lambdaclass/lambdaworks/tree/main/math)
- [Crypto primitives](https://github.com/lambdaclass/lambdaworks/tree/main/crypto)
- [STARK Prover](https://github.com/lambdaclass/lambdaworks/tree/main/provers/stark)
- [Plonk Prover](https://github.com/lambdaclass/lambdaworks/tree/main/provers/plonk)
- [Cairo Prover](https://github.com/lambdaclass/lambdaworks/tree/main/provers/cairo)
- [Groth 16](https://github.com/lambdaclass/lambdaworks/tree/main/provers/groth16)

### Crypto
- [Elliptic curves](https://github.com/lambdaclass/lambdaworks/tree/main/math/src/elliptic_curve)
- [Multiscalar multiplication](https://github.com/lambdaclass/lambdaworks/tree/main/math/src/msm)
- [Hashes](https://github.com/lambdaclass/lambdaworks/tree/main/crypto/src/hash)

Most of math and crypto crates supports no-std without allocation with `no-default-features`. A few functions and modules require the `alloc` feature.

Both Math and Crypto support wasm with target `wasm32-unknown-unknown`. To see an example of how to use this to deploy a verifier in a browser, check the Cairo Prover wasm-pack verifier.

## Examples - mini apps
- [Merkle Tree CLI](https://github.com/lambdaclass/lambdaworks/tree/main/examples/merkle-tree-cli)
- [Proving Miden](https://github.com/lambdaclass/lambdaworks/tree/main/examples/prove-miden)
- [Shamir's secret sharing](https://github.com/lambdaclass/lambdaworks/tree/main/examples/shamir_secret_sharing)

## Exercises and Challenges
- [lambdaworks exercises and challenges](https://github.com/lambdaclass/lambdaworks_exercises/tree/main)
- [Roadmap for Sparkling Water Bootcamp](https://github.com/lambdaclass/sparkling_water_bootcamp/blob/main/README.md)

## Citing lambdaworks

If you use ```lambdaworks``` libraries in your research projects, please cite them using the following template:

``` bibtex
@software{lambdaworks,
  author={lambdaworks contributors},
  title={lambdaworks},
  url={https://github.com/lambdaclass/lambdaworks},
  year={2023}
}
```

## List of features

Disclaimer: This list contains cryptographic primitives and mathematical structures that we want to support in lambdaworks. It can be expanded later to include new primitives. If you find there is a mistake or there has been an update in another library, please let us know.

List of symbols:
- :heavy_check_mark: means the feature is currently supported.
- 🏗️ means that the feature is partially implemented or is under active construction.
- :x: means that the feature is not currently supported.

| Finite Fields  | Lambdaworks        | Arkworks           | Halo2    | gnark              | Constantine |
| -------------- | ------------------ | ------------------ | -------- | ------------------ | ----------- |
| StarkField 252 | :heavy_check_mark: | :heavy_check_mark: | :x:      | :heavy_check_mark: | :x:         |
| Mersenne 31    | :heavy_check_mark: | :x:                | :x:      | :x:                | :x:         |
| Baby Bear      | :heavy_check_mark: | :x:                | :x:      | :x:                | :x:         |
| MiniGoldilocks | :heavy_check_mark: | :x:                | :x:      | :heavy_check_mark: | :x:         |
| **ZK friendly Hash function** | **Lambdaworks** | **Arkworks**       | **Halo2**          | **gnark** | **Constantine** |
| Poseidon                      | 🏗️              | :heavy_check_mark: | :heavy_check_mark: | :x:       | :x:             |
| Pedersen                      | 🏗️              | :heavy_check_mark: | :heavy_check_mark: | :x:       | :x:             |
| Rescue Prime XLIX             | :x:             | :x:                | :x:                | :x:       | :x:             |
| **Elliptic Curves** | **Lambdaworks** | **Arkworks**          | **Halo2**          | **gnark**          | **Constantine**    |
| BLS12-381           | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| BLS12-377           | 🏗️                 | :heavy_check_mark: | :x:                | :heavy_check_mark: | :heavy_check_mark: |
| BN-254              | 🏗️                 | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| Pallas              | :heavy_check_mark: | :heavy_check_mark: | :x:                | :x:                | :heavy_check_mark: |
| Vesta               | :heavy_check_mark: | :heavy_check_mark: | :x:                | :x:                | :heavy_check_mark: |
| Bandersnatch        | 🏗️                 | :heavy_check_mark: | :x:                | :heavy_check_mark:  | :heavy_check_mark: |
| **STARKs**       | **Lambdaworks**     | **Arkworks** | **Halo2** | **gnark** | **Constantine** |
| STARK Prover     | :heavy_check_mark:  | :x:          | :x:       | :x:       | :x:             |
| CAIRO Prover     | 🏗️                  | :x:          | :x:       | :x:       | :x:             |
| **SNARKs** | **Lambdaworks**    | **Arkworks**       | **Halo2** | **gnark**          | **Constantine** |
| Groth16    | :heavy_check_mark: | :heavy_check_mark: | :x:       | :heavy_check_mark: | :x:             |
| Plonk      | 🏗️                 | :heavy_check_mark: | ✔️         | :heavy_check_mark: | :x:             |
| Spartan    | :x:                | :heavy_check_mark: | :x:       | :x:                | :x:             |
| Marlin     | :x:                | :heavy_check_mark: | :x:       | :x:                | :x:             |
| GKR        | :x:                | :heavy_check_mark: | :x:       | :heavy_check_mark: | :x:             |
| **Polynomial Commitment Schemes** | **Lambdaworks**    | **Arkworks**       | **Halo2**          | **gnark**          | **Constantine** |
| KZG10                             | :heavy_check_mark: | ✔️                  | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark:             |
| FRI                               | 🏗️                 | :x:                | :x:                | :heavy_check_mark: | :x:             |
| IPA                               | 🏗️                 | ✔️                  | :heavy_check_mark: | :x:                | :x:             |
| Brakedown                         | :x:                | :x: | :x:                | :x:                | :x:             |
| Basefold                          | :x:                | :x: | :x:                | :x:                | :x:             |
| **Folding Schemes** | **Lambdaworks** | **Arkworks**       | **Halo2** | **gnark** | **Constantine** |
| Nova                | :x:             | :heavy_check_mark: | :x:       | :x:       | :x:             |
| Supernova           | :x:             | :x:                | :x:       | :x:       | :x:             |
| Protostar           | :x:             | :x:                | :x:       | :x:       | :x:             |
| Protogalaxy         | :x:             | :heavy_check_mark: | :x:       | :x:       | :x:             |

Additionally, provers are compatible with the following frontends and VMs:

| Backend | Frontend | Status |
|---------|----------|--------|
| Groth16 | Arkworks | :heavy_check_mark: |
| Groth16 | Gnark    | :x: |
| Groth16 | Circom   | 🏗️  |
| Plonk   | Gnark    | 🏗️  |
| Plonk   | Noir    | :x: |
| Stark   | Winterfell | :heavy_check_mark: | 
| Stark   | Miden | :heavy_check_mark: |
| Stark   | Cairo | :heavy_check_mark: |

This can be used in a multi prover setting for extra security, or as a standalone to be used with Rust. 

## Additional tooling usage

### Fuzzers

Fuzzers are divided between the ones that use only the CPU, the ones that use Metal, and the ones that use CUDA.

CPU Fuzzers can be run with the command ```bash make run-fuzzer FUZZER=fuzzer_name```

For example:

```bash
make run-fuzzer FUZZER=field_from_hex
```

The list of fuzzers can be found in `fuzz/no_gpu_fuzz`

Fuzzers for FTT in Metal and Cuda can be run with `make run-metal-fuzzer` and `make run-cuda-fuzzer`


Run a specific fuzzer from the ones contained in **fuzz/fuzz_targets/** folder with`cargo`, for example to run the one for the target `field_from_hex`:

```bash
make run-fuzzer FUZZER=field_from_hex
```

### Documentation building

To serve the documentation locally, first install both [mdbook](https://rust-lang.github.io/mdBook/guide/installation.html) and the [Katex preprocessor](https://github.com/lzanini/mdbook-katex#getting-started) to render LaTeX, then run

``` shell
make docs
```

## 📚 References

The following links, repos and projects have been important in the development of this library and we want to thank and acknowledge them. 

- [Starkware](https://starkware.co/)
- [Winterfell](https://github.com/facebook/winterfell)
- [Anatomy of a Stark](https://aszepieniec.github.io/stark-anatomy/overview)
- [Giza](https://github.com/maxgillett/giza)
- [Ministark](https://github.com/andrewmilson/ministark)
- [Sandstorm](https://github.com/andrewmilson/sandstorm)
- [STARK-101](https://starkware.co/stark-101/)
- [starknet-rs](https://github.com/xJonathanLEI/starknet-rs/)
- [Risc0](https://github.com/risc0/risc0)
- [Neptune](https://github.com/Neptune-Crypto)
- [Summary on FRI low degree test](https://eprint.iacr.org/2022/1216)
- [STARKs paper](https://eprint.iacr.org/2018/046)
- [DEEP FRI](https://eprint.iacr.org/2019/336)
- [BrainSTARK](https://aszepieniec.github.io/stark-brainfuck/)
- [Plonky2](https://github.com/mir-protocol/plonky2)
- [Aztec](https://github.com/AztecProtocol)
- [Arkworks](https://github.com/arkworks-rs)
- [Thank goodness it's FRIday](https://vitalik.ca/general/2017/11/22/starks_part_2.html)
- [Diving DEEP FRI](https://blog.lambdaclass.com/diving-deep-fri/)
- [Periodic constraints](https://blog.lambdaclass.com/periodic-constraints-and-recursion-in-zk-starks/)
- [Chiplets Miden VM](https://wiki.polygon.technology/docs/miden/design/chiplets/main/)
- [Valida](https://github.com/valida-xyz/valida/tree/main)
- [Solidity Verifier](https://github.com/starkware-libs/starkex-contracts/tree/master/evm-verifier/solidity/contracts/cpu)
- [CAIRO verifier](https://github.com/starkware-libs/cairo-lang/tree/master/src/starkware/cairo/stark_verifier)
- [EthSTARK](https://github.com/starkware-libs/ethSTARK/tree/master)
- [CAIRO whitepaper](https://eprint.iacr.org/2021/1063.pdf)
- [Gnark](https://github.com/Consensys/gnark)
