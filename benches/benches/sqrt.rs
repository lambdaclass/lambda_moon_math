use ark_ff::{BigInt, Field};
use ark_std::{test_rng, UniformRand};
use ark_test_curves::starknet_fp::Fq as F;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use lambdaworks_math::unsigned_integer::element::UnsignedInteger;

const BENCHMARK_NAME: &str = "sqrt";

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = test_rng();

    let mut arkworks_vec = Vec::new();
    for _i in 0..1000 {
        let a = F::rand(&mut rng);
        let square = a * a;
        arkworks_vec.push(square);
    }

    // arkworks-ff
    {
        c.bench_function(
            &format!(
                "{} | ark-ff - branch: faster-benchmarks-and-starknet-field",
                BENCHMARK_NAME
            ),
            |b| {
                b.iter(|| {
                    let mut iter = arkworks_vec.iter();

                    for _i in 0..1000 {
                        let a = iter.next().unwrap();
                        black_box(black_box(a).sqrt());
                    }
                });
            },
        );
    }

    // lambdaworks-math
    {
        use lambdaworks_math::field::{
            element::FieldElement, fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
        };
        let mut lambdaworks_vec = Vec::new();
        for arkworks_felt in arkworks_vec {
            let big_int: BigInt<4> = arkworks_felt.into();
            let mut limbs = big_int.0;
            limbs.reverse();

            let a: FieldElement<Stark252PrimeField> =
                FieldElement::from(&UnsignedInteger { limbs });

            assert_eq!(a.representative().limbs, limbs);

            lambdaworks_vec.push(a);
        }

        c.bench_function(&format!("{} | lambdaworks", BENCHMARK_NAME,), |b| {
            b.iter(|| {
                let mut iter = lambdaworks_vec.iter();

                for _i in 0..1000 {
                    let a = iter.next().unwrap();
                    black_box(black_box(a).sqrt());
                }
            });
        });
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
