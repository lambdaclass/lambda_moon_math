use std::{mem::offset_of, time::Duration};

use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::mersenne31::field::Mersenne31Field;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use p3_mersenne_31::Mersenne31;
use rand::random;

pub mod utils;

pub type F = FieldElement<Mersenne31Field>;

const BENCHMARK_NAME: &str = "Mersenne 31";

fn create_random_lambdaworks_vector(num: usize) -> Vec<(F, F)> {
    let mut result = Vec::with_capacity(num);
    for _ in 0..result.capacity() {
        result.push((F::new(random()), F::new(random())));
    }
    result
}

fn to_plonk3_vec(lambdaworks_vec: &Vec<(F, F)>) -> Vec<(Mersenne31, Mersenne31)> {
    let mut p3_vec = Vec::with_capacity(lambdaworks_vec.len());
    for lambdaworks_pair in lambdaworks_vec {
        let a: Mersenne31 = Mersenne31::new(*lambdaworks_pair.0.value());
        let b: Mersenne31 = Mersenne31::new(*lambdaworks_pair.1.value());

        p3_vec.push((a, b));
    }
    p3_vec
}

// fn to_stwo_vec() -> {}

fn criterion_benchmark(c: &mut Criterion) {
    let lambdaworks_vec = create_random_lambdaworks_vector(1_000_000);
    let plonky_3_vec = to_plonk3_vec(&lambdaworks_vec);

    {
        c.bench_function(
            &format!("{} - Add 1M elements | Lambdaworks", BENCHMARK_NAME),
            |b| {
                b.iter(|| {
                    let mut iter = lambdaworks_vec.iter();

                    for _i in 0..1000000 {
                        let a = iter.next().unwrap();
                        black_box(black_box(a.0) + (black_box(a.1)));
                    }
                });
            },
        );
    }

    {

        c.bench_function(
            &format!("{} - Add 1M elements | Polnky3", BENCHMARK_NAME),
            |b| {
                b.iter(|| {
                    let mut iter = plonky_3_vec.iter();

                    for _i in 0..1000000 {
                        let a = iter.next().unwrap();
                        black_box(black_box(a.0) + (black_box(a.1)));
                    }
                });
            },
        );
    }

    {
        c.bench_function(
            &format!("{} - Mul 1M elements | Lambdaworks", BENCHMARK_NAME),
            |b| {
                b.iter(|| {
                    let mut iter = lambdaworks_vec.iter();

                    for _i in 0..1000000 {
                        let a = iter.next().unwrap();
                        black_box(black_box(a.0) * (black_box(a.1)));
                    }
                });
            },
        );
    }

    {

        c.bench_function(
            &format!("{} - Mul 1M elements | Polnky3", BENCHMARK_NAME),
            |b| {
                b.iter(|| {
                    let mut iter = plonky_3_vec.iter();

                    for _i in 0..1000000 {
                        let a = iter.next().unwrap();
                        black_box(black_box(a.0) * (black_box(a.1)));
                    }
                });
            },
        );
    }

    {
        c.bench_function(
            &format!("{} - Div 1M elements | Lambdaworks", BENCHMARK_NAME),
            |b| {
                b.iter(|| {
                    let mut iter = lambdaworks_vec.iter();

                    for _i in 0..1000000 {
                        let a = iter.next().unwrap();
                        black_box(black_box(a.0) / (black_box(a.1)));
                    }
                });
            },
        );
    }

    {

        c.bench_function(
            &format!("{} - Div 1M elements | Polnky3", BENCHMARK_NAME),
            |b| {
                b.iter(|| {
                    let mut iter = plonky_3_vec.iter();

                    for _i in 0..1000000 {
                        let a = iter.next().unwrap();
                        black_box(black_box(a.0) / (black_box(a.1)));
                    }
                });
            },
        );
    }


}

criterion_group! {
    name = benches;
    // This can be any expression that returns a `Criterion` object.
    config = Criterion::default()
        .significance_level(0.01)
        .measurement_time(Duration::from_secs(15))
        .sample_size(300);
    targets = criterion_benchmark
}
criterion_main!(benches);
