use std::time::Duration;

use lambdaworks_math::field::{element::FieldElement, fields::mersenne31::{extension::{Degree2ExtensionField, Degree4ExtensionField}, field::Mersenne31Field}};

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use p3_mersenne_31::Mersenne31;
use p3_field::extension::{Complex, BinomialExtensionField};

use rand::random;

pub mod utils;

pub type FpE = FieldElement<Mersenne31Field>;
pub type Fp2E = FieldElement<Degree2ExtensionField>;
pub type Fp4E = FieldElement<Degree4ExtensionField>;

pub type PlonkyFpE = Mersenne31;
pub type PlonkyFp2E = Complex<Mersenne31>;
pub type PlonkyFp4E = BinomialExtensionField<Complex<Mersenne31>, 2>;


const BENCHMARK_NAME: &str = "Mersenne 31 Degree 4 extension";

fn create_random_lambdaworks_vector(num: usize) -> Vec<(Fp4E, Fp4E)> {
    let mut result = Vec::with_capacity(num);
    for _ in 0..result.capacity() {
        result.push((
            Fp4E::new([
                Fp2E::new([FpE::new(random()), FpE::new(random())]), 
                Fp2E::new([FpE::new(random()), FpE::new(random())])]), 
            Fp4E::new([
                Fp2E::new([FpE::new(random()), FpE::new(random())]), 
                Fp2E::new([FpE::new(random()), FpE::new(random())])])
        ));
    }
    result
}

fn to_plonk3_vec(lambdaworks_vec: &Vec<(Fp4E, Fp4E)>) -> Vec<(PlonkyFp4E, PlonkyFp4E)> {
    let mut p3_vec = Vec::with_capacity(lambdaworks_vec.len());
    for lambdaworks_pair in lambdaworks_vec {
        let a: PlonkyFp4E = PlonkyFp4E::new(
            PlonkyFp2E::new(
                PlonkyFpE::new(*lambdaworks_pair.0.value()[0].value()[0].value()), 
                PlonkyFpE::new(*lambdaworks_pair.0.value()[0].value()[1].value())
            ),
            PlonkyFp2E::new(
                PlonkyFpE::new(*lambdaworks_pair.0.value()[1].value()[0].value()),
                PlonkyFpE::new(*lambdaworks_pair.0.value()[1].value()[1].value())
            )
        );
        let b: PlonkyFp4E = PlonkyFp4E::new(
            PlonkyFp2E::new(
                PlonkyFpE::new(*lambdaworks_pair.1.value()[0].value()[0].value()), 
                PlonkyFpE::new(*lambdaworks_pair.1.value()[0].value()[1].value())
            ),
            PlonkyFp2E::new(
                PlonkyFpE::new(*lambdaworks_pair.1.value()[1].value()[0].value()),
                PlonkyFpE::new(*lambdaworks_pair.1.value()[1].value()[1].value())
            )
        );
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
                        black_box(black_box(a.0.clone()) + (black_box(a.1.clone())));
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
                        black_box(black_box(a.0.clone()) + (black_box(a.1.clone())));
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
                        black_box(black_box(a.0.clone()) * (black_box(a.1.clone())));
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
                        black_box(black_box(a.0.clone()) * (black_box(a.1.clone())));
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
                        black_box(black_box(a.0.clone()) / (black_box(a.1.clone())));
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
                        black_box(black_box(a.0.clone()) / (black_box(a.1.clone())));
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
        .measurement_time(Duration::from_secs(40))
        .sample_size(200);
    targets = criterion_benchmark
}
criterion_main!(benches);
