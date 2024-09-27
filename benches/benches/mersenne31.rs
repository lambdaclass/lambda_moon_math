use criterion::{black_box, criterion_group, criterion_main, Criterion};
use lambdaworks_math::field::{
    element::FieldElement,
    fields::mersenne31::{
        extension::{Degree2ExtensionField, Degree4ExtensionField},
        field::Mersenne31Field,
    },
};
use rand::random;
use std::{ops::Add, time::Duration};

pub mod utils;

const BENCHMARK_NAME: &str = "add";

pub fn criterion_benchmark(c: &mut Criterion) {}

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
