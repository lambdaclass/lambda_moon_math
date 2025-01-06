use criterion::Criterion;
use std::hint::black_box;

use lambdaworks_math::field::{
    element::FieldElement,
    fields::{
        fft_friendly::babybear::Babybear31PrimeField,
        u32_montgomery_backend_prime_field::U32MontgomeryBackendPrimeField,
    },
};

use p3_baby_bear::BabyBear;
use p3_field::{Field, FieldAlgebra};

use rand::random;
use rand::Rng;

pub type U32Babybear31PrimeField = U32MontgomeryBackendPrimeField<2013265921>;
pub type F = FieldElement<U32Babybear31PrimeField>;
pub type F64 = FieldElement<Babybear31PrimeField>;

pub fn rand_field_elements(num: usize) -> Vec<(F, F)> {
    let mut result = Vec::with_capacity(num);
    for _ in 0..result.capacity() {
        result.push((F::from(random::<u64>()), F::from(random::<u64>())));
    }
    result
}

fn rand_babybear_elements_p3(num: usize) -> Vec<(BabyBear, BabyBear)> {
    let mut rng = rand::thread_rng();
    (0..num)
        .map(|_| (rng.gen::<BabyBear>(), rng.gen::<BabyBear>()))
        .collect()
}

pub fn babybear_ops_benchmarks(c: &mut Criterion) {
    let input: Vec<Vec<(F, F)>> = [1000000]
        .into_iter()
        .map(rand_field_elements)
        .collect::<Vec<_>>();
    let mut group = c.benchmark_group("BabyBear operations using Lambdaworks u32");

    for i in input.clone().into_iter() {
        group.bench_with_input(format!("Addition {:?}", &i.len()), &i, |bench, i| {
            bench.iter(|| {
                for (x, y) in i {
                    black_box(black_box(x) + black_box(y));
                }
            });
        });
    }

    for i in input.clone().into_iter() {
        group.bench_with_input(format!("Multiplication {:?}", &i.len()), &i, |bench, i| {
            bench.iter(|| {
                for (x, y) in i {
                    black_box(black_box(x) * black_box(y));
                }
            });
        });
    }

    for i in input.clone().into_iter() {
        group.bench_with_input(format!("Square {:?}", &i.len()), &i, |bench, i| {
            bench.iter(|| {
                for (x, _) in i {
                    black_box(black_box(x).square());
                }
            });
        });
    }

    for i in input.clone().into_iter() {
        group.bench_with_input(format!("Inverse {:?}", &i.len()), &i, |bench, i| {
            bench.iter(|| {
                for (x, _) in i {
                    black_box(black_box(x).inv().unwrap());
                }
            });
        });
    }

    for i in input.clone().into_iter() {
        group.bench_with_input(format!("Division {:?}", &i.len()), &i, |bench, i| {
            bench.iter(|| {
                for (x, y) in i {
                    black_box(black_box(x) / black_box(y));
                }
            });
        });
    }
}

pub fn rand_field_elements_u64(num: usize) -> Vec<(F64, F64)> {
    let mut result = Vec::with_capacity(num);
    for _ in 0..result.capacity() {
        result.push((F64::from(random::<u64>()), F64::from(random::<u64>())));
    }
    result
}
pub fn babybear_ops_benchmarks_f64(c: &mut Criterion) {
    let input: Vec<Vec<(F64, F64)>> = [1000000]
        .into_iter()
        .map(rand_field_elements_u64)
        .collect::<Vec<_>>();
    let mut group = c.benchmark_group("BabyBear operations using Lambdaworks u64");

    for i in input.clone().into_iter() {
        group.bench_with_input(format!("Addition {:?}", &i.len()), &i, |bench, i| {
            bench.iter(|| {
                for (x, y) in i {
                    black_box(black_box(x) + black_box(y));
                }
            });
        });
    }

    for i in input.clone().into_iter() {
        group.bench_with_input(format!("Multiplication {:?}", &i.len()), &i, |bench, i| {
            bench.iter(|| {
                for (x, y) in i {
                    black_box(black_box(x) * black_box(y));
                }
            });
        });
    }

    for i in input.clone().into_iter() {
        group.bench_with_input(format!("Square {:?}", &i.len()), &i, |bench, i| {
            bench.iter(|| {
                for (x, _) in i {
                    black_box(black_box(x).square());
                }
            });
        });
    }

    for i in input.clone().into_iter() {
        group.bench_with_input(format!("Inverse {:?}", &i.len()), &i, |bench, i| {
            bench.iter(|| {
                for (x, _) in i {
                    black_box(black_box(x).inv().unwrap());
                }
            });
        });
    }

    for i in input.clone().into_iter() {
        group.bench_with_input(format!("Division {:?}", &i.len()), &i, |bench, i| {
            bench.iter(|| {
                for (x, y) in i {
                    black_box(black_box(x) / black_box(y));
                }
            });
        });
    }
}

pub fn babybear_p3_ops_benchmarks(c: &mut Criterion) {
    let input: Vec<Vec<(BabyBear, BabyBear)>> = [1000000]
        .into_iter()
        .map(rand_babybear_elements_p3)
        .collect::<Vec<_>>();

    let mut group = c.benchmark_group("BabyBear operations using Plonky3");

    for i in input.clone().into_iter() {
        group.bench_with_input(format!("Addition {:?}", &i.len()), &i, |bench, i| {
            bench.iter(|| {
                for (x, y) in i {
                    black_box(black_box(*x) + black_box(*y));
                }
            });
        });
    }

    for i in input.clone().into_iter() {
        group.bench_with_input(format!("Multiplication {:?}", &i.len()), &i, |bench, i| {
            bench.iter(|| {
                for (x, y) in i {
                    black_box(black_box(*x) * black_box(*y));
                }
            });
        });
    }

    for i in input.clone().into_iter() {
        group.bench_with_input(format!("Square {:?}", &i.len()), &i, |bench, i| {
            bench.iter(|| {
                for (x, _) in i {
                    black_box(black_box(x).square());
                }
            });
        });
    }
    for i in input.clone().into_iter() {
        group.bench_with_input(format!("Inverse {:?}", &i.len()), &i, |bench, i| {
            bench.iter(|| {
                for (x, _) in i {
                    black_box(black_box(x).inverse());
                }
            });
        });
    }

    for i in input.clone().into_iter() {
        group.bench_with_input(format!("Division {:?}", &i.len()), &i, |bench, i| {
            bench.iter(|| {
                for (x, y) in i {
                    black_box(black_box(*x) / black_box(*y));
                }
            });
        });
    }
}
