use criterion::{black_box, criterion_group, criterion_main, Criterion};
use lambdaworks_math::{
    cyclic_group::IsGroup,
    elliptic_curve::{
        short_weierstrass::curves::{
            bls12_377::curve::BLS12377Curve, bls12_381::curve::BLS12381Curve,
            bn_254::curve::BN254Curve,
        },
        traits::IsEllipticCurve,
    },
    field::{element::FieldElement, traits::IsField},
};

use lambdaworks_math::gpu::icicle::{
    bls12_377::bls12_377_g1_msm, bls12_381::bls12_381_g1_msm, bn254::bn254_g1_msm,
};
use rand::{rngs::StdRng, Rng, SeedableRng};

pub fn generate_cs_and_points<C: IsEllipticCurve>(
    msm_size: usize,
) -> (Vec<FieldElement<C::BaseField>>, Vec<C::PointRepresentation>)
where
    <C::BaseField as IsField>::BaseType: From<u64>,
{
    // We use a seeded rng so the benchmarks are reproducible.
    let mut rng = StdRng::seed_from_u64(42);

    let g = C::generator();

    let cs: Vec<_> = (0..msm_size)
        .map(|_| FieldElement::<C::BaseField>::new(rng.gen::<u64>().into()))
        .collect();

    let points: Vec<_> = (0..msm_size)
        .map(|_| g.operate_with_self(rng.gen::<u64>()))
        .collect();

    (cs, points)
}

pub fn msm_benchmarks_with_size(c: &mut Criterion, msm_size: usize) {
    let mut group = c.benchmark_group(format!("MSM benchmarks with size {msm_size}"));

    let (cs, points) = generate_cs_and_points::<BLS12381Curve>(msm_size);
    group.bench_function("BLS12_381", |bench| {
        bench.iter(|| black_box(bls12_381_g1_msm(&cs, &points, None)));
    });

    let (cs, points) = generate_cs_and_points::<BLS12377Curve>(msm_size);
    group.bench_function("BLS12_377", |bench| {
        bench.iter(|| black_box(bls12_377_g1_msm(&cs, &points, None)));
    });

    let (cs, points) = generate_cs_and_points::<BN254Curve>(msm_size);
    group.bench_function("BN_254", |bench| {
        bench.iter(|| black_box(bn254_g1_msm(&cs, &points, None)));
    });
}

pub fn run_benchmarks(c: &mut Criterion) {
    let exponents = 1..=18;

    for exp in exponents {
        let msm_size = 1 << exp;

        msm_benchmarks_with_size(c, msm_size);
    }
}

criterion_group!(icicle_msm, run_benchmarks);
criterion_main!(icicle_msm);
