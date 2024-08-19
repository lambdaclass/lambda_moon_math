use criterion::{black_box, Criterion};
use lambdaworks_math::{
    cyclic_group::IsGroup,
    elliptic_curve::{
        short_weierstrass::curves::bls12_381::{
            curve::BLS12381Curve,
            default_types::{FrElement, FrField},
        },
        traits::IsEllipticCurve,
        wnaf::WnafTable,
    },
    unsigned_integer::element::U256,
};
use rand::{Rng, SeedableRng};

#[allow(dead_code)]
pub fn wnaf_bls12_381_benchmarks(c: &mut Criterion) {
    let scalar_size = 1000;

    let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(9001);
    let mut scalars = Vec::new();
    for _i in 0..scalar_size {
        scalars.push(FrElement::new(U256::from(rng.gen::<u128>())));
    }

    let g = BLS12381Curve::generator();

    let mut group = c.benchmark_group("BLS12-381 WNAF");
    group.significance_level(0.1).sample_size(100);
    group.throughput(criterion::Throughput::Elements(1));

    group.bench_function(
        format!(
            "Naive BLS12-381 vector multiplication with size {}",
            scalar_size
        ),
        |bencher| {
            bencher.iter(|| {
                black_box(
                    scalars
                        .clone()
                        .iter()
                        .map(|scalar| {
                            black_box(
                                black_box(g.clone())
                                    .operate_with_self(black_box(scalar.clone().representative()))
                                    .to_affine(),
                            )
                        })
                        .collect::<Vec<_>>(),
                )
            });
        },
    );

    group.bench_function(
        format!(
            "WNAF BLS12-381 vector multiplication with size {}",
            scalar_size
        ),
        |bencher| {
            bencher.iter(|| {
                black_box(
                    black_box(WnafTable::<BLS12381Curve, FrField>::new(
                        black_box(&g.clone()),
                        scalar_size,
                    ))
                    .multi_scalar_mul(&black_box(scalars.clone())),
                )
            });
        },
    );
}
