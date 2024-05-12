use criterion::{black_box, criterion_group, criterion_main, Criterion};
use lambdaworks_crypto::hash::poseidon::starknet::PoseidonCairoStark252;
use lambdaworks_crypto::hash::poseidon::Poseidon;
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::fft_friendly::stark_252_prime_field::Stark252PrimeField;
use lambdaworks_math::traits::ByteConversion;
use pathfinder_crypto::MontFelt;
use rand::{Rng, RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;

fn poseidon_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("poseidon");

    let mut rng = ChaCha8Rng::seed_from_u64(1);

    group.bench_function("lambdaworks", |b| {
        b.iter_batched(
            || {
                let mut lw_felt_x_bytes: [u8; 32] = Default::default();
                rng.fill_bytes(&mut lw_felt_x_bytes);
                let mut lw_felt_y_bytes: [u8; 32] = Default::default();
                rng.fill_bytes(&mut lw_felt_y_bytes);
                let lw_x =
                    FieldElement::<Stark252PrimeField>::from_bytes_be(&lw_felt_x_bytes).unwrap();
                let lw_y =
                    FieldElement::<Stark252PrimeField>::from_bytes_be(&lw_felt_y_bytes).unwrap();
                (lw_x, lw_y)
            },
            |(a, b)| black_box(PoseidonCairoStark252::hash(&a, &b)),
            criterion::BatchSize::SmallInput,
        )
    });

    group.bench_function("starknet-rs", |b| {
        b.iter_batched(
            || {
                let mut mont_x: [u64; 4] = Default::default();
                rng.fill(&mut mont_x);
                let mut mont_y: [u64; 4] = Default::default();
                rng.fill(&mut mont_y);
                let sn_ff_x = starknet_ff::FieldElement::from_mont(mont_x);
                let sn_ff_y = starknet_ff::FieldElement::from_mont(mont_y);

                (sn_ff_x, sn_ff_y)
            },
            |(a, b)| black_box(starknet_crypto::poseidon_hash(a, b)),
            criterion::BatchSize::SmallInput,
        )
    });

    group.bench_function("pathfinder", |b| {
        b.iter_batched(
            || (MontFelt::random(&mut rng), MontFelt::random(&mut rng)),
            |(a, b)| black_box(pathfinder_crypto::hash::poseidon_hash(a, b)),
            criterion::BatchSize::SmallInput,
        )
    });
}
criterion_group!(poseidon, poseidon_benchmarks);
criterion_main!(poseidon);
