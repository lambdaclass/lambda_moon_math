use criterion::{
    black_box, criterion_group, criterion_main, measurement::WallTime, BenchmarkGroup, Criterion,
};
use lambdaworks_math::{
    field::fields::fft_friendly::stark_252_prime_field::Stark252PrimeField, traits::Deserializable,
};
use platinum_prover::layouts::plain::air::{verify_cairo_proof, PublicInputs};
use stark_platinum_prover::proof::{
    options::{ProofOptions, SecurityLevel},
    stark::StarkProof,
};
pub mod functions;

fn load_proof_and_pub_inputs(
    input_path: &str,
) -> (
    StarkProof<Stark252PrimeField, Stark252PrimeField>,
    PublicInputs,
) {
    let program_content = std::fs::read(input_path).unwrap();
    let mut bytes = program_content.as_slice();
    let proof_len = usize::from_be_bytes(bytes[0..8].try_into().unwrap());
    bytes = &bytes[8..];
    let proof: StarkProof<Stark252PrimeField, Stark252PrimeField> =
        serde_cbor::from_slice(&bytes[0..proof_len]).unwrap();
    bytes = &bytes[proof_len..];

    let public_inputs = PublicInputs::deserialize(bytes).unwrap();

    (proof, public_inputs)
}

fn verifier_benches(c: &mut Criterion) {
    #[cfg(feature = "parallel")]
    {
        let num_threads: usize = std::env::var("NUM_THREADS")
            .unwrap_or("8".to_string())
            .parse()
            .unwrap();
        println!("Running benchmarks using {} threads", num_threads);
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .unwrap();
    };

    let mut group = c.benchmark_group("VERIFIER");
    run_verifier_bench(
        &mut group,
        "fibonacci/500",
        &cairo0_proof_path("fibonacci_500.proof"),
    );
    run_verifier_bench(
        &mut group,
        "fibonacci/1000",
        &cairo0_proof_path("fibonacci_1000.proof"),
    );
}

fn cairo0_proof_path(program_name: &str) -> String {
    const CARGO_DIR: &str = env!("CARGO_MANIFEST_DIR");
    const PROGRAM_BASE_REL_PATH: &str = "/benches/proofs/";
    let program_base_path = CARGO_DIR.to_string() + PROGRAM_BASE_REL_PATH;
    program_base_path + program_name
}

fn run_verifier_bench(
    group: &mut BenchmarkGroup<'_, WallTime>,
    benchname: &str,
    program_path: &str,
) {
    let (proof, pub_inputs) = load_proof_and_pub_inputs(program_path);
    let proof_options = ProofOptions::new_secure(SecurityLevel::Conjecturable80Bits, 3);
    group.bench_function(benchname, |bench| {
        bench.iter(|| {
            assert!(black_box(verify_cairo_proof(
                &proof,
                &pub_inputs,
                &proof_options
            )))
        });
    });
}

criterion_group!(benches, verifier_benches);
criterion_main!(benches);
