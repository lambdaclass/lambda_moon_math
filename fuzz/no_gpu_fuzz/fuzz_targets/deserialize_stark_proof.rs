#![no_main]
use lambdaworks_math::field::fields::fft_friendly::stark_252_prime_field::Stark252PrimeField;
use lambdaworks_math::traits::Deserializable;
use libfuzzer_sys::fuzz_target;
use stark_platinum_prover::proof::stark::StarkProof;

fuzz_target!(|data: Vec<u8>| {
    let _proof = StarkProof::<Stark252PrimeField>::deserialize(&data);
});
