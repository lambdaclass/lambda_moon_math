use super::layouts::plain::air::CairoAIR;
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::fft_friendly::stark_252_prime_field::Stark252PrimeField;
use serde::{Deserialize, Serialize};
use stark_platinum_prover::proof::options::ProofOptions;
use stark_platinum_prover::proof::options::SecurityLevel;
use stark_platinum_prover::proof::stark::StarkProof;
use stark_platinum_prover::transcript::StoneProverTranscript;
use stark_platinum_prover::verifier::{IsStarkVerifier, Verifier};
use std::collections::HashMap;
use wasm_bindgen::prelude::wasm_bindgen;

#[wasm_bindgen]
pub struct Stark252PrimeFieldProof(StarkProof<Stark252PrimeField, Stark252PrimeField>);

#[wasm_bindgen]
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Eq, PartialEq, Hash)]
pub struct FE(FieldElement<Stark252PrimeField>);

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MemoryMap(HashMap<FE, FE>);

#[wasm_bindgen]
pub fn verify_cairo_proof_wasm(proof_bytes: &[u8], proof_options: &ProofOptions) -> bool {
    let bytes = proof_bytes;

    // This logic is the same as main verify, with only error handling changing. In wasm, we simply return a false if the proof is invalid, instead of rising an error.

    // Proof len was stored as an u32, 4u8 needs to be read
    let proof_len = u32::from_le_bytes(bytes[0..4].try_into().unwrap()) as usize;

    let bytes = &bytes[4..];
    if bytes.len() < proof_len {
        return false;
    }

    let Ok((proof, _)) =
        bincode::serde::decode_from_slice(&bytes[0..proof_len], bincode::config::standard())
    else {
        return false;
    };
    let bytes = &bytes[proof_len..];

    let Ok((pub_inputs, _)) = bincode::serde::decode_from_slice(bytes, bincode::config::standard())
    else {
        return false;
    };

    Verifier::<CairoAIR>::verify(
        &proof,
        &pub_inputs,
        proof_options,
        StoneProverTranscript::new(&[]),
    )
}

/// WASM Function for verifying a proof with default 100 bits of security
#[wasm_bindgen]
pub fn verify_cairo_proof_wasm_100_bits(proof_bytes: &[u8]) -> bool {
    let proof_options = ProofOptions::new_secure(SecurityLevel::Conjecturable100Bits, 3);
    verify_cairo_proof_wasm(proof_bytes, &proof_options)
}

#[wasm_bindgen]
pub fn new_proof_options(
    blowup_factor: u8,
    fri_number_of_queries: usize,
    coset_offset: usize,
    grinding_factor: u8,
) -> ProofOptions {
    ProofOptions {
        blowup_factor,
        fri_number_of_queries,
        coset_offset: coset_offset as u64,
        grinding_factor,
    }
}
