use lambdaworks_math::field::{
    element::FieldElement, fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
};

pub mod constraints;
pub mod context;
pub mod debug;
pub mod domain;
pub mod examples;
pub mod frame;
pub mod fri;
pub mod grinding;
pub mod proof;
pub mod prover;
pub mod table;
pub mod trace;
pub mod traits;
pub mod transcript;
pub mod utils;
pub mod verifier;

#[cfg(test)]
pub mod tests;

/// Configurations of the Prover available in compile time
pub mod config;

pub type PrimeField = Stark252PrimeField;
pub type Felt252 = FieldElement<PrimeField>;
