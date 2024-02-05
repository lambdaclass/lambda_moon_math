use lambdaworks_math::field::{
    element::FieldElement, fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
};

pub mod cairo_layout;
pub mod cairo_mem;
pub mod decode;
pub mod errors;
pub mod execution_trace;
pub mod layouts;
pub mod register_states;
pub mod runner;
pub mod transition_constraints;

#[cfg(test)]
pub mod tests;

#[cfg(feature = "wasm")]
pub mod wasm_wrappers;

pub type PrimeField = Stark252PrimeField;
pub type Felt252 = FieldElement<PrimeField>;
