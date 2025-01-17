use lambdaworks_math::field::{
    element::FieldElement, fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
};

use crate::{
    examples::{
        bit_flags::{self, BitFlagsAIR},
        dummy_air::{self, DummyAIR},
        fibonacci_2_cols_shifted::{self, Fibonacci2ColsShifted},
        fibonacci_2_columns::{self, Fibonacci2ColsAIR},
        fibonacci_rap::{fibonacci_rap_trace, FibonacciRAP, FibonacciRAPPublicInputs},
        quadratic_air::{self, QuadraticAIR, QuadraticPublicInputs},
        read_only_memory::{sort_rap_trace, ReadOnlyPublicInputs, ReadOnlyRAP},
        simple_fibonacci::{self, FibonacciAIR, FibonacciPublicInputs},
        simple_periodic_cols::{self, SimplePeriodicAIR, SimplePeriodicPublicInputs}, //         simple_periodic_cols::{self, SimplePeriodicAIR, SimplePeriodicPublicInputs},
    },
    proof::options::ProofOptions,
    prover::{IsStarkProver, Prover},
    transcript::StoneProverTranscript,
    verifier::{IsStarkVerifier, Verifier},
    Felt252,
};

#[test_log::test]
fn test_prove_fib() {
    let mut trace = simple_fibonacci::fibonacci_trace([Felt252::from(1), Felt252::from(1)], 8);

    let proof_options = ProofOptions::default_test_options();

    let pub_inputs = FibonacciPublicInputs {
        a0: Felt252::one(),
        a1: Felt252::one(),
    };

    let proof = Prover::<FibonacciAIR<Stark252PrimeField>>::prove(
        &mut trace,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[]),
    )
    .unwrap();
    assert!(Verifier::<FibonacciAIR<Stark252PrimeField>>::verify(
        &proof,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[]),
    ));
}

#[test_log::test]
fn test_prove_simple_periodic_8() {
    let mut trace = simple_periodic_cols::simple_periodic_trace::<Stark252PrimeField>(8);

    let proof_options = ProofOptions::default_test_options();

    let pub_inputs = SimplePeriodicPublicInputs {
        a0: Felt252::one(),
        a1: Felt252::from(8),
    };

    let proof = Prover::<SimplePeriodicAIR<Stark252PrimeField>>::prove(
        &mut trace,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[]),
    )
    .unwrap();
    assert!(Verifier::<SimplePeriodicAIR<Stark252PrimeField>>::verify(
        &proof,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[]),
    ));
}

#[test_log::test]
fn test_prove_simple_periodic_32() {
    let mut trace = simple_periodic_cols::simple_periodic_trace::<Stark252PrimeField>(32);

    let proof_options = ProofOptions::default_test_options();

    let pub_inputs = SimplePeriodicPublicInputs {
        a0: Felt252::one(),
        a1: Felt252::from(32768),
    };

    let proof = Prover::<SimplePeriodicAIR<Stark252PrimeField>>::prove(
        &mut trace,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[]),
    )
    .unwrap();
    assert!(Verifier::<SimplePeriodicAIR<Stark252PrimeField>>::verify(
        &proof,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[]),
    ));
}

#[test_log::test]
fn test_prove_fib_2_cols() {
    let mut trace = fibonacci_2_columns::compute_trace([Felt252::from(1), Felt252::from(1)], 16);

    let proof_options = ProofOptions::default_test_options();
    let pub_inputs = FibonacciPublicInputs {
        a0: Felt252::one(),
        a1: Felt252::one(),
    };

    let proof = Prover::<Fibonacci2ColsAIR<Stark252PrimeField>>::prove(
        &mut trace,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[]),
    )
    .unwrap();

    assert!(Verifier::<Fibonacci2ColsAIR<Stark252PrimeField>>::verify(
        &proof,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[])
    ));
}

#[test_log::test]
fn test_prove_fib_2_cols_shifted() {
    let mut trace = fibonacci_2_cols_shifted::compute_trace(FieldElement::one(), 16);

    let claimed_index = 14;
    let claimed_value = trace.main_table.get_row(claimed_index)[0];
    let proof_options = ProofOptions::default_test_options();

    let pub_inputs = fibonacci_2_cols_shifted::PublicInputs {
        claimed_value,
        claimed_index,
    };

    let proof = Prover::<Fibonacci2ColsShifted<_>>::prove(
        &mut trace,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[]),
    )
    .unwrap();
    assert!(Verifier::<Fibonacci2ColsShifted<_>>::verify(
        &proof,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[])
    ));
}

#[test_log::test]
fn test_prove_quadratic() {
    let mut trace = quadratic_air::quadratic_trace(Felt252::from(3), 32);

    let proof_options = ProofOptions::default_test_options();

    let pub_inputs = QuadraticPublicInputs {
        a0: Felt252::from(3),
    };

    let proof = Prover::<QuadraticAIR<Stark252PrimeField>>::prove(
        &mut trace,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[]),
    )
    .unwrap();
    assert!(Verifier::<QuadraticAIR<Stark252PrimeField>>::verify(
        &proof,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[])
    ));
}

#[test_log::test]
fn test_prove_rap_fib() {
    let steps = 16;
    let mut trace = fibonacci_rap_trace([Felt252::from(1), Felt252::from(1)], steps);

    let proof_options = ProofOptions::default_test_options();

    let pub_inputs = FibonacciRAPPublicInputs {
        steps,
        a0: Felt252::one(),
        a1: Felt252::one(),
    };

    let proof = Prover::<FibonacciRAP<Stark252PrimeField>>::prove(
        &mut trace,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[]),
    )
    .unwrap();
    assert!(Verifier::<FibonacciRAP<Stark252PrimeField>>::verify(
        &proof,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[])
    ));
}

#[test_log::test]
fn test_prove_dummy() {
    let trace_length = 16;
    let mut trace = dummy_air::dummy_trace(trace_length);

    let proof_options = ProofOptions::default_test_options();

    let proof = Prover::<DummyAIR>::prove(
        &mut trace,
        &(),
        &proof_options,
        StoneProverTranscript::new(&[]),
    )
    .unwrap();
    assert!(Verifier::<DummyAIR>::verify(
        &proof,
        &(),
        &proof_options,
        StoneProverTranscript::new(&[])
    ));
}

#[test_log::test]
fn test_prove_bit_flags() {
    let mut trace = bit_flags::bit_prefix_flag_trace(32);
    let proof_options = ProofOptions::default_test_options();

    let proof = Prover::<BitFlagsAIR>::prove(
        &mut trace,
        &(),
        &proof_options,
        StoneProverTranscript::new(&[]),
    )
    .unwrap();

    assert!(Verifier::<BitFlagsAIR>::verify(
        &proof,
        &(),
        &proof_options,
        StoneProverTranscript::new(&[]),
    ));
}

#[test_log::test]
fn test_prove_read_only_memory() {
    let address_col = vec![
        FieldElement::<Stark252PrimeField>::from(3), // a0
        FieldElement::<Stark252PrimeField>::from(2), // a1
        FieldElement::<Stark252PrimeField>::from(2), // a2
        FieldElement::<Stark252PrimeField>::from(3), // a3
        FieldElement::<Stark252PrimeField>::from(4), // a4
        FieldElement::<Stark252PrimeField>::from(5), // a5
        FieldElement::<Stark252PrimeField>::from(1), // a6
        FieldElement::<Stark252PrimeField>::from(3), // a7
    ];
    let value_col = vec![
        FieldElement::<Stark252PrimeField>::from(10), // v0
        FieldElement::<Stark252PrimeField>::from(5),  // v1
        FieldElement::<Stark252PrimeField>::from(5),  // v2
        FieldElement::<Stark252PrimeField>::from(10), // v3
        FieldElement::<Stark252PrimeField>::from(25), // v4
        FieldElement::<Stark252PrimeField>::from(25), // v5
        FieldElement::<Stark252PrimeField>::from(7),  // v6
        FieldElement::<Stark252PrimeField>::from(10), // v7
    ];

    let pub_inputs = ReadOnlyPublicInputs {
        a0: FieldElement::<Stark252PrimeField>::from(3),
        v0: FieldElement::<Stark252PrimeField>::from(10),
        a_sorted0: FieldElement::<Stark252PrimeField>::from(1), // a6
        v_sorted0: FieldElement::<Stark252PrimeField>::from(7), // v6
    };
    let mut trace = sort_rap_trace(address_col, value_col);
    let proof_options = ProofOptions::default_test_options();
    let proof = Prover::<ReadOnlyRAP<Stark252PrimeField>>::prove(
        &mut trace,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[]),
    )
    .unwrap();
    assert!(Verifier::<ReadOnlyRAP<Stark252PrimeField>>::verify(
        &proof,
        &pub_inputs,
        &proof_options,
        StoneProverTranscript::new(&[])
    ));
}
