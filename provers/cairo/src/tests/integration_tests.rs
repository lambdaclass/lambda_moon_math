use crate::{
    cairo_layout::CairoLayout,
    layouts::plain::air::{generate_cairo_proof, verify_cairo_proof, CairoAIR},
    runner::run::generate_prover_args,
    tests::utils::{cairo0_program_path, test_prove_cairo_program},
    Felt252,
};
use lambdaworks_math::field::fields::fft_friendly::stark_252_prime_field::Stark252PrimeField;
use stark_platinum_prover::{
    debug::validate_trace,
    domain::Domain,
    proof::{
        options::{ProofOptions, SecurityLevel},
        stark::StarkProof,
    },
    traits::AIR,
    transcript::StoneProverTranscript,
};

#[test_log::test]
fn test_prove_cairo_simple_program() {
    let layout = CairoLayout::Plain;
    test_prove_cairo_program(&cairo0_program_path("simple_program.json"), layout);
}

#[test_log::test]
fn test_prove_cairo_fibonacci_5() {
    let layout = CairoLayout::Plain;
    test_prove_cairo_program(&cairo0_program_path("fibonacci_5.json"), layout);
}

// #[test_log::test]
// fn test_prove_cairo_fibonacci_5_from_trace() {
//     test_prove_cairo_program_from_trace(
//         &cairo0_program_path("fibonacci_5_trace.bin"),
//         &cairo0_program_path("fibonacci_5_memory.bin"),
//     );
// }

#[test_log::test]
fn test_verifier_rejects_wrong_authentication_paths() {
    // Setup
    let proof_options = ProofOptions::default_test_options();
    let program_content = std::fs::read(cairo0_program_path("fibonacci_5.json")).unwrap();
    let (mut main_trace, pub_inputs) =
        generate_prover_args(&program_content, CairoLayout::Plain).unwrap();

    // Generate the proof
    let mut proof = generate_cairo_proof(&mut main_trace, &pub_inputs, &proof_options).unwrap();

    // Change order of authentication path hashes
    let query = 0;
    let mut original_path = proof.deep_poly_openings[query]
        .main_trace_polys
        .proof
        .merkle_path
        .clone();
    original_path.swap(0, 1);
    // For the test to make sense, we have to make sure
    // that the two hashes are different.
    assert_ne!(original_path[0], original_path[1]);
    proof.deep_poly_openings[query]
        .main_trace_polys
        .proof
        .merkle_path = original_path;

    // Verifier should reject the proof
    assert!(!verify_cairo_proof(&proof, &pub_inputs, &proof_options));
}

#[ignore = "too much time"]
#[test_log::test]
fn test_prove_cairo_fibonacci_1000() {
    let layout = CairoLayout::Plain;
    test_prove_cairo_program(&cairo0_program_path("fibonacci_1000.json"), layout);
}

// #[cfg_attr(feature = "metal", ignore)]
// #[test_log::test]
// fn test_prove_cairo_fibonacci_casm() {
//     let layout = CairoLayout::Plain;
//     test_prove_cairo1_program(&cairo1_program_path("fibonacci_cairo1_mod.casm"), layout);
// }

#[test_log::test]
fn test_verifier_rejects_proof_of_a_slightly_different_program() {
    let program_content = std::fs::read(cairo0_program_path("simple_program.json")).unwrap();
    let (mut main_trace, mut pub_input) =
        generate_prover_args(&program_content, CairoLayout::Plain).unwrap();

    let proof_options = ProofOptions::default_test_options();

    let proof = generate_cairo_proof(&mut main_trace, &pub_input, &proof_options).unwrap();

    // We modify the original program and verify using this new "corrupted" version
    let mut corrupted_program = pub_input.public_memory.clone();
    corrupted_program.insert(Felt252::one(), Felt252::from(6));
    corrupted_program.insert(Felt252::from(3), Felt252::from(5));

    // Here we use the corrupted version of the program in the public inputs
    pub_input.public_memory = corrupted_program;
    assert!(!verify_cairo_proof(&proof, &pub_input, &proof_options));
}

#[test_log::test]
fn test_verifier_rejects_proof_with_different_range_bounds() {
    let program_content = std::fs::read(cairo0_program_path("simple_program.json")).unwrap();
    let (mut main_trace, mut pub_inputs) =
        generate_prover_args(&program_content, CairoLayout::Plain).unwrap();

    let proof_options = ProofOptions::default_test_options();
    let proof = generate_cairo_proof(&mut main_trace, &pub_inputs, &proof_options).unwrap();

    pub_inputs.range_check_min = Some(pub_inputs.range_check_min.unwrap() + 1);
    assert!(!verify_cairo_proof(&proof, &pub_inputs, &proof_options));

    pub_inputs.range_check_min = Some(pub_inputs.range_check_min.unwrap() - 1);
    pub_inputs.range_check_max = Some(pub_inputs.range_check_max.unwrap() - 1);
    assert!(!verify_cairo_proof(&proof, &pub_inputs, &proof_options));
}

#[test_log::test]
fn test_verifier_rejects_proof_with_different_security_params() {
    let program_content = std::fs::read(cairo0_program_path("fibonacci_5.json")).unwrap();
    let (mut main_trace, pub_inputs) =
        generate_prover_args(&program_content, CairoLayout::Plain).unwrap();

    let proof_options_prover = ProofOptions::new_secure(SecurityLevel::Conjecturable80Bits, 3);

    let proof = generate_cairo_proof(&mut main_trace, &pub_inputs, &proof_options_prover).unwrap();

    let proof_options_verifier = ProofOptions::new_secure(SecurityLevel::Conjecturable128Bits, 3);

    assert!(!verify_cairo_proof(
        &proof,
        &pub_inputs,
        &proof_options_verifier
    ));
}

#[test_log::test]
fn check_simple_cairo_trace_evaluates_to_zero() {
    let program_content = std::fs::read(cairo0_program_path("simple_program.json")).unwrap();
    let (mut trace, public_input) =
        generate_prover_args(&program_content, CairoLayout::Plain).unwrap();
    let main_trace_polys = trace.compute_trace_polys_main::<Stark252PrimeField>();
    let mut transcript = StoneProverTranscript::new(&[]);

    let proof_options = ProofOptions::default_test_options();
    let cairo_air = CairoAIR::new(trace.num_rows(), &public_input, &proof_options);
    let rap_challenges = cairo_air.build_rap_challenges(&mut transcript);

    cairo_air.build_auxiliary_trace(&mut trace, &rap_challenges);

    let aux_trace_polys = trace.compute_trace_polys_aux::<Stark252PrimeField>();

    let domain = Domain::new(&cairo_air);

    assert!(validate_trace(
        &cairo_air,
        &main_trace_polys,
        &aux_trace_polys,
        &domain,
        &rap_challenges
    ));
}

#[test]
fn deserialize_and_verify() {
    let program_content = std::fs::read(cairo0_program_path("fibonacci_10.json")).unwrap();
    let (mut main_trace, pub_inputs) =
        generate_prover_args(&program_content, CairoLayout::Plain).unwrap();

    let proof_options = ProofOptions::default_test_options();

    // The proof is generated and serialized.
    let proof = generate_cairo_proof(&mut main_trace, &pub_inputs, &proof_options).unwrap();
    let proof_bytes: Vec<u8> = serde_cbor::to_vec(&proof).unwrap();

    // The trace and original proof are dropped to show that they are decoupled from
    // the verifying process.
    drop(main_trace);
    drop(proof);

    // At this point, the verifier only knows about the serialized proof, the proof options
    // and the public inputs.
    let proof: StarkProof<Stark252PrimeField, Stark252PrimeField> =
        serde_cbor::from_slice(&proof_bytes).unwrap();

    // The proof is verified successfully.
    assert!(verify_cairo_proof(&proof, &pub_inputs, &proof_options));
}
