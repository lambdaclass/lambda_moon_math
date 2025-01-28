use std::fs;

use crate::*;
use lambdaworks_groth16_bn_254::*;

const TEST_DIR: &str = "test_files";

// Proves & verifies a Poseidon circuit with 1 input and 2 outputs. The input is decimal 100.
// Converts the following circuit: https://github.com/iden3/circomlib/blob/master/circuits/poseidon.circom
// Poseidon constants: https://github.com/iden3/circomlib/blob/master/circuits/poseidon_constants.circom
#[test]
fn poseidon_parse_prove_verify() {
    let test_dir = format!("{TEST_DIR}/poseidon");

    let (qap, w) = circom_to_lambda(
        &fs::read_to_string(format!("{test_dir}/test.r1cs.json")).expect("Error reading the file"),
        &fs::read_to_string(format!("{test_dir}/witness.json")).expect("Error reading the file"),
    );

    let (pk, vk) = setup(&qap);

    let accept = verify(
        &vk,
        &Prover::prove(&w, &qap, &pk),
        &w[..qap.num_of_public_inputs],
    );
    assert!(accept);
}

// Converts following Circom circuit and inputs into Lambdaworks-compatible QAP and witness assignments
//
// template Test() {
//	signal input x;
//	signal output out;
//
//	signal sym_1;
//	signal y;
//
//	sym_1 <== x * x;
//	out <== (sym_1 * x) + (x + 5);
// }
//
// { "x": 3 }
#[test]
fn vitalik_w_and_qap() {
    let test_dir = format!("{TEST_DIR}/vitalik_example");

    let (qap, w) = circom_to_lambda(
        &fs::read_to_string(format!("{test_dir}/test.r1cs.json")).expect("Error reading the file"),
        &fs::read_to_string(format!("{test_dir}/witness.json")).expect("Error reading the file"),
    );

    // Circom witness contains outputs before circuit inputs where Lambdaworks puts inputs before the output. Freshly generated
    // witness assignment "w" must be in form ["1", "x", "~out", "sym_1"]
    assert_eq!(
        w,
        ["1", "3", "23", "9"]
            .map(FrElement::from_hex_unchecked)
            .to_vec()
    );

    // Regarding QAP, we expect 2 constraints in following form
    // -sym_1 = -x * x
    // x + 5 - ~out = -sym_1 * x
    //
    // Same ordering difference exists for variable matrices, too. Circom adapter changes the
    // order of the rows in the same way it rearranges the witness ordering.

    const R_MINUS_ONE_STR: &str =
        "30644E72E131A029B85045B68181585D2833E84879B9709143E1F593F0000000";
    const ZERO_STR: &str = "0x0";
    const ONE_STR: &str = "0x1";

    let expected_num_of_public_inputs = 1;
    let [temp_l, temp_r, temp_o] = [
        [
            [ZERO_STR, ZERO_STR],        // 1
            [R_MINUS_ONE_STR, ZERO_STR], // x
            [ZERO_STR, ZERO_STR],        // ~out
            [ZERO_STR, R_MINUS_ONE_STR], // sym_1
        ],
        [
            [ZERO_STR, ZERO_STR], // 1
            [ONE_STR, ONE_STR],   // x
            [ZERO_STR, ZERO_STR], // ~out
            [ZERO_STR, ZERO_STR], // sym_1
        ],
        [
            [ZERO_STR, "5"],             // 1
            [ZERO_STR, ONE_STR],         // x
            [ZERO_STR, R_MINUS_ONE_STR], // ~out
            [R_MINUS_ONE_STR, ZERO_STR], // sym_1
        ],
    ]
    .map(|matrix| matrix.map(|row| row.map(FrElement::from_hex_unchecked).to_vec()));
    let expected_qap =
        QAP::from_variable_matrices(expected_num_of_public_inputs, &temp_l, &temp_r, &temp_o);

    let expected_l = expected_qap.l;
    let expected_r = expected_qap.r;
    let expected_o = expected_qap.o;

    assert_eq!(qap.l, expected_l);
    assert_eq!(qap.r, expected_r);
    assert_eq!(qap.o, expected_o);
}

#[test]
fn fibonacci_verify() {
    // Define the directory containing the R1CS and witness files
    let test_dir = format!("{TEST_DIR}/fibonacci");

    // Step 1: Parse R1CS and witness from JSON files
    let (qap, w) = circom_to_lambda(
        &fs::read_to_string(format!("{test_dir}/fibonacci.r1cs.json"))
            .expect("Error reading the R1CS file"),
        &fs::read_to_string(format!("{test_dir}/witness.json"))
            .expect("Error reading the witness file"),
    );

    // Step 2: Generate the proving and verifying keys using the QAP
    let (pk, vk) = setup(&qap);

    // Step 3: Generate the proof using the proving key and witness
    let proof = Prover::prove(&w, &qap, &pk);

    // Step 4: Verify the proof using the verifying key and public inputs
    let accept = verify(&vk, &proof, &w[..qap.num_of_public_inputs]);

    // Step 5: Assert that the verification is successful
    assert!(accept, "Proof verification failed.");

    println!("Proof verification succeeded. All steps completed.");
}
