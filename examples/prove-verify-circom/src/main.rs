use std::fs;

use lambdaworks_circom_adapter_bls12_381::*;
use lambdaworks_groth16_bls12_381::*;

const TEST_DIR: &str = "input_files/";

fn main() {
    println!("\nReading input files");
    let r1cs_file_content =
        &fs::read_to_string(format!("{TEST_DIR}test.r1cs.json")).expect("Error reading the file");
    let witness_file_content =
        &fs::read_to_string(format!("{TEST_DIR}witness.json")).expect("Error reading the file");

    println!("\nConverting to Lambdaworks-compatible QAP and witness assignments");
    let (qap, w) = circom_to_lambda(r1cs_file_content, witness_file_content);

    println!("\nPerforming trusted setup");
    let (pk, vk) = setup(&qap);

    println!("\nProving");
    let proof = Prover::prove(&w, &qap, &pk);

    println!("\nVerifying");
    let accept = verify(&vk, &proof, &w[..qap.num_of_public_inputs]);

    assert!(accept);
    println!("Proof verified!");
}
