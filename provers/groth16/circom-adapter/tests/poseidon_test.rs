use lambdaworks_circom_adapter::{circom_to_lambda, read_circom_r1cs, read_circom_witness};

// Proves & verifies a Poseidon circuit with 1 input and 2 outputs. The input is decimal 100.
#[test]
fn poseidon_parse_prove_verify() {
    let circom_r1cs =
        read_circom_r1cs("./tests/poseidon/test.r1cs.json").expect("could not read r1cs");
    let circom_wtns =
        read_circom_witness("./tests/poseidon/witness.json").expect("could not read witness");

    let (qap, wtns, pubs) = circom_to_lambda(circom_r1cs, circom_wtns);

    let (pk, vk) = lambdaworks_groth16::setup(&qap);
    let proof = lambdaworks_groth16::Prover::prove(&wtns, &qap, &pk);
    let accept = lambdaworks_groth16::verify(&vk, &proof, &pubs);
    assert!(accept, "proof verification failed");
}
