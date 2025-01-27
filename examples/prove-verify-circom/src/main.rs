use lambdaworks_circom_adapter::*;
use lambdaworks_groth16::*;

fn main() {
    println!("Reading input files");
    let circom_r1cs =
        read_circom_r1cs("./examples/prove-verify-circom/input_files/test.r1cs.json").unwrap();
    let circom_witness =
        read_circom_witness("./examples/prove-verify-circom/input_files/witness.json").unwrap();

    println!("Converting to Lambdaworks-compatible QAP and witness assignments");
    let (qap, witness, _) = circom_to_lambda(circom_r1cs, circom_witness);

    println!("Performing trusted setup");
    let (pk, vk) = setup(&qap);

    println!("Proving");
    let proof = Prover::prove(&witness, &qap, &pk);

    println!("Verifying");
    let accept = verify(&vk, &proof, &witness[..qap.num_of_public_inputs]);

    assert!(accept, "Proof verification failed!");
    println!("Proof verified!");
}
