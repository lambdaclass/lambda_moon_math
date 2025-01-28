#[cfg(test)]
mod integration_tests;

use lambdaworks_groth16_bls12_381::{common::FrElement, QuadraticArithmeticProgram as QAP};
use lambdaworks_math::unsigned_integer::element::UnsignedInteger;
use serde_json::Value;

pub fn circom_to_lambda(
    r1cs_file_content: &str,
    witness_file_content: &str,
) -> (QAP, Vec<FrElement>) {
    let circom_r1cs: Value = serde_json::from_str(r1cs_file_content).expect("Error parsing JSON");
    let [mut l, mut r, mut o] = build_lro_from_circom_r1cs(&circom_r1cs);

    let mut witness: Vec<_> = serde_json::from_str::<Vec<String>>(witness_file_content)
        .expect("Error parsing JSON")
        .iter()
        .map(|num_str| circom_str_to_lambda_field_element(num_str))
        .collect();
    adjust_lro_and_witness(&circom_r1cs, &mut l, &mut r, &mut o, &mut witness);

    // Lambdaworks considers "1" a public input, so compensate for it
    let num_of_pub_inputs = circom_r1cs["nPubInputs"].as_u64().unwrap() as usize + 1;

    (
        QAP::from_variable_matrices(num_of_pub_inputs, &l, &r, &o),
        witness,
    )
}

/// Takes as input circom.r1cs.json file and outputs LRO matrices
#[inline]
fn build_lro_from_circom_r1cs(circom_r1cs: &Value) -> [Vec<Vec<FrElement>>; 3] {
    let num_of_vars = circom_r1cs["nVars"].as_u64().unwrap() as usize; // Includes "1"
    let num_of_gates = circom_r1cs["nConstraints"].as_u64().unwrap() as usize;

    let mut l = vec![vec![FrElement::zero(); num_of_gates]; num_of_vars];
    let mut r = vec![vec![FrElement::zero(); num_of_gates]; num_of_vars];
    let mut o = vec![vec![FrElement::zero(); num_of_gates]; num_of_vars];

    for (constraint_idx, constraint) in circom_r1cs["constraints"]
        .as_array()
        .unwrap()
        .iter()
        .enumerate()
    {
        let constraint = constraint.as_array().unwrap();
        for (var_idx, str_val) in constraint[0].as_object().unwrap() {
            l[var_idx.parse::<usize>().unwrap()][constraint_idx] =
                circom_str_to_lambda_field_element(str_val.as_str().unwrap());
        }
        for (var_idx, str_val) in constraint[1].as_object().unwrap() {
            r[var_idx.parse::<usize>().unwrap()][constraint_idx] =
                circom_str_to_lambda_field_element(str_val.as_str().unwrap());
        }
        for (var_idx, str_val) in constraint[2].as_object().unwrap() {
            o[var_idx.parse::<usize>().unwrap()][constraint_idx] =
                circom_str_to_lambda_field_element(str_val.as_str().unwrap());
        }
    }

    [l, r, o]
}

/// Circom witness ordering: ["1", ..outputs, ...inputs, ...other_signals]
/// Lambda witness ordering: ["1", ...inputs, ..outputs,  ...other_signals]
/// Same applies to rows of LRO (each representing a variable)
/// This function compensates this difference
#[inline]
fn adjust_lro_and_witness(
    circom_r1cs: &Value,
    l: &mut [Vec<FrElement>],
    r: &mut [Vec<FrElement>],
    o: &mut [Vec<FrElement>],
    witness: &mut [FrElement],
) {
    let num_of_private_inputs = circom_r1cs["nPrvInputs"].as_u64().unwrap() as usize;
    let num_of_pub_inputs = circom_r1cs["nPubInputs"].as_u64().unwrap() as usize;
    let num_of_inputs = num_of_pub_inputs + num_of_private_inputs;
    let num_of_outputs = circom_r1cs["nOutputs"].as_u64().unwrap() as usize;

    let mut temp_l = Vec::with_capacity(num_of_inputs);
    let mut temp_r = Vec::with_capacity(num_of_inputs);
    let mut temp_o = Vec::with_capacity(num_of_inputs);
    let mut temp_witness = Vec::with_capacity(num_of_inputs);

    for i in 0..num_of_inputs {
        temp_l.push(l[num_of_outputs + 1 + i].clone());
        temp_r.push(r[num_of_outputs + 1 + i].clone());
        temp_o.push(o[num_of_outputs + 1 + i].clone());
        temp_witness.push(witness[num_of_outputs + 1 + i].clone());
    }

    for i in 0..num_of_inputs {
        let temp_l_i = l[1 + i].clone();
        l[1 + i].clone_from(&temp_l[i]);
        l[num_of_outputs + 1 + i].clone_from(&temp_l_i);

        let temp_r_i = r[1 + i].clone();
        r[1 + i].clone_from(&temp_r[i]);
        r[num_of_outputs + 1 + i].clone_from(&temp_r_i);

        let temp_o_i = o[1 + i].clone();
        o[1 + i].clone_from(&temp_o[i]);
        o[num_of_outputs + 1 + i].clone_from(&temp_o_i);

        let temp_witness_i = witness[1 + i].clone();
        witness[1 + i].clone_from(&temp_witness[i]);
        witness[num_of_outputs + 1 + i].clone_from(&temp_witness_i);
    }
}

#[inline]
fn circom_str_to_lambda_field_element(value: &str) -> FrElement {
    FrElement::from(&UnsignedInteger::<4>::from_dec_str(value).unwrap())
}
