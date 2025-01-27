// allowing unused mut until we solve the signal ordering thing
#![allow(unused_mut)]

use lambdaworks_groth16::{common::FrElement, QuadraticArithmeticProgram};

mod readers;
pub use readers::*;

/// Given a Circom R1CS and witness it returns a QAP, witness, and public signals; all compatible with Lambdaworks.
pub fn circom_to_lambda(
    circom_r1cs: CircomR1CS,
    mut witness: CircomWitness,
) -> (QuadraticArithmeticProgram, Vec<FrElement>, Vec<FrElement>) {
    let num_of_outputs = circom_r1cs.num_outputs;
    // let num_of_private_inputs = circom_r1cs.num_priv_inputs;
    let num_of_pub_inputs = circom_r1cs.num_pub_inputs;

    let [mut l, mut r, mut o] = build_lro_from_circom_r1cs(circom_r1cs);
    // adjust_lro_and_witness(
    //     num_of_outputs,
    //     num_of_private_inputs,
    //     num_of_pub_inputs,
    //     &mut l,
    //     &mut r,
    //     &mut o,
    //     &mut witness,
    // );

    // we could get a slice using the QAP but the QAP does not keep track of the number of private inputs;
    // so instead we get the public signals here
    let mut public_inputs = witness[..num_of_pub_inputs + num_of_outputs + 1].to_vec();
    // public_inputs.insert(0, witness[0].clone()); // this is usually 1

    (
        // Lambdaworks considers "1" a public input, so compensate for it
        QuadraticArithmeticProgram::from_variable_matrices(public_inputs.len(), &l, &r, &o),
        witness,
        public_inputs,
    )
}

/// Takes as input circom.r1cs.json file and outputs LRO matrices
#[inline]
fn build_lro_from_circom_r1cs(circom_r1cs: CircomR1CS) -> [Vec<Vec<FrElement>>; 3] {
    let mut l = vec![vec![FrElement::zero(); circom_r1cs.num_constraints]; circom_r1cs.num_vars];
    let mut r = l.clone(); // same as above
    let mut o = l.clone(); // same as above

    // assign each constraint from the R1CS hash-maps to LRO matrices
    for (constraint_idx, constraint) in circom_r1cs.constraints.into_iter().enumerate() {
        // destructuring here to avoid clones
        let [lc, rc, oc] = constraint;

        for (var_idx, val) in lc {
            l[var_idx][constraint_idx] = val;
        }
        for (var_idx, val) in rc {
            r[var_idx][constraint_idx] = val;
        }
        for (var_idx, val) in oc {
            o[var_idx][constraint_idx] = val;
        }
    }

    [l, r, o]
}

/// Change the ordering of private-inputs and public-outputs from Circom to Lambdaworks style,
/// for both LRO matrices and witness.
#[inline]
#[deprecated = "seems to be working without this?"]
#[allow(unused)]
fn adjust_lro_and_witness(
    num_of_outputs: usize,
    num_of_private_inputs: usize,
    num_of_pub_inputs: usize,
    l: &mut [Vec<FrElement>],
    r: &mut [Vec<FrElement>],
    o: &mut [Vec<FrElement>],
    witness: &mut [FrElement],
) {
    let num_of_inputs = num_of_pub_inputs + num_of_private_inputs;

    let mut temp_l = Vec::with_capacity(num_of_inputs);
    let mut temp_r = Vec::with_capacity(num_of_inputs);
    let mut temp_o = Vec::with_capacity(num_of_inputs);
    let mut temp_witness = Vec::with_capacity(num_of_inputs);

    temp_l.extend_from_slice(&l[num_of_outputs + 1..num_of_outputs + 1 + num_of_inputs]);
    temp_r.extend_from_slice(&r[num_of_outputs + 1..num_of_outputs + 1 + num_of_inputs]);
    temp_o.extend_from_slice(&o[num_of_outputs + 1..num_of_outputs + 1 + num_of_inputs]);
    temp_witness
        .extend_from_slice(&witness[num_of_outputs + 1..num_of_outputs + 1 + num_of_inputs]);

    let temp_l_i = &l[1..=num_of_inputs].to_vec();
    l[1..=num_of_inputs].clone_from_slice(&temp_l[..num_of_inputs]);
    l[num_of_outputs + 1..num_of_inputs + num_of_outputs + 1].clone_from_slice(&temp_l_i);

    let temp_r_i = &r[1..=num_of_inputs].to_vec();
    r[1..=num_of_inputs].clone_from_slice(&temp_r[..num_of_inputs]);
    r[num_of_outputs + 1..num_of_inputs + num_of_outputs + 1].clone_from_slice(&temp_r_i);

    let temp_o_i = &o[1..=num_of_inputs].to_vec();
    o[1..=num_of_inputs].clone_from_slice(&temp_o[..num_of_inputs]);
    o[num_of_outputs + 1..num_of_inputs + num_of_outputs + 1].clone_from_slice(&temp_o_i);

    let temp_witness_i = &witness[1..=num_of_inputs].to_vec();
    witness[1..=num_of_inputs].clone_from_slice(&temp_witness[..num_of_inputs]);
    witness[num_of_outputs + 1..num_of_inputs + num_of_outputs + 1]
        .clone_from_slice(&temp_witness_i);
}

// #[inline]
// fn circom_str_to_lambda_field_element(value: &str) -> FrElement {
//     FrElement::from(&UnsignedInteger::<4>::from_dec_str(value).unwrap())
// }
