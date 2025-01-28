use crate::arkworks_cs_to_lambda_cs;
use ark_bls12_381::Fr;
use ark_relations::{lc, r1cs::ConstraintSystem, r1cs::Variable};
use lambdaworks_groth16_bls12_381::{setup, verify, Prover, QuadraticArithmeticProgram};
use rand::Rng;

#[test]
fn pinocchio_paper_example() {
    /*
        pub inp a, b, c, d
        pub out result
        sig e

        c * d = e
        (a + b) + e = result
    */
    let cs = ConstraintSystem::<Fr>::new_ref();

    let _a = Fr::from(1555);
    let _b = Fr::from(25555);
    let _c = Fr::from(35555);
    let _d = Fr::from(45555);

    let a = cs.new_input_variable(|| Ok(_a)).unwrap();
    let b = cs.new_input_variable(|| Ok(_b)).unwrap();
    let c = cs.new_input_variable(|| Ok(_c)).unwrap();
    let d = cs.new_input_variable(|| Ok(_d)).unwrap();

    let e = cs.new_witness_variable(|| Ok(_c * _d)).unwrap();
    cs.enforce_constraint(lc!() + c, lc!() + d, lc!() + e)
        .unwrap();

    let result = cs.new_input_variable(|| Ok(_c * _d * (_a + _b))).unwrap();

    cs.enforce_constraint(lc!() + a + b, lc!() + e, lc!() + result)
        .unwrap();

    let lambda_cs = arkworks_cs_to_lambda_cs(&cs);

    let qap = QuadraticArithmeticProgram::from_r1cs(lambda_cs.constraints);

    let (pk, vk) = setup(&qap);

    let accept = verify(
        &vk,
        &Prover::prove(&lambda_cs.witness, &qap, &pk),
        &lambda_cs.witness[..qap.num_of_public_inputs],
    );
    assert!(accept);
}

#[test]
fn vitalik_example() {
    /*
        pub out ~out
        sig x, sym_1, y, sym_2

        x * x = sym_1;
        sym_1 * x = y;
        (y + x) * 1 = sym_2
        (sym_2 + 5) * 1 = ~out
    */
    let cs = ConstraintSystem::<Fr>::new_ref();

    //["0x1", "0x3", "0x23", "0x9", "0x1b", "0x1e"]
    let _x = Fr::from(3);
    let _sym_1 = Fr::from(9);
    let _y = Fr::from(27);
    let _sym_2 = Fr::from(30);

    let _out = Fr::from(35);

    let x = cs.new_witness_variable(|| Ok(_x)).unwrap();
    let sym_1 = cs.new_witness_variable(|| Ok(_sym_1)).unwrap();
    let y = cs.new_witness_variable(|| Ok(_y)).unwrap();
    let sym_2 = cs.new_witness_variable(|| Ok(_sym_2)).unwrap();

    let out = cs.new_input_variable(|| Ok(_out)).unwrap();

    cs.enforce_constraint(lc!() + x, lc!() + x, lc!() + sym_1)
        .unwrap();
    cs.enforce_constraint(lc!() + sym_1, lc!() + x, lc!() + y)
        .unwrap();
    cs.enforce_constraint(lc!() + y + x, lc!() + Variable::One, lc!() + sym_2)
        .unwrap();
    cs.enforce_constraint(
        lc!() + sym_2 + (Fr::from(5), Variable::One),
        lc!() + Variable::One,
        lc!() + out,
    )
    .unwrap();

    let lambda_cs = arkworks_cs_to_lambda_cs(&cs);

    let qap = QuadraticArithmeticProgram::from_r1cs(lambda_cs.constraints);

    let (pk, vk) = setup(&qap);

    let accept = verify(
        &vk,
        &Prover::prove(&lambda_cs.witness, &qap, &pk),
        &lambda_cs.witness[..qap.num_of_public_inputs],
    );
    assert!(accept);
}

#[test]
fn failing_vitalik() {
    // Same circuit as vitalik_example, but with an incorrect witness assignment.
    let cs = ConstraintSystem::<Fr>::new_ref();

    //["0x1", "0x3", "0x23", "0x9", "0x1b", "0x1e"]
    let _x = Fr::from(3);
    let _sym_1 = Fr::from(10); // should be have been 9
    let _y = Fr::from(27);
    let _sym_2 = Fr::from(30);

    let _out = Fr::from(35);

    let x = cs.new_witness_variable(|| Ok(_x)).unwrap();
    let sym_1 = cs.new_witness_variable(|| Ok(_sym_1)).unwrap();
    let y = cs.new_witness_variable(|| Ok(_y)).unwrap();
    let sym_2 = cs.new_witness_variable(|| Ok(_sym_2)).unwrap();

    let out = cs.new_input_variable(|| Ok(_out)).unwrap();

    cs.enforce_constraint(lc!() + x, lc!() + x, lc!() + sym_1)
        .unwrap();
    cs.enforce_constraint(lc!() + sym_1, lc!() + x, lc!() + y)
        .unwrap();
    cs.enforce_constraint(lc!() + y + x, lc!() + Variable::One, lc!() + sym_2)
        .unwrap();
    cs.enforce_constraint(
        lc!() + sym_2 + (Fr::from(5), Variable::One),
        lc!() + Variable::One,
        lc!() + out,
    )
    .unwrap();

    let lambda_cs = arkworks_cs_to_lambda_cs(&cs);

    let qap = QuadraticArithmeticProgram::from_r1cs(lambda_cs.constraints);

    let (pk, vk) = setup(&qap);

    let accept = verify(
        &vk,
        &Prover::prove(&lambda_cs.witness, &qap, &pk),
        &lambda_cs.witness[..qap.num_of_public_inputs],
    );
    assert!(!accept);
}

#[test]
fn exponentiation_example() {
    /*
        Generates a "linear exponentiation" circuit with a random base and a random exponent.
        Only the output ~out is public input.
    */
    let cs = ConstraintSystem::<Fr>::new_ref();

    let mut rng = rand::thread_rng();
    let x = rng.gen::<u64>();
    let exp = rng.gen::<u8>(); // Bigger data types take too much time for a test

    let x = Fr::from(x);
    let mut _x = cs.new_witness_variable(|| Ok(x)).unwrap();

    let mut acc = Fr::from(x);
    let mut _acc = cs.new_witness_variable(|| Ok(x)).unwrap();

    for _ in 0..exp - 1 {
        acc *= x;
        let _new_acc = cs.new_witness_variable(|| Ok(acc)).unwrap();
        cs.enforce_constraint(lc!() + _acc, lc!() + _x, lc!() + _new_acc)
            .unwrap();
        _acc = _new_acc;
    }

    let _out = cs.new_input_variable(|| Ok(acc)).unwrap();
    cs.enforce_constraint(lc!() + _out, lc!() + Variable::One, lc!() + _acc)
        .unwrap();

    let lambda_cs = arkworks_cs_to_lambda_cs(&cs);

    let qap = QuadraticArithmeticProgram::from_r1cs(lambda_cs.constraints);

    let (pk, vk) = setup(&qap);

    let proof = Prover::prove(&lambda_cs.witness, &qap, &pk);

    let public_inputs = &lambda_cs.witness[..qap.num_of_public_inputs];
    let accept = verify(&vk, &proof, public_inputs);
    assert!(accept);
}
