use lambdaworks_circom_adapter::{circom_to_lambda, read_circom_r1cs, read_circom_witness};
use lambdaworks_groth16::{common::FrElement, QuadraticArithmeticProgram};

/// Converts following Circom circuit and inputs into Lambdaworks-compatible QAP and witness assignments.
///
/// ```csharp
/// template Test() {
///	signal input x;
///	signal output out;
///
///	signal sym_1;
///	signal y;
///
///	sym_1 <== x * x;
///	out <== (sym_1 * x) + (x + 5);
/// }
///
/// // { "x": 3 }
/// ```
///
#[test]
fn vitalik_w_and_qap() {
    let circom_wtns = read_circom_witness("./tests/vitalik_example/witness.json")
        .expect("could not read witness");
    let circom_r1cs =
        read_circom_r1cs("./tests/vitalik_example/test.r1cs.json").expect("could not read r1cs");

    let (qap, wtns, pubs) = circom_to_lambda(circom_r1cs, circom_wtns);

    // Circom witness contains outputs before circuit inputs where Lambdaworks puts inputs before the output.
    // Freshly generated witness assignment "w" must be in form ["1", "x", "~out", "sym_1"]
    assert_eq!(
        wtns,
        // ["1", "3", "23", "9"]
        // 1, out, x, sym_1
        ["1", "23", "3", "9"]
            .map(FrElement::from_hex_unchecked)
            .to_vec(),
        "incorrect witness"
    );
    assert_eq!(
        pubs,
        ["1", "23"].map(FrElement::from_hex_unchecked).to_vec(),
        "incorrect public signals"
    );

    // Regarding QAP, we expect 2 constraints in following form
    // -sym_1 = -x * x
    // x + 5 - ~out = -sym_1 * x
    //
    // Same ordering difference exists for variable matrices, too. Circom adapter changes the
    // order of the rows in the same way it rearranges the witness ordering.

    const _M1_: &str = "0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000000";
    const _0_: &str = "0x0";
    const _1_: &str = "0x1";

    #[rustfmt::skip]
    let [temp_l, temp_r, temp_o] = [
        // L //
        [
            [_0_,  _0_],  // 1
            [_0_,  _0_],  // ~out
            [_M1_, _0_],  // x
            [_0_,  _M1_], // sym_1
        ],
        // R //
        [
            [_0_, _0_],   // 1
            [_0_, _0_],   // ~out
            [_1_, _1_],   // x
            [_0_, _0_],   // sym_1
        ],
        // O //
        [
            [_0_, "5"],   // 1
            [_0_, _M1_],  // ~out
            [_0_, _1_],   // x
            [_M1_, _0_],  // sym_1
        ],
    ]
    .map(|matrix| matrix.map(|row| row.map(FrElement::from_hex_unchecked).to_vec()));

    let expected_qap =
        QuadraticArithmeticProgram::from_variable_matrices(pubs.len(), &temp_l, &temp_r, &temp_o);

    assert_eq!(qap.l, expected_qap.l);
    assert_eq!(qap.r, expected_qap.r);
    assert_eq!(qap.o, expected_qap.o);

    // check proofs
    let (pk, vk) = lambdaworks_groth16::setup(&qap);
    let proof = lambdaworks_groth16::Prover::prove(&wtns, &qap, &pk);
    let accept = lambdaworks_groth16::verify(&vk, &proof, &pubs);
    assert!(accept, "proof verification failed");
}
