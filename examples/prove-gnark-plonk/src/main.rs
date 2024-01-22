
use std::fs;

use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::traits::IsFFTField;
use lambdaworks_math::{
    elliptic_curve::short_weierstrass::curves::bls12_381::default_types::{FrElement, FrField},
    polynomial::Polynomial,
};

use lambdaworks_plonk::test_utils::circuit_json::common_preprocessed_input_from_json;

use serde::{Deserialize, Serialize};

fn main() {
        // This is the circuit for x * e + 5 == y
        let data = fs::read_to_string("go_exporter/example/frontend_precomputed_values.json").expect("Unable to read file");

        let srs = common_preprocessed_input_from_json(&data);

        // Public input
        let x = FieldElement::from(2_u64);
        let y = FieldElement::from(11_u64);

        // Private variable
        let e = FieldElement::from(3_u64);

        let public_input = vec![x.clone(), y];
        let witness = test_witness_2(x, e);

        let kzg = KZG::new(srs);
        let verifying_key = setup(&common_preprocessed_input, &kzg);
        let random_generator = TestRandomFieldGenerator {};

        let prover = Prover::new(kzg.clone(), random_generator);
        let proof = prover.prove(
            &witness,
            &public_input,
            &common_preprocessed_input,
            &verifying_key,
        );

        let verifier = Verifier::new(kzg);
        assert!(verifier.verify(
            &proof,
            &public_input,
            &common_preprocessed_input,
            &verifying_key
        ));
    }
    fn test_happy_path_from_json() {
        let (witness, common_preprocessed_input, public_input) =
            common_preprocessed_input_from_json(
                r#"{
            "N": 4,
            "N_Padded": 4,
            "Omega": "8d51ccce760304d0ec030002760300000001000000000000",
             "Input": [
             "2",
             "4"
            ],
            "Ql": [
             "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000000",
             "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000000",
             "0",
             "1"
            ],
            "Qr": [
             "0",
             "0",
             "0",
             "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000000"
            ],
            "Qm": [
             "0",
             "0",
             "1",
             "0"
            ],
            "Qo": [
             "0",
             "0",
             "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000000",
             "0"
            ],
            "Qc": [
             "0",
             "0",
             "0",
             "0"
            ],
            "A": [
             "2",
             "4",
             "2",
             "4"
            ],
            "B": [
             "2",
             "2",
             "2",
             "4"
            ],
            "C": [
             "2",
             "2",
             "4",
             "2"
            ],
            "Permutation": [
             11,
             3,
             2,
             1,
             0,
             4,
             5,
             10,
             6,
             8,
             7,
             9
            ]
           }"#,
            );
        let srs = test_srs(common_preprocessed_input.n);

        let kzg = KZG::new(srs);
        let verifying_key = setup(&common_preprocessed_input, &kzg);
        let random_generator = TestRandomFieldGenerator {};

        let prover = Prover::new(kzg.clone(), random_generator);
        let proof = prover.prove(
            &witness,
            &public_input,
            &common_preprocessed_input,
            &verifying_key,
        );

        let verifier = Verifier::new(kzg);
        assert!(verifier.verify(
            &proof,
            &public_input,
            &common_preprocessed_input,
            &verifying_key
        ));
}
