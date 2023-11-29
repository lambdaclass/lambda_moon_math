use lambdaworks_groth16::{common::*, setup, test_circuits::*, verify, Proof, Prover};

// #[test]
// fn vitalik() {
//     let qap = vitalik_qap(); // x^3 + x + 5 = 35
//     let (pk, vk) = setup(&qap);

//     for w in [
//         ["0x1", "0x3", "0x23", "0x9", "0x1b", "0x1e"], // x = 3
//         ["0x1", "0x1", "0x7", "0x1", "0x1", "0x2"],    // x = 1
//     ] {
//         let w = w.map(FrElement::from_hex_unchecked).to_vec();
//         let accept = verify(
//             &vk,
//             &Prover::prove(&w, &qap, &pk),
//             &w[..qap.num_of_public_inputs],
//         );
//         assert!(accept);
//     }
// }

// #[test]
// fn vitalik_from_csv() {
//     let qap = vitalik_qap_from_csv(); // x^3 + x + 5 = 35
//     let (pk, vk) = setup(&qap);

//     for w in [
//         ["0x1", "0x3", "0x23", "0x9", "0x1b", "0x1e"], // x = 3
//         ["0x1", "0x1", "0x7", "0x1", "0x1", "0x2"],    // x = 1
//     ] {
//         let w = w.map(FrElement::from_hex_unchecked).to_vec();
//         let accept = verify(
//             &vk,
//             &Prover::prove(&w, &qap, &pk),
//             &w[..qap.num_of_public_inputs],
//         );
//         assert!(accept);
//     }
// }

#[test]
fn rsa_from_csv() {
    let qap = rsa_qap_from_csv();
    let (pk, vk) = setup(&qap);

    for w in [[
        "0x1",
        "0x756BA64",
        "0x933BFB764",
        "0x6f7739dd2e43813d06d687794654a135732ee2d6029cd727ed3661b36187d682", // inv(p-1)
        "0x16d2c9bdf8fe34aa36e29687fb7f10210c7c8916fd9946e96733333300cccccd", // inv(q-1)
        "0x141",
        "0x933BFB764",
    ]] {
        let w = w.map(FrElement::from_hex_unchecked).to_vec();
        let accept = verify(
            &vk,
            &Prover::prove(&w, &qap, &pk),
            &w[..qap.num_of_public_inputs],
        );
        assert!(accept);
    }
}

#[test]
fn qap_2() {
    let qap = test_qap_2();
    let (pk, vk) = setup(&qap);

    // 1, x, y, ~out, sym_1, sym_2, sym_3, sym_4
    let w = ["0x1", "0x5", "0x3", "0x0", "0x19", "0x9", "0x0", "0x0"] // x = 3
        .map(FrElement::from_hex_unchecked)
        .to_vec();

    let serialized_proof = Prover::prove(&w, &qap, &pk).serialize();
    let deserialized_proof = Proof::deserialize(&serialized_proof).unwrap();

    let accept = verify(&vk, &deserialized_proof, &w[..qap.num_of_public_inputs]);
    assert!(accept);
}
