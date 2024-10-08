use lambdaworks_math::{elliptic_curve::traits::IsPairing, msm::pippenger::msm};

use crate::common::{FrElement, G1Point, Pairing};
use crate::prover::Proof;
use crate::setup::VerifyingKey;

pub fn verify(vk: &VerifyingKey, proof: &Proof, pub_inputs: &[FrElement]) -> bool {
    // Extraer los G1Points originales de los wrappers y clonarlos
    let verifier_k_tau_g1: Vec<G1Point> = vk
        .verifier_k_tau_g1
        .iter()
        .map(|wrapper| wrapper.0.clone())
        .collect();

    // Realizar la MSM pasando el slice de G1Point
    let k_tau_assigned_verifier_g1 = msm(
        &pub_inputs
            .iter()
            .map(|elem| elem.representative())
            .collect::<Vec<_>>(),
        &verifier_k_tau_g1, // Aquí pasamos &[G1Point]
    )
    .unwrap();

    // Extraer los G2Points y PairingOutput originales
    let delta_g2 = &vk.delta_g2.0;
    let gamma_g2 = &vk.gamma_g2.0;
    let alpha_g1_times_beta_g2 = vk.alpha_g1_times_beta_g2.0.clone();

    // Calcular los emparejamientos
    let pairing_pi3_delta_g2 = Pairing::compute(&proof.pi3, delta_g2).unwrap();
    let pairing_k_tau_gamma_g2 = Pairing::compute(&k_tau_assigned_verifier_g1, gamma_g2).unwrap();
    let pairing_pi1_pi2 = Pairing::compute(&proof.pi1, &proof.pi2).unwrap();

    // Verificar la ecuación final
    pairing_pi3_delta_g2 * alpha_g1_times_beta_g2 * pairing_k_tau_gamma_g2 == pairing_pi1_pi2
}
