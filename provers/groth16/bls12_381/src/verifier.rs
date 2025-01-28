use lambdaworks_math::{elliptic_curve::traits::IsPairing, msm::pippenger::msm};

use crate::common::{FrElement, Pairing};
use crate::prover::Proof;
use crate::setup::VerifyingKey;

pub fn verify(vk: &VerifyingKey, proof: &Proof, pub_inputs: &[FrElement]) -> bool {
    // [γ^{-1} * (β*l(τ) + α*r(τ) + o(τ))]_1
    let k_tau_assigned_verifier_g1 = msm(
        &pub_inputs
            .iter()
            .map(|elem| elem.representative())
            .collect::<Vec<_>>(),
        &vk.verifier_k_tau_g1,
    )
    .unwrap();

    Pairing::compute(&proof.pi3, &vk.delta_g2).unwrap()
        * vk.alpha_g1_times_beta_g2.clone()
        * Pairing::compute(&k_tau_assigned_verifier_g1, &vk.gamma_g2).unwrap()
        == Pairing::compute(&proof.pi1, &proof.pi2).unwrap()
}
