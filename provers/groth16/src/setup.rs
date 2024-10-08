use crate::{common::*, QuadraticArithmeticProgram};
use lambdaworks_math::{
    cyclic_group::IsGroup,
    elliptic_curve::{
        short_weierstrass::{point::ShortWeierstrassProjectivePoint, traits::IsShortWeierstrass},
        traits::{IsEllipticCurve, IsPairing},
    },
    field::traits::IsField,
    traits::AsBytes,
};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};

/// Wrapper para G1Point que permite implementar Serialize
#[derive(Debug)]

pub struct SerializableG1Point(pub G1Point);

impl Serialize for SerializableG1Point {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        // Suponiendo que G1Point tiene un método `to_bytes()`
        let bytes = self.0.as_bytes();
        serializer.serialize_bytes(&bytes)
    }
}

/// Wrapper para G2Point que permite implementar Serialize
#[derive(Debug)]
pub struct SerializableG2Point(pub G2Point);

impl Serialize for SerializableG2Point {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let bytes = self.0.as_bytes();
        serializer.serialize_bytes(&bytes)
    }
}

/// Wrapper para PairingOutput que permite implementar Serialize
#[derive(Debug)]
pub struct SerializablePairingOutput(pub PairingOutput);

impl Serialize for SerializablePairingOutput {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        // Assuming that FieldElement implements AsBytes
        self.0.serialize(serializer)
    }
}

#[derive(Debug, Serialize)]
pub struct VerifyingKey {
    /// e([alpha]_1, [beta]_2) calculado durante el setup
    pub alpha_g1_times_beta_g2: SerializablePairingOutput,
    pub delta_g2: SerializableG2Point,
    pub gamma_g2: SerializableG2Point,
    /// [K_0(τ)]_1, [K_1(τ)]_1, ..., [K_k(τ)]_1
    /// donde K_i(τ) = γ^{-1} * (β*l(τ) + α*r(τ) + o(τ))
    /// y "k" es el número de entradas públicas
    pub verifier_k_tau_g1: Vec<SerializableG1Point>,
}
#[derive(Debug)] // needed?
pub struct ProvingKey {
    pub alpha_g1: G1Point,
    pub beta_g1: G1Point,
    pub beta_g2: G2Point,
    pub delta_g1: G1Point,
    pub delta_g2: G2Point,
    // [A_0(τ)]_1, [A_1(τ)]_1, ..., [A_n(τ)]_1
    pub l_tau_g1: Vec<G1Point>,
    // [B_0(τ)]_1, [B_1(τ)]_1, ..., [B_n(τ)]_1
    pub r_tau_g1: Vec<G1Point>,
    // [B_0(τ)]_2, [B_1(τ)]_2, ..., [B_n(τ)]_2
    pub r_tau_g2: Vec<G2Point>,
    // [K_{k+1}(τ)]_1, [K_{k+2}(τ)]_1, ..., [K_n(τ)]_1
    // where K_i(τ) = ƍ^{-1} * (β*l(τ) + α*r(τ) + o(τ))
    // and "k" is the number of public inputs
    pub prover_k_tau_g1: Vec<G1Point>,
    // [delta^{-1} * t(τ) * tau^0]_1, [delta^{-1} * t(τ) * τ^1]_1, ..., [delta^{-1} * t(τ) * τ^m]_1
    pub z_powers_of_tau_g1: Vec<G1Point>,
}

struct ToxicWaste {
    tau: FrElement,
    alpha: FrElement,
    beta: FrElement,
    gamma: FrElement,
    delta: FrElement,
}

impl ToxicWaste {
    pub fn new() -> Self {
        Self {
            tau: sample_fr_elem(),
            alpha: sample_fr_elem(),
            beta: sample_fr_elem(),
            gamma: sample_fr_elem(),
            delta: sample_fr_elem(),
        }
    }
}

pub fn setup(qap: &QuadraticArithmeticProgram) -> (ProvingKey, VerifyingKey) {
    let g1: G1Point = Curve::generator();
    let g2: G2Point = TwistedCurve::generator();

    let tw = ToxicWaste::new();

    let l_tau: Vec<_> = qap.l.iter().map(|p| p.evaluate(&tw.tau)).collect();
    let r_tau: Vec<_> = qap.r.iter().map(|p| p.evaluate(&tw.tau)).collect();

    let mut to_be_inversed = [tw.delta.clone(), tw.gamma.clone()];
    FrElement::inplace_batch_inverse(&mut to_be_inversed).unwrap();
    let [delta_inv, gamma_inv] = to_be_inversed;

    let k_tau: Vec<_> = l_tau
        .iter()
        .zip(&r_tau)
        .enumerate()
        .map(|(i, (l, r))| {
            let unshifted = &tw.beta * l + &tw.alpha * r + &qap.o[i].evaluate(&tw.tau);
            if i < qap.num_of_public_inputs {
                &gamma_inv * &unshifted
            } else {
                &delta_inv * &unshifted
            }
        })
        .collect();

    let alpha_g1 = g1.operate_with_self(tw.alpha.representative());
    let beta_g2 = g2.operate_with_self(tw.beta.representative());

    let alpha_g1_times_beta_g2 = Pairing::compute(&alpha_g1, &beta_g2).unwrap();

    let delta_g2 = g2.operate_with_self(tw.delta.representative());

    (
        ProvingKey {
            alpha_g1,
            beta_g1: g1.operate_with_self(tw.beta.representative()),
            beta_g2,
            delta_g1: g1.operate_with_self(tw.delta.representative()),
            delta_g2: delta_g2.clone(),
            l_tau_g1: batch_operate(&l_tau, &g1),
            r_tau_g1: batch_operate(&r_tau, &g1),
            r_tau_g2: batch_operate(&r_tau, &g2),
            prover_k_tau_g1: batch_operate(&k_tau[qap.num_of_public_inputs..], &g1),
            z_powers_of_tau_g1: batch_operate(
                &core::iter::successors(
                    // Start from delta^{-1} * t(τ)
                    // Note that t(τ) = (τ^N - 1) because our domain is roots of unity
                    Some(&delta_inv * (&tw.tau.pow(qap.num_of_gates) - FrElement::one())),
                    |prev| Some(prev * &tw.tau),
                )
                .take(qap.num_of_gates * 2)
                .collect::<Vec<_>>(),
                &g1,
            ),
        },
        VerifyingKey {
            alpha_g1_times_beta_g2: SerializablePairingOutput(alpha_g1_times_beta_g2),
            delta_g2: SerializableG2Point(delta_g2),
            gamma_g2: SerializableG2Point(g2.operate_with_self(tw.gamma.representative())),
            verifier_k_tau_g1: batch_operate(&k_tau[..qap.num_of_public_inputs], &g1)
                .into_iter()
                .map(SerializableG1Point)
                .collect(),
        },
    )
}

fn batch_operate<E: IsShortWeierstrass>(
    elems: &[FrElement],
    point: &ShortWeierstrassProjectivePoint<E>,
) -> Vec<ShortWeierstrassProjectivePoint<E>> {
    elems
        .iter()
        .map(|elem| point.operate_with_self(elem.representative()))
        .collect()
}
