use crate::{common::*, ProvingKey, QuadraticArithmeticProgram};
use lambdaworks_math::errors::DeserializationError;
use lambdaworks_math::traits::{AsBytes, Deserializable};
use lambdaworks_math::{cyclic_group::IsGroup, msm::pippenger::msm};
use std::mem::size_of;

#[derive(Debug)]
pub struct Proof {
    pub pi1: G1Point,
    pub pi2: G2Point,
    pub pi3: G1Point,
}

impl Proof {
    pub fn serialize(&self) -> Vec<u8> {
        let mut bytes: Vec<u8> = Vec::new();
        [
            Self::serialize_commitment(&self.pi1),
            Self::serialize_commitment(&self.pi2),
            Self::serialize_commitment(&self.pi3),
        ]
        .iter()
        .for_each(|serialized| {
            bytes.extend_from_slice(&(serialized.len() as u32).to_be_bytes());
            bytes.extend_from_slice(serialized);
        });
        bytes
    }

    pub fn deserialize(bytes: &[u8]) -> Result<Self, DeserializationError>
    where
        Self: Sized,
    {
        let (offset, pi1) = Self::deserialize_commitment::<G1Point>(bytes, 0)?;
        let (offset, pi2) = Self::deserialize_commitment::<G2Point>(bytes, offset)?;
        let (_, pi3) = Self::deserialize_commitment::<G1Point>(bytes, offset)?;
        Ok(Self { pi1, pi2, pi3 })
    }

    fn serialize_commitment<Commitment: AsBytes>(cm: &Commitment) -> Vec<u8> {
        cm.as_bytes()
    }

    // Repetitive. Same as in plonk/src/prover.rs
    fn deserialize_commitment<Commitment: Deserializable>(
        bytes: &[u8],
        offset: usize,
    ) -> Result<(usize, Commitment), DeserializationError> {
        let mut offset = offset;
        let element_size_bytes: [u8; size_of::<u32>()] = bytes
            .get(offset..offset + size_of::<u32>())
            .ok_or(DeserializationError::InvalidAmountOfBytes)?
            .try_into()
            .map_err(|_| DeserializationError::InvalidAmountOfBytes)?;
        let element_size = u32::from_be_bytes(element_size_bytes) as usize;
        offset += size_of::<u32>();
        let commitment = Commitment::deserialize(
            bytes
                .get(offset..offset + element_size)
                .ok_or(DeserializationError::InvalidAmountOfBytes)?,
        )?;
        offset += element_size;
        Ok((offset, commitment))
    }
}

pub struct Prover;
impl Prover {
    pub fn prove(w: &[FrElement], qap: &QuadraticArithmeticProgram, pk: &ProvingKey) -> Proof {
        let h_coefficients = qap
            .calculate_h_coefficients(w)
            .iter()
            .map(|elem| elem.representative())
            .collect::<Vec<_>>();

        let w = w
            .iter()
            .map(|elem| elem.representative())
            .collect::<Vec<_>>();

        // Sample randomness for hiding
        let r = sample_fr_elem();
        let s = sample_fr_elem();

        // [π_1]_1
        let pi1 = msm(&w, &pk.l_tau_g1)
            .unwrap()
            .operate_with(&pk.alpha_g1)
            .operate_with(&pk.delta_g1.operate_with_self(r.representative()));

        // [π_2]_2
        let pi2 = msm(&w, &pk.r_tau_g2)
            .unwrap()
            .operate_with(&pk.beta_g2)
            .operate_with(&pk.delta_g2.operate_with_self(s.representative()));

        // [ƍ^{-1} * t(τ)*h(τ)]_1
        let t_tau_h_tau_assigned_g1 = msm(
            &h_coefficients,
            &pk.z_powers_of_tau_g1[..h_coefficients.len()],
        )
        .unwrap();

        // [ƍ^{-1} * (β*l(τ) + α*r(τ) + o(τ))]_1
        let k_tau_assigned_prover_g1 = msm(
            &w[qap.num_of_public_inputs..],
            &pk.prover_k_tau_g1[..qap.num_of_private_inputs()],
        )
        .unwrap();

        // [π_2]_1
        let pi2_g1 = msm(&w, &pk.r_tau_g1)
            .unwrap()
            .operate_with(&pk.beta_g1)
            .operate_with(&pk.delta_g1.operate_with_self(s.representative()));

        // [π_3]_1
        let pi3 = k_tau_assigned_prover_g1
            .operate_with(&t_tau_h_tau_assigned_g1)
            // s[π_1]_1
            .operate_with(&pi1.operate_with_self(s.representative()))
            // r[π_2]_1
            .operate_with(&pi2_g1.operate_with_self(r.representative()))
            // -rs[ƍ]_1
            .operate_with(&pk.delta_g1.operate_with_self((-(&r * &s)).representative()));

        Proof { pi1, pi2, pi3 }
    }
}

#[cfg(test)]
mod tests {
    use lambdaworks_math::elliptic_curve::traits::IsEllipticCurve;

    use super::*;

    #[test]
    fn serde() {
        let proof = Proof {
            pi1: Curve::generator().operate_with_self(sample_fr_elem().representative()),
            pi2: TwistedCurve::generator().operate_with_self(sample_fr_elem().representative()),
            pi3: Curve::generator().operate_with_self(sample_fr_elem().representative()),
        };
        let deserialized_proof = Proof::deserialize(&proof.serialize()).unwrap();

        assert_eq!(proof.pi1, deserialized_proof.pi1);
        assert_eq!(proof.pi2, deserialized_proof.pi2);
        assert_eq!(proof.pi3, deserialized_proof.pi3);
    }
}
