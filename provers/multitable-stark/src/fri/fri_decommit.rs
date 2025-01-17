use lambdaworks_crypto::merkle_tree::proof::Proof;
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::traits::IsField;

use crate::config::Commitment;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FriDecommitment<F: IsField> {
    pub layers_auth_paths: Vec<Proof<Commitment>>,
    pub layers_evaluations_sym: Vec<FieldElement<F>>,
}
