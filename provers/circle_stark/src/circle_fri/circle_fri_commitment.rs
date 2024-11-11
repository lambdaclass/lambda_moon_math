use super::circle_fri_functions::{fold_circle_to_univariate, fold_univariate};
use crate::config::{BatchedMerkleTree, BatchedMerkleTreeBackend};
use lambdaworks_crypto::fiat_shamir::is_transcript::IsTranscript;
use lambdaworks_crypto::merkle_tree::{merkle::MerkleTree, traits::IsMerkleTreeBackend};
use lambdaworks_math::field::fields::mersenne31::field::Mersenne31Field;
use lambdaworks_math::traits::AsBytes;
use lambdaworks_math::{
    circle::cosets::Coset,
    field::{
        element::FieldElement,
        traits::{IsFFTField, IsField},
    },
    polynomial::Polynomial,
};
type F = Mersenne31Field;
type FE = FieldElement<F>;

/// Represents a layer in the Circle FRI protocol.
#[derive(Clone)]
pub struct CircleFriLayer {
    pub evaluation: Vec<FE>,
    pub merkle_tree: MerkleTree<BatchedMerkleTreeBackend<F>>,
    pub coset_shift: FE,
    pub domain_size: usize,
}

impl CircleFriLayer {
    pub fn new(
        evaluation: &[FE],
        merkle_tree: MerkleTree<BatchedMerkleTreeBackend<F>>,
        coset_shift: FE,
        domain_size: usize,
    ) -> Self {
        Self {
            evaluation: evaluation.to_vec(),
            merkle_tree,
            coset_shift,
            domain_size,
        }
    }
}
