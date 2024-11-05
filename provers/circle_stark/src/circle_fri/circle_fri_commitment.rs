use super::circle_fri_functions::{fold_circle_to_univariate, fold_univariate};
use crate::config::{BatchedMerkleTree, BatchedMerkleTreeBackend};
use lambdaworks_crypto::fiat_shamir::is_transcript::IsTranscript;
use lambdaworks_crypto::merkle_tree::{merkle::MerkleTree, traits::IsMerkleTreeBackend};
use lambdaworks_math::traits::AsBytes;
use lambdaworks_math::{
    circle::cosets::Coset,
    field::{
        element::FieldElement,
        traits::{IsFFTField, IsField},
    },
    polynomial::Polynomial,
};

/// Represents a layer in the Circle FRI protocol.
#[derive(Clone)]
pub struct CircleFriLayer<F, B>
where
    F: IsField,
    FieldElement<F>: AsBytes,
    B: IsMerkleTreeBackend,
{
    pub evaluation: Vec<FieldElement<F>>,
    pub merkle_tree: MerkleTree<B>,
    pub coset_shift: FieldElement<F>,
    pub domain_size: usize,
}

impl<F, B> CircleFriLayer<F, B>
where
    F: IsField,
    FieldElement<F>: AsBytes,
    B: IsMerkleTreeBackend,
{
    pub fn new(
        evaluation: &[FieldElement<F>],
        merkle_tree: MerkleTree<B>,
        coset_shift: FieldElement<F>,
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
