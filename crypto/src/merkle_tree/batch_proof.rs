use super::{
    traits::IsMerkleTreeBackend,
    utils::{get_parent_pos, get_sibling_pos},
};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
extern crate std;
use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BatchProof<T>
where
    T: PartialEq + Eq,
{
    pub auth: HashMap<usize, T>,
}

#[cfg(feature = "serde")]
impl<T> BatchProof<T>
where
    T: PartialEq + Eq + Serialize + for<'a> Deserialize<'a>,
{
    pub fn serialize(&self) -> Vec<u8> {
        bincode::serde::encode_to_vec(self, bincode::config::standard()).unwrap()
    }

    pub fn deserialize(bytes: &[u8]) -> Self {
        let (decoded, _): (BatchProof<T>, usize) =
            bincode::serde::decode_from_slice(bytes, bincode::config::standard()).unwrap();
        decoded
    }
}

impl<T> BatchProof<T>
where
    T: PartialEq + Eq,
{
    pub fn verify<B>(&self, root_hash: B::Node, hashed_leaves: HashMap<usize, B::Node>) -> bool
    where
        B: IsMerkleTreeBackend<Node = T>,
    {
        let root_pos = 0;

        // Iterate the levels starting from the leaves, and build the upper level only using
        // the provided leaves and the auth map.
        // Return true if the constructed root matches the given one.
        let mut current_level = hashed_leaves;
        loop {
            let mut parent_level = HashMap::<usize, B::Node>::new();

            for (pos, node) in current_level.iter() {
                // Levels are expected to have tuples of nodes. If the first one was
                // already processed and parent was set, skip the sibling.
                let parent_pos = get_parent_pos(*pos);
                if parent_level.contains_key(&parent_pos) {
                    continue;
                }

                // Get the sibling node from the current level.
                // If doesn't exist, then it must have been provided in the batch auth.
                // If neither, then verification fails.
                let sibling_pos = get_sibling_pos(*pos);
                let sibling_node = if current_level.contains_key(&sibling_pos) {
                    current_level.get(&sibling_pos).unwrap()
                } else if self.auth.contains_key(&sibling_pos) {
                    self.auth.get(&sibling_pos).unwrap()
                } else {
                    panic!("Leaf with position {pos} has sibling {sibling_pos}, but it's not included in the auth map. ");
                };

                let parent_node = B::hash_new_parent(node, sibling_node);

                // Root must match the provided root hash.
                if parent_pos == root_pos {
                    return parent_node == root_hash;
                }

                // Create a new element for the next, upper level
                parent_level.insert(parent_pos, parent_node);
            }

            // We didn't create any parents, and we didn't reach the root neither. Verification fails.
            if parent_level.is_empty() {
                return false;
            }

            // Issue the next level in the next iteration
            current_level = parent_level;
        }
    }
}

#[cfg(test)]
mod tests {

    // #[cfg(feature = "alloc")]
    use crate::merkle_tree::{merkle::MerkleTree, test_merkle::TestBackend as TB};
    use alloc::vec::Vec;
    use lambdaworks_math::field::{element::*, fields::u64_prime_field::U64PrimeField};
    extern crate std;
    use std::collections::HashMap;

    /// Goldilocks
    pub type Ecgfp5 = U64PrimeField<0xFFFF_FFFF_0000_0001_u64>;
    pub type Ecgfp5FE = FieldElement<Ecgfp5>;
    pub type TestBackend = TB<Ecgfp5>;
    pub type TestMerkleTreeEcgfp = MerkleTree<TestBackend>;

    //
    //          20
    //       /      \
    //      6       14
    //     / \      / \
    //    2   4    6   8
    //
    // Proves inclusion of leaves whose indices are passed into 'proven_leaf_indices'
    // array. These leaf indices start from 0 for the first leaf, 2 in the example above.
    // If leaf_indices is [0, 3], then the test will create proof and verify inclusion
    // of leaves with indices 0 and 3, that are, 2 and 8.
    #[test]
    fn batch_proof_pen_and_paper_example() {
        let leaves_values: Vec<Ecgfp5FE> = (1..u64::pow(2, 2) + 1).map(Ecgfp5FE::new).collect();
        let merkle_tree: MerkleTree<TestBackend> = TestMerkleTreeEcgfp::build(&leaves_values);

        let proven_leaves_indices = [0, 3];
        let first_leaf_pos = merkle_tree.nodes_len() / 2;
        let proven_leaves_positions: Vec<usize> = proven_leaves_indices
            .iter()
            .map(|leaf_index| leaf_index + first_leaf_pos)
            .collect();

        let batch_proof = merkle_tree
            .get_batch_proof(&proven_leaves_positions)
            .unwrap();

        let proven_leaves_values_hashed: HashMap<usize, Ecgfp5FE> = proven_leaves_positions
            .iter()
            .map(|pos| (*pos, merkle_tree.get_leaf(*pos - first_leaf_pos)))
            .collect();

        assert!(batch_proof.verify::<TestBackend>(merkle_tree.root, proven_leaves_values_hashed));
    }

    #[test]
    fn batch_proof_big_tree_one_leaf() {
        let leaves_values: Vec<Ecgfp5FE> = (1..u64::pow(2, 16) + 1).map(Ecgfp5FE::new).collect();
        let merkle_tree: MerkleTree<TestBackend> = TestMerkleTreeEcgfp::build(&leaves_values);

        let proven_leaves_indices = [76]; // Only prove one of the leaves
        let first_leaf_pos = merkle_tree.nodes_len() / 2;
        let proven_leaves_positions: Vec<usize> = proven_leaves_indices
            .iter()
            .map(|leaf_index| leaf_index + first_leaf_pos)
            .collect();

        let batch_proof = merkle_tree
            .get_batch_proof(&proven_leaves_positions)
            .unwrap();

        let proven_leaves_values_hashed: HashMap<usize, Ecgfp5FE> = proven_leaves_positions
            .iter()
            .map(|pos| (*pos, merkle_tree.get_leaf(*pos - first_leaf_pos)))
            .collect();

        assert!(batch_proof.verify::<TestBackend>(merkle_tree.root, proven_leaves_values_hashed));
    }

    #[test]
    fn batch_proof_big_tree_many_leaves() {
        // Just add -18 to make the test case a little more complex
        let all_leaves_values: Vec<Ecgfp5FE> =
            (1..u64::pow(2, 16) - 18).map(Ecgfp5FE::new).collect();
        let merkle_tree: MerkleTree<TestBackend> = TestMerkleTreeEcgfp::build(&all_leaves_values);

        let proven_leaves_indices = usize::pow(2, 4) + 5..(usize::pow(2, 13) + 7);
        let first_leaf_pos = merkle_tree.nodes_len() / 2;
        let proven_leaves_positions: Vec<usize> = proven_leaves_indices
            .map(|leaf_index| leaf_index + first_leaf_pos)
            .collect();

        let batch_proof = merkle_tree
            .get_batch_proof(&proven_leaves_positions)
            .unwrap();

        let proven_leaves_values_hashed: HashMap<usize, Ecgfp5FE> = proven_leaves_positions
            .iter()
            .map(|pos| (*pos, merkle_tree.get_leaf(*pos - first_leaf_pos)))
            .collect();

        assert!(batch_proof.verify::<TestBackend>(merkle_tree.root, proven_leaves_values_hashed));
    }

    #[test]
    fn create_a_merkle_tree_with_10000_elements_and_verify_that_a_series_of_elements_belong_to_it()
    {
        let all_leaves_values: Vec<Ecgfp5FE> = (1..10000).map(Ecgfp5FE::new).collect();
        let merkle_tree = TestMerkleTreeEcgfp::build(&all_leaves_values);

        let proven_leaves_indices = [0].iter();
        let first_leaf_pos = merkle_tree.nodes_len() / 2;
        let proven_leaves_positions: Vec<usize> = proven_leaves_indices
            .clone()
            .map(|leaf_index| leaf_index + first_leaf_pos)
            .collect();

        let batch_proof = merkle_tree
            .get_batch_proof(&proven_leaves_positions)
            .unwrap();

        let proven_leaves_values_hashed: HashMap<usize, Ecgfp5FE> = proven_leaves_positions
            .iter()
            .map(|pos| (*pos, merkle_tree.get_leaf(*pos - first_leaf_pos)))
            .collect();

        assert!(batch_proof.verify::<TestBackend>(merkle_tree.root, proven_leaves_values_hashed));
    }

    #[cfg(feature = "serde")]
    #[test]
    fn proof_is_same_after_serde() {
        let all_leaves_values: Vec<Ecgfp5FE> = (1..10000).map(Ecgfp5FE::new).collect();
        let merkle_tree = TestMerkleTreeEcgfp::build(&all_leaves_values);

        let proven_leaves_indices = [0].iter();
        let first_leaf_pos = merkle_tree.nodes_len() / 2;
        let proven_leaves_positions: Vec<usize> = proven_leaves_indices
            .clone()
            .map(|leaf_index| leaf_index + first_leaf_pos)
            .collect();

        let batch_proof = merkle_tree
            .get_batch_proof(&proven_leaves_positions)
            .unwrap();

        let serialized = batch_proof.serialize();
        let batch_proof_2 = super::BatchProof::<Ecgfp5FE>::deserialize(&serialized);

        assert_eq!(batch_proof, batch_proof_2);
    }
}
