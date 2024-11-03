use core::fmt::Display;

extern crate std;
use std::collections::HashSet;

use alloc::vec::Vec;

use super::{batch_proof::BatchProof, proof::Proof, traits::IsMerkleTreeBackend, utils::*};

const ROOT: usize = 0;

#[derive(Debug)]
pub enum Error {
    OutOfBounds,
}
impl Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Accessed node was out of bound")
    }
}

#[cfg(feature = "std")]
impl std::error::Error for Error {}

#[derive(Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MerkleTree<B: IsMerkleTreeBackend> {
    pub root: B::Node,
    nodes: Vec<B::Node>,
}

impl<B> MerkleTree<B>
where
    B: IsMerkleTreeBackend,
{
    pub fn build(unhashed_leaves: &[B::Data]) -> Option<Self> {
        if unhashed_leaves.is_empty() {
            return None;
        }

        let hashed_leaves: Vec<B::Node> = B::hash_leaves(unhashed_leaves);

        //The leaf must be a power of 2 set
        let hashed_leaves = complete_until_power_of_two(hashed_leaves);
        let leaves_len = hashed_leaves.len();

        //The length of leaves minus one inner node in the merkle tree
        //The first elements are overwritten by build function, it doesn't matter what it's there
        let mut nodes = vec![hashed_leaves[0].clone(); leaves_len - 1];
        nodes.extend(hashed_leaves);

        //Build the inner nodes of the tree
        build::<B>(&mut nodes, leaves_len);

        // #[cfg(test)]
        // print_positions(nodes.len(), HashSet::new());

        Some(MerkleTree {
            root: nodes[ROOT].clone(),
            nodes,
        })
    }

    /// Returns the leaf at the given index.
    /// leaf_index starts from 0 for the first leaf.
    pub fn get_leaf(&self, leaf_index: usize) -> B::Node {
        let first_leaf_pos = self.nodes_len() / 2;
        self.nodes[leaf_index + first_leaf_pos].clone()
    }

    pub fn get_proof(&self, pos: usize) -> Option<Proof<B::Node>> {
        let first_leaf_index = self.nodes_len() / 2;
        let pos = pos + first_leaf_index;
        let Ok(merkle_path) = self.build_merkle_path(pos) else {
            return None;
        };

        Some(Proof { merkle_path })
    }

    /// Builds and returns a batch proof for when a Merkle tree is used to prove inclusion of multiple leaves.
    ///
    /// pos_list is a list of leaf positions (within all tree) to create a batch inclusion proof for.
    /// pos_list need not be continuous, but the resulting proof becomes the smallest when so.
    pub fn get_batch_proof(&self, pos_list: &[usize]) -> Option<BatchProof<B::Node>> {
        let batch_auth_path_positions = self.get_batch_auth_path_positions(pos_list);

        let batch_auth_path_nodes_iter = batch_auth_path_positions
            .iter()
            .map(|pos| (*pos, self.nodes[*pos].clone()).clone());

        Some(BatchProof::<B::Node> {
            auth: batch_auth_path_nodes_iter.collect(),
        })
    }

    pub fn nodes_len(&self) -> usize {
        self.nodes.len()
    }

    /// Builds and returns a proof of inclusion for the leaf whose position is passed as an argument.
    ///
    /// pos parameter is the index in overall Merkle tree, including the inner nodes
    fn build_merkle_path(&self, pos: usize) -> Result<Vec<B::Node>, Error> {
        let mut merkle_path = Vec::new();
        let mut pos = pos;

        while pos != ROOT {
            let Some(node) = self.nodes.get(get_sibling_pos(pos)) else {
                // out of bounds, exit returning the current merkle_path
                return Err(Error::OutOfBounds);
            };
            merkle_path.push(node.clone());

            pos = get_parent_pos(pos);
        }

        Ok(merkle_path)
    }

    /// Batch Merkle proofs require multiple authentication paths to be computed, and some nodes in these paths
    /// can be obtained from the leaves that are subject to the batch Merkle proof.
    /// This function returns a set of node positions that are supposedly just enough to satisfy all authentication
    /// paths in the batch Merkle proof.
    ///
    /// See the following Merkle tree, where we build a batch authentication path for leaves [15,24] inclusively.
    /// We'd only need nodes (12) and (6), because all the other nodes that would be needed in the authentication
    /// path of any leaf can be obtained from just the leaves (which are public input).
    /// If we were to build a batch authentication path for leaves [15,26] inclusively, then we'd need (6) alone,
    /// because we could use nodes (25) and (26) to build (12), to be combined with (11) to obtain (5).
    ///
    /// 0
    /// 1                               2
    /// 3               4               5               6
    /// 7       8       9       10      11      12      13      14
    /// 15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30
    ///
    /// leaf_positions is a list of leaves to create a batch inclusion proof for.
    /// For the example above, it would be [15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
    /// leaf_positions need not be continuous, but the resulting auth_set becomes the smallest when so.
    fn get_batch_auth_path_positions(&self, leaf_positions: &[usize]) -> Vec<usize> {
        let mut auth_set = HashSet::<usize>::new();
        // Add all the leaves to the set of obtainable nodes, because we already have them.
        let mut obtainable_nodes_by_level: HashSet<usize> =
            leaf_positions.iter().cloned().collect();

        // Iterate levels starting from the leaves up to the root
        for _ in (1..self.levels()).rev() {
            let mut parent_level_obtainable_positions = HashSet::new();
            for pos in &obtainable_nodes_by_level {
                let sibling_pos = get_sibling_pos(*pos);

                let sibling_is_obtainable = obtainable_nodes_by_level.contains(&sibling_pos)
                    || auth_set.contains(&sibling_pos);
                if !sibling_is_obtainable {
                    auth_set.insert(sibling_pos);
                }
                parent_level_obtainable_positions.insert(get_parent_pos(*pos));
            }

            obtainable_nodes_by_level = parent_level_obtainable_positions;
        }

        auth_set.into_iter().collect()
    }

    fn levels(&self) -> usize {
        (self.nodes_len() as f32).log2().ceil() as usize
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use lambdaworks_math::field::{element::FieldElement, fields::u64_prime_field::U64PrimeField};

    use crate::merkle_tree::{merkle::MerkleTree, test_merkle::TestBackend};

    const MODULUS: u64 = 13;
    type U64PF = U64PrimeField<MODULUS>;
    type FE = FieldElement<U64PF>;
    type TestTree = MerkleTree<TestBackend<U64PF>>;

    #[test]
    fn build_merkle_tree_from_a_power_of_two_list_of_values() {
        let values: Vec<FE> = (1..5).map(FE::new).collect();
        let merkle_tree = TestTree::build(&values);
        assert_eq!(merkle_tree.root, FE::new(20));
    }

    #[test]
    // expected | 8 | 7 | 1 | 6 | 1 | 7 | 7 | 2 | 4 | 6 | 8 | 10 | 10 | 10 | 10 |
    fn build_merkle_tree_from_an_odd_set_of_leaves() {
        const MODULUS: u64 = 13;
        type U64PF = U64PrimeField<MODULUS>;
        type FE = FieldElement<U64PF>;

        let values: Vec<FE> = (1..6).map(FE::new).collect();

        let merkle_tree = TestTree::build(&values);
        assert_eq!(merkle_tree.root, FE::new(8));
    }

    #[test]
    fn build_merkle_tree_from_a_single_value() {
        const MODULUS: u64 = 13;
        type U64PF = U64PrimeField<MODULUS>;
        type FE = FieldElement<U64PF>;

        let values: Vec<FE> = vec![FE::new(1)]; // Single element
        let merkle_tree = MerkleTree::<TestBackend<U64PF>>::build(&values).unwrap();
        assert_eq!(merkle_tree.root, FE::new(2)); // Adjusted expected value
    }

    #[test]
    fn build_empty_tree_should_not_panic() {
        assert!(MerkleTree::<TestBackend<U64PF>>::build(&[]).is_none());
    }
}
