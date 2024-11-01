use lambdaworks_crypto::merkle_tree::{
    backends::types::{BatchKeccak256Backend, Keccak256Backend},
    merkle::MerkleTree,
};

// Merkle Trees configuration

// Security of both hashes should match

pub type FriMerkleTreeBackend<F> = Keccak256Backend<F>;
pub type FriMerkleTree<F> = MerkleTree<FriMerkleTreeBackend<F>>;

pub const COMMITMENT_SIZE: usize = 32;
pub type Commitment = [u8; COMMITMENT_SIZE];

pub type BatchedMerkleTreeBackend<F> = BatchKeccak256Backend<F>;
pub type BatchedMerkleTree<F> = MerkleTree<BatchedMerkleTreeBackend<F>>;
