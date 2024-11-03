use core::marker::PhantomData;

use crate::hash::poseidon::Poseidon;
use crate::merkle_tree::traits::IsMerkleTreeBackend;
use alloc::vec::Vec;
use lambdaworks_math::{
    field::{element::FieldElement, traits::IsField},
    traits::AsBytes,
};
use sha3::{
    digest::{generic_array::GenericArray, OutputSizeUser},
    Digest,
};

#[derive(Clone)]
pub struct FieldElementVectorBackend<F, D: Digest, const NUM_BYTES: usize> {
    phantom1: PhantomData<F>,
    phantom2: PhantomData<D>,
}

impl<F, D: Digest, const NUM_BYTES: usize> Default for FieldElementVectorBackend<F, D, NUM_BYTES> {
    fn default() -> Self {
        Self {
            phantom1: PhantomData,
            phantom2: PhantomData,
        }
    }
}

impl<F, D: Digest, const NUM_BYTES: usize> IsMerkleTreeBackend
    for FieldElementVectorBackend<F, D, NUM_BYTES>
where
    F: IsField,
    FieldElement<F>: AsBytes,
    [u8; NUM_BYTES]: From<GenericArray<u8, <D as OutputSizeUser>::OutputSize>>,
    Vec<FieldElement<F>>: Sync + Send,
{
    type Node = [u8; NUM_BYTES];
    type Data = Vec<FieldElement<F>>;

    fn hash_data(input: &Vec<FieldElement<F>>) -> [u8; NUM_BYTES] {
        let mut hasher = D::new();
        for element in input.iter() {
            hasher.update(element.as_bytes());
        }
        let mut result_hash = [0_u8; NUM_BYTES];
        result_hash.copy_from_slice(&hasher.finalize());
        result_hash
    }

    fn hash_new_parent(left: &[u8; NUM_BYTES], right: &[u8; NUM_BYTES]) -> [u8; NUM_BYTES] {
        let mut hasher = D::new();
        hasher.update(left);
        hasher.update(right);
        let mut result_hash = [0_u8; NUM_BYTES];
        result_hash.copy_from_slice(&hasher.finalize());
        result_hash
    }
}

#[derive(Clone, Default)]
pub struct BatchPoseidonTree<P: Poseidon + Default> {
    _poseidon: PhantomData<P>,
}

impl<P> IsMerkleTreeBackend for BatchPoseidonTree<P>
where
    P: Poseidon + Default,
    Vec<FieldElement<P::F>>: Sync + Send,
    FieldElement<P::F>: Sync + Send,
{
    type Node = FieldElement<P::F>;
    type Data = Vec<FieldElement<P::F>>;

    fn hash_data(input: &Vec<FieldElement<P::F>>) -> FieldElement<P::F> {
        P::hash_many(input)
    }

    fn hash_new_parent(
        left: &FieldElement<P::F>,
        right: &FieldElement<P::F>,
    ) -> FieldElement<P::F> {
        P::hash(left, right)
    }
}

#[cfg(test)]
mod tests {
    use lambdaworks_math::field::{
        element::FieldElement, fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
    };
    use sha2::Sha512;
    use sha3::{Keccak256, Keccak512, Sha3_256, Sha3_512};

    use crate::merkle_tree::{
        backends::field_element_vector::FieldElementVectorBackend, merkle::MerkleTree,
    };

    type F = Stark252PrimeField;
    type FE = FieldElement<F>;

    #[test]
    fn hash_data_field_element_backend_works_with_sha3_256() {
        let values = [
            vec![FE::from(2u64), FE::from(11u64)],
            vec![FE::from(3u64), FE::from(14u64)],
            vec![FE::from(4u64), FE::from(7u64)],
            vec![FE::from(5u64), FE::from(3u64)],
            vec![FE::from(6u64), FE::from(5u64)],
            vec![FE::from(7u64), FE::from(16u64)],
            vec![FE::from(8u64), FE::from(19u64)],
            vec![FE::from(9u64), FE::from(21u64)],
        ];
        let merkle_tree = MerkleTree::<FieldElementVectorBackend<F, Sha3_256, 32>>::build(&values);
        let proof = merkle_tree.get_proof(0).unwrap();
        assert!(proof.verify::<FieldElementVectorBackend<F, Sha3_256, 32>>(
            &merkle_tree.root,
            0,
            &values[0]
        ));
    }

    #[test]
    fn hash_data_field_element_backend_works_with_keccak256() {
        let values = [
            vec![FE::from(2u64), FE::from(11u64)],
            vec![FE::from(3u64), FE::from(14u64)],
            vec![FE::from(4u64), FE::from(7u64)],
            vec![FE::from(5u64), FE::from(3u64)],
            vec![FE::from(6u64), FE::from(5u64)],
            vec![FE::from(7u64), FE::from(16u64)],
            vec![FE::from(8u64), FE::from(19u64)],
            vec![FE::from(9u64), FE::from(21u64)],
        ];
        let merkle_tree = MerkleTree::<FieldElementVectorBackend<F, Keccak256, 32>>::build(&values);
        let proof = merkle_tree.get_proof(0).unwrap();
        assert!(proof.verify::<FieldElementVectorBackend<F, Keccak256, 32>>(
            &merkle_tree.root,
            0,
            &values[0]
        ));
    }

    #[test]
    fn hash_data_field_element_backend_works_with_sha3_512() {
        let values = [
            vec![FE::from(2u64), FE::from(11u64)],
            vec![FE::from(3u64), FE::from(14u64)],
            vec![FE::from(4u64), FE::from(7u64)],
            vec![FE::from(5u64), FE::from(3u64)],
            vec![FE::from(6u64), FE::from(5u64)],
            vec![FE::from(7u64), FE::from(16u64)],
            vec![FE::from(8u64), FE::from(19u64)],
            vec![FE::from(9u64), FE::from(21u64)],
        ];
        let merkle_tree = MerkleTree::<FieldElementVectorBackend<F, Sha3_512, 64>>::build(&values);
        let proof = merkle_tree.get_proof(0).unwrap();
        assert!(proof.verify::<FieldElementVectorBackend<F, Sha3_512, 64>>(
            &merkle_tree.root,
            0,
            &values[0]
        ));
    }

    #[test]
    fn hash_data_field_element_backend_works_with_keccak512() {
        let values = [
            vec![FE::from(2u64), FE::from(11u64)],
            vec![FE::from(3u64), FE::from(14u64)],
            vec![FE::from(4u64), FE::from(7u64)],
            vec![FE::from(5u64), FE::from(3u64)],
            vec![FE::from(6u64), FE::from(5u64)],
            vec![FE::from(7u64), FE::from(16u64)],
            vec![FE::from(8u64), FE::from(19u64)],
            vec![FE::from(9u64), FE::from(21u64)],
        ];
        let merkle_tree = MerkleTree::<FieldElementVectorBackend<F, Keccak512, 64>>::build(&values);
        let proof = merkle_tree.get_proof(0).unwrap();
        assert!(proof.verify::<FieldElementVectorBackend<F, Keccak512, 64>>(
            &merkle_tree.root,
            0,
            &values[0]
        ));
    }

    #[test]
    fn hash_data_field_element_backend_works_with_sha2_512() {
        let values = [
            vec![FE::from(2u64), FE::from(11u64)],
            vec![FE::from(3u64), FE::from(14u64)],
            vec![FE::from(4u64), FE::from(7u64)],
            vec![FE::from(5u64), FE::from(3u64)],
            vec![FE::from(6u64), FE::from(5u64)],
            vec![FE::from(7u64), FE::from(16u64)],
            vec![FE::from(8u64), FE::from(19u64)],
            vec![FE::from(9u64), FE::from(21u64)],
        ];
        let merkle_tree = MerkleTree::<FieldElementVectorBackend<F, Sha512, 64>>::build(&values);
        let proof = merkle_tree.get_proof(0).unwrap();
        assert!(proof.verify::<FieldElementVectorBackend<F, Sha512, 64>>(
            &merkle_tree.root,
            0,
            &values[0]
        ));
    }
}
