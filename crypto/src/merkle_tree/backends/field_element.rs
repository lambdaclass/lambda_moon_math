use crate::hash::poseidon::Poseidon;

use crate::merkle_tree::traits::IsMerkleTreeBackend;
use core::marker::PhantomData;
use lambdaworks_math::{
    field::{element::FieldElement, traits::IsField},
    traits::AsBytes,
};
use sha3::{
    digest::{generic_array::GenericArray, OutputSizeUser},
    Digest,
};

#[derive(Clone)]
pub struct FieldElementBackend<F, D: Digest, const NUM_BYTES: usize> {
    phantom1: PhantomData<F>,
    phantom2: PhantomData<D>,
}

impl<F, D: Digest, const NUM_BYTES: usize> Default for FieldElementBackend<F, D, NUM_BYTES> {
    fn default() -> Self {
        Self {
            phantom1: PhantomData,
            phantom2: PhantomData,
        }
    }
}

impl<F, D: Digest, const NUM_BYTES: usize> IsMerkleTreeBackend
    for FieldElementBackend<F, D, NUM_BYTES>
where
    F: IsField,
    FieldElement<F>: AsBytes + Sync + Send,
    [u8; NUM_BYTES]: From<GenericArray<u8, <D as OutputSizeUser>::OutputSize>>,
{
    type Node = [u8; NUM_BYTES];
    type Data = FieldElement<F>;

    fn hash_data(input: &FieldElement<F>) -> [u8; NUM_BYTES] {
        let mut hasher = D::new();
        hasher.update(input.as_bytes());
        hasher.finalize().into()
    }

    fn hash_new_parent(left: &[u8; NUM_BYTES], right: &[u8; NUM_BYTES]) -> [u8; NUM_BYTES] {
        let mut hasher = D::new();
        hasher.update(left);
        hasher.update(right);
        hasher.finalize().into()
    }
}

#[derive(Clone, Default)]
pub struct TreePoseidon<P: Poseidon + Default> {
    _poseidon: PhantomData<P>,
}

impl<P> IsMerkleTreeBackend for TreePoseidon<P>
where
    P: Poseidon + Default,
    FieldElement<P::F>: Sync + Send,
{
    type Node = FieldElement<P::F>;
    type Data = FieldElement<P::F>;

    fn hash_data(input: &FieldElement<P::F>) -> FieldElement<P::F> {
        P::hash_single(input)
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
    use alloc::vec::Vec;
    use lambdaworks_math::field::{
        element::FieldElement, fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
    };
    use sha3::{Keccak256, Keccak512, Sha3_256, Sha3_512};

    use crate::merkle_tree::{backends::field_element::FieldElementBackend, merkle::MerkleTree};

    type F = Stark252PrimeField;
    type FE = FieldElement<F>;

    #[test]
    fn hash_data_field_element_backend_works_with_keccak_256() {
        let values: Vec<FE> = (1..6).map(FE::from).collect();

        let merkle_tree = MerkleTree::<FieldElementBackend<F, Keccak256, 32>>::build(&values);
        let proof = merkle_tree.get_proof(0).unwrap();
        assert!(proof.verify::<FieldElementBackend<F, Keccak256, 32>>(
            &merkle_tree.root,
            0,
            &values[0]
        ));
    }

    #[test]
    fn hash_data_field_element_backend_works_with_sha3_256() {
        let values: Vec<FE> = (1..6).map(FE::from).collect();
        let merkle_tree = MerkleTree::<FieldElementBackend<F, Sha3_256, 32>>::build(&values);
        let proof = merkle_tree.get_proof(0).unwrap();
        assert!(proof.verify::<FieldElementBackend<F, Sha3_256, 32>>(
            &merkle_tree.root,
            0,
            &values[0]
        ));
    }

    #[test]
    fn hash_data_field_element_backend_works_with_keccak_512() {
        let values: Vec<FE> = (1..6).map(FE::from).collect();
        let merkle_tree = MerkleTree::<FieldElementBackend<F, Keccak512, 64>>::build(&values);
        let proof = merkle_tree.get_proof(0).unwrap();
        assert!(proof.verify::<FieldElementBackend<F, Keccak512, 64>>(
            &merkle_tree.root,
            0,
            &values[0]
        ));
    }

    #[test]
    fn hash_data_field_element_backend_works_with_sha3_512() {
        let values: Vec<FE> = (1..6).map(FE::from).collect();
        let merkle_tree = MerkleTree::<FieldElementBackend<F, Sha3_512, 64>>::build(&values);
        let proof = merkle_tree.get_proof(0).unwrap();
        assert!(proof.verify::<FieldElementBackend<F, Sha3_512, 64>>(
            &merkle_tree.root,
            0,
            &values[0]
        ));
    }
}
