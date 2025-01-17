use lambdaworks_crypto::merkle_tree::proof::Proof;
use lambdaworks_math::errors::DeserializationError;

use super::config::Commitment;

pub fn serialize_proof(proof: &Proof<Commitment>) -> Vec<u8> {
    let mut bytes = vec![];
    bytes.extend(proof.merkle_path.len().to_be_bytes());
    for commitment in &proof.merkle_path {
        bytes.extend(commitment);
    }
    bytes
}

pub fn deserialize_proof(bytes: &[u8]) -> Result<(Proof<Commitment>, &[u8]), DeserializationError> {
    let mut bytes = bytes;
    let mut merkle_path = vec![];
    let merkle_path_len = usize::from_be_bytes(
        bytes
            .get(..8)
            .ok_or(DeserializationError::InvalidAmountOfBytes)?
            .try_into()
            .map_err(|_| DeserializationError::InvalidAmountOfBytes)?,
    );
    bytes = &bytes[8..];

    for _ in 0..merkle_path_len {
        let commitment = bytes
            .get(..32)
            .ok_or(DeserializationError::InvalidAmountOfBytes)?
            .try_into()
            .map_err(|_| DeserializationError::InvalidAmountOfBytes)?;
        merkle_path.push(commitment);
        bytes = &bytes[32..];
    }

    Ok((Proof { merkle_path }, bytes))
}
