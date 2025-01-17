#[cfg(feature = "parallel")]
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use sha3::{Digest, Keccak256};

const PREFIX: [u8; 8] = [0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xed];

/// Checks if the bit-string `Hash(Hash(prefix || seed || grinding_factor) || nonce)`
/// has at least `grinding_factor` zeros to the left.
/// `prefix` is the bit-string `0x123456789abcded`
///
/// # Parameters
///
/// * `seed`: the input seed,
/// * `nonce`: the value to be tested,
/// * `grinding_factor`: the number of leading zeros needed.
///
/// # Returns
///
/// `true` if the number of leading zeros is at least `grinding_factor`, and `false` otherwise.
pub fn is_valid_nonce(seed: &[u8; 32], nonce: u64, grinding_factor: u8) -> bool {
    let inner_hash = get_inner_hash(seed, grinding_factor);
    let limit = 1 << (64 - grinding_factor);
    is_valid_nonce_for_inner_hash(&inner_hash, nonce, limit)
}

/// Performs grinding, returning a new nonce for the proof.
/// The nonce generated is such that:
/// Hash(Hash(prefix || seed || grinding_factor) || nonce) has at least `grinding_factor` zeros
/// to the left.
/// `prefix` is the bit-string `0x123456789abcded`
///
/// # Parameters
///
/// * `seed`: the input seed,
/// * `grinding_factor`: the number of leading zeros needed.
///
/// # Returns
///
/// A `nonce` satisfying the required condition.
pub fn generate_nonce(seed: &[u8; 32], grinding_factor: u8) -> Option<u64> {
    let inner_hash = get_inner_hash(seed, grinding_factor);
    let limit = 1 << (64 - grinding_factor);

    #[cfg(not(feature = "parallel"))]
    return (0..u64::MAX).find(|&candidate_nonce| {
        is_valid_nonce_for_inner_hash(&inner_hash, candidate_nonce, limit)
    });

    #[cfg(feature = "parallel")]
    return (0..u64::MAX).into_par_iter().find_any(|&candidate_nonce| {
        is_valid_nonce_for_inner_hash(&inner_hash, candidate_nonce, limit)
    });
}

/// Checks if the leftmost 8 bytes of `Hash(inner_hash || candidate_nonce)` are less than `limit`
/// when interpreted as `u64`.
#[inline(always)]
fn is_valid_nonce_for_inner_hash(inner_hash: &[u8; 32], candidate_nonce: u64, limit: u64) -> bool {
    let mut data = [0; 40];
    data[..32].copy_from_slice(inner_hash);
    data[32..].copy_from_slice(&candidate_nonce.to_be_bytes());

    let digest = Keccak256::digest(data);

    let seed_head = u64::from_be_bytes(digest[..8].try_into().unwrap());
    seed_head < limit
}

/// Returns the bit-string constructed as
/// Hash(prefix || seed || grinding_factor)
/// `prefix` is the bit-string `0x123456789abcded`
fn get_inner_hash(seed: &[u8; 32], grinding_factor: u8) -> [u8; 32] {
    let mut inner_data = [0u8; 41];
    inner_data[0..8].copy_from_slice(&PREFIX);
    inner_data[8..40].copy_from_slice(seed);
    inner_data[40] = grinding_factor;

    let digest = Keccak256::digest(inner_data);
    digest[..32].try_into().unwrap()
}

#[cfg(test)]
mod test {
    use crate::grinding::is_valid_nonce;

    #[test]
    fn test_invalid_nonce_grinding_factor_6() {
        // This setting produces a hash with 5 leading zeros, therefore not enough for grinding
        // factor 6.
        let seed = [
            174, 187, 26, 134, 6, 43, 222, 151, 140, 48, 52, 67, 69, 181, 177, 165, 111, 222, 148,
            92, 130, 241, 171, 2, 62, 34, 95, 159, 37, 116, 155, 217,
        ];
        let nonce = 4;
        let grinding_factor = 6;
        assert!(!is_valid_nonce(&seed, nonce, grinding_factor));
    }

    #[test]
    fn test_invalid_nonce_grinding_factor_9() {
        // This setting produces a hash with 8 leading zeros, therefore not enough for grinding
        // factor 9.
        let seed = [
            174, 187, 26, 134, 6, 43, 222, 151, 140, 48, 52, 67, 69, 181, 177, 165, 111, 222, 148,
            92, 130, 241, 171, 2, 62, 34, 95, 159, 37, 116, 155, 217,
        ];
        let nonce = 287;
        let grinding_factor = 9;
        assert!(!is_valid_nonce(&seed, nonce, grinding_factor));
    }

    #[test]
    fn test_is_valid_nonce_grinding_factor_10() {
        let seed = [
            37, 68, 26, 150, 139, 142, 66, 175, 33, 47, 199, 160, 9, 109, 79, 234, 135, 254, 39,
            11, 225, 219, 206, 108, 224, 165, 25, 72, 189, 96, 218, 95,
        ];
        let nonce = 0x5ba;
        let grinding_factor = 10;
        assert!(is_valid_nonce(&seed, nonce, grinding_factor));
    }

    #[test]
    fn test_is_valid_nonce_grinding_factor_20() {
        let seed = [
            37, 68, 26, 150, 139, 142, 66, 175, 33, 47, 199, 160, 9, 109, 79, 234, 135, 254, 39,
            11, 225, 219, 206, 108, 224, 165, 25, 72, 189, 96, 218, 95,
        ];
        let nonce = 0x2c5db8;
        let grinding_factor = 20;
        assert!(is_valid_nonce(&seed, nonce, grinding_factor));
    }

    #[test]
    fn test_invalid_nonce_grinding_factor_19() {
        // This setting would pass for grinding factor 20 instead of 19. The nonce is invalid
        // here because the grinding factor is part of the inner hash, changing the outer hash
        // and the resulting number of leading zeros.
        let seed = [
            37, 68, 26, 150, 139, 142, 66, 175, 33, 47, 199, 160, 9, 109, 79, 234, 135, 254, 39,
            11, 225, 219, 206, 108, 224, 165, 25, 72, 189, 96, 218, 95,
        ];
        let nonce = 0x2c5db8;
        let grinding_factor = 19;
        assert!(!is_valid_nonce(&seed, nonce, grinding_factor));
    }

    #[test]
    fn test_is_valid_nonce_grinding_factor_30() {
        let seed = [
            37, 68, 26, 150, 139, 142, 66, 175, 33, 47, 199, 160, 9, 109, 79, 234, 135, 254, 39,
            11, 225, 219, 206, 108, 224, 165, 25, 72, 189, 96, 218, 95,
        ];
        let nonce = 0x1ae839e1;
        let grinding_factor = 30;
        assert!(is_valid_nonce(&seed, nonce, grinding_factor));
    }

    #[test]
    fn test_is_valid_nonce_grinding_factor_33() {
        let seed = [
            37, 68, 26, 150, 139, 142, 66, 175, 33, 47, 199, 160, 9, 109, 79, 234, 135, 254, 39,
            11, 225, 219, 206, 108, 224, 165, 25, 72, 189, 96, 218, 95,
        ];
        let nonce = 0x4cc3123f;
        let grinding_factor = 33;
        assert!(is_valid_nonce(&seed, nonce, grinding_factor));
    }
}
