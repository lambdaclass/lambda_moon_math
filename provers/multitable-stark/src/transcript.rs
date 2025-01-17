use lambdaworks_crypto::fiat_shamir::is_transcript::IsTranscript;
use lambdaworks_math::{
    field::{
        element::FieldElement, fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
        traits::IsFFTField,
    },
    traits::{AsBytes, ByteConversion},
    unsigned_integer::element::U256,
};
use sha3::{Digest, Keccak256};

/// A transcript implementing `IsTranscript` and compatible with Stone (https://github.com/starkware-libs/stone-prover).
pub struct StoneProverTranscript {
    state: [u8; 32],
    seed_increment: U256,
    counter: u32,
    spare_bytes: Vec<u8>,
}

impl StoneProverTranscript {
    /// The maximum multiple of the modulus of `p` that fits in 256 bits, where
    /// `p = 0x800000000000011000000000000000000000000000000000000000000000001`
    const MODULUS_MAX_MULTIPLE: U256 = U256::from_hex_unchecked(
        "f80000000000020f00000000000000000000000000000000000000000000001f",
    );

    /// The value of `2^{-256} mod p`, where
    /// `p = 0x800000000000011000000000000000000000000000000000000000000000001`
    const R_INV: U256 = U256::from_hex_unchecked(
        "0x40000000000001100000000000012100000000000000000000000000000000",
    );
    pub fn new(public_input_data: &[u8]) -> Self {
        StoneProverTranscript {
            state: Self::keccak_hash(public_input_data),
            seed_increment: U256::from_hex_unchecked("1"),
            counter: 0,
            spare_bytes: vec![],
        }
    }

    pub fn sample_block(&mut self, used_bytes: usize) -> Vec<u8> {
        let mut first_part: Vec<u8> = self.state.to_vec();
        let mut counter_bytes: Vec<u8> = vec![0; 28]
            .into_iter()
            .chain(self.counter.to_be_bytes().to_vec())
            .collect();
        self.counter += 1;
        first_part.append(&mut counter_bytes);
        let block = Self::keccak_hash(&first_part);
        self.spare_bytes.extend(&block[used_bytes..]);
        block[..used_bytes].to_vec()
    }

    pub fn sample(&mut self, num_bytes: usize) -> Vec<u8> {
        let num_blocks = num_bytes / 32;
        let mut result: Vec<u8> = Vec::new();

        for _ in 0..num_blocks {
            let mut block = self.sample_block(32);
            result.append(&mut block);
        }

        let rest = num_bytes % 32;
        if rest <= self.spare_bytes.len() {
            result.append(&mut self.spare_bytes[..rest].to_vec());
            self.spare_bytes.drain(..rest);
        } else {
            let mut block = self.sample_block(rest);
            result.append(&mut block);
        }
        result
    }

    pub fn sample_big_int(&mut self) -> U256 {
        U256::from_bytes_be(&self.sample(32)).unwrap()
    }

    fn keccak_hash(data: &[u8]) -> [u8; 32] {
        let mut hasher = Keccak256::new();
        hasher.update(data);
        let mut result_hash = [0_u8; 32];
        result_hash.copy_from_slice(&hasher.finalize_reset());
        result_hash
    }
}

impl IsTranscript<Stark252PrimeField> for StoneProverTranscript {
    fn append_field_element(&mut self, element: &FieldElement<Stark252PrimeField>) {
        let limbs = element.value().limbs;
        let mut bytes: [u8; 32] = [0; 32];

        for i in (0..4).rev() {
            let limb_bytes = limbs[i].to_be_bytes();
            for j in 0..8 {
                bytes[i * 8 + j] = limb_bytes[j]
            }
        }
        self.append_bytes(&bytes);
    }

    fn append_bytes(&mut self, new_bytes: &[u8]) {
        let mut result_hash = [0_u8; 32];
        result_hash.copy_from_slice(&self.state);
        result_hash.reverse();

        let digest = U256::from_bytes_be(&self.state).unwrap();
        let new_seed = (digest + self.seed_increment).to_bytes_be();
        self.state = Self::keccak_hash(&[&new_seed, new_bytes].concat());
        self.counter = 0;
        self.spare_bytes.clear();
    }

    fn state(&self) -> [u8; 32] {
        self.state
    }

    fn sample_field_element(&mut self) -> FieldElement<Stark252PrimeField> {
        let mut result = self.sample_big_int();
        while result >= Self::MODULUS_MAX_MULTIPLE {
            result = self.sample_big_int();
        }
        FieldElement::<Stark252PrimeField>::new(result) * FieldElement::new(Self::R_INV)
    }

    fn sample_u64(&mut self, upper_bound: u64) -> u64 {
        // assert!(upper_bound < (1 << 12));
        let mut bytes = [0u8; 8];
        bytes.copy_from_slice(&self.sample(8));
        let u64_val: u64 = u64::from_be_bytes(bytes);
        u64_val % upper_bound
    }
}

/// Returns a batch of size `size` of field elements sampled from the transcript `transcript`.
pub fn batch_sample_challenges<F: IsFFTField>(
    size: usize,
    transcript: &mut impl IsTranscript<F>,
) -> Vec<FieldElement<F>>
where
    FieldElement<F>: AsBytes,
{
    (0..size)
        .map(|_| transcript.sample_field_element())
        .collect()
}

#[cfg(test)]
mod tests {
    use lambdaworks_math::field::{
        element::FieldElement, fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
    };

    use crate::transcript::{IsTranscript, StoneProverTranscript};

    use std::num::ParseIntError;

    type FE = FieldElement<Stark252PrimeField>;

    pub fn decode_hex(s: &str) -> Result<Vec<u8>, ParseIntError> {
        (0..s.len())
            .step_by(2)
            .map(|i| u8::from_str_radix(&s[i..i + 2], 16))
            .collect()
    }

    #[test]
    fn sample_bytes_from_stone_prover_channel() {
        let mut transcript = StoneProverTranscript::new(&[0x01, 0x02, 0x03]);
        transcript.append_bytes(&[0x04, 0x05, 0x06]);
        assert_eq!(
            transcript.sample(32),
            vec![
                0x8a, 0x3a, 0x67, 0xd1, 0x25, 0xa5, 0xa5, 0xea, 0x57, 0xc3, 0xfb, 0xe2, 0xc2, 0x55,
                0xb6, 0x0d, 0x0c, 0x89, 0x13, 0xa6, 0x27, 0x13, 0xe0, 0x99, 0xb3, 0x77, 0xc6, 0xc2,
                0x9a, 0x21, 0x85, 0x97
            ]
        );
        assert_eq!(
            transcript.sample(64),
            vec![
                0x56, 0xde, 0x56, 0x2a, 0xfd, 0x98, 0x19, 0xb9, 0xaa, 0xa0, 0x1b, 0x16, 0xf4, 0xeb,
                0x33, 0x71, 0xd5, 0xd8, 0x0f, 0x35, 0x29, 0xd8, 0xc1, 0x7a, 0x4b, 0xf4, 0x10, 0xe3,
                0x19, 0xb7, 0x64, 0x4a, 0xd2, 0x1c, 0xff, 0x14, 0x3d, 0xfd, 0xca, 0x32, 0x2c, 0x59,
                0xa3, 0x47, 0x5d, 0xd0, 0x34, 0xdf, 0x6d, 0xa7, 0x0c, 0xf5, 0xd2, 0x6a, 0xdd, 0x65,
                0xe0, 0x6d, 0x1e, 0x4f, 0xc7, 0x39, 0x52, 0x32
            ]
        );
        assert_eq!(
            transcript.sample(48),
            vec![
                0xe4, 0xb6, 0x3c, 0xfc, 0x03, 0xc9, 0x82, 0x8b, 0x63, 0x53, 0xb9, 0xad, 0x73, 0x6d,
                0x23, 0x88, 0x4c, 0x07, 0xb4, 0x9d, 0xf1, 0x1d, 0xef, 0xb9, 0x53, 0xfa, 0x02, 0xb5,
                0x3c, 0x43, 0xcf, 0xa3, 0x30, 0x5a, 0x02, 0x7e, 0xa6, 0x5e, 0x3c, 0x86, 0x3d, 0xdb,
                0x48, 0xea, 0x73, 0xbf, 0xdf, 0xab
            ]
        );
        assert_eq!(
            transcript.sample(32),
            vec![
                0x82, 0xe1, 0xd4, 0xf8, 0xf0, 0x61, 0xa4, 0x17, 0x4b, 0xed, 0x58, 0x4e, 0xb5, 0x73,
                0x26, 0xb7, 0x63, 0x10, 0x37, 0x97, 0xbe, 0x0b, 0x57, 0xaf, 0x74, 0xfe, 0x33, 0x19,
                0xbd, 0xe5, 0x53, 0x21,
            ]
        );
        assert_eq!(
            transcript.sample(16),
            vec![
                0xb0, 0xc6, 0x7a, 0x04, 0x19, 0x0a, 0x25, 0x72, 0xa8, 0x2e, 0xfa, 0x97, 0x92, 0x44,
                0x73, 0xe9
            ]
        );
        assert_eq!(
            transcript.sample(8),
            vec![0xbd, 0x41, 0x28, 0xdd, 0x3a, 0xbc, 0x66, 0x18]
        );
        assert_eq!(
            transcript.sample(32),
            vec![
                0xcb, 0x66, 0xc9, 0x72, 0x39, 0x85, 0xe8, 0x7c, 0x30, 0xe1, 0xc7, 0x1d, 0x2f, 0x83,
                0x4a, 0xcd, 0x33, 0x85, 0xfb, 0xd5, 0x40, 0x69, 0x22, 0x6e, 0xc0, 0xf1, 0x8c, 0x40,
                0x26, 0x2f, 0x5f, 0x7c,
            ]
        );
        transcript.append_bytes(&[0x03, 0x02]);
        assert_eq!(
            transcript.sample(32),
            vec![
                0x69, 0x63, 0x72, 0x01, 0x84, 0x8b, 0x22, 0x82, 0xa6, 0x14, 0x6d, 0x47, 0xbb, 0xa9,
                0xa3, 0xc8, 0xdc, 0x1b, 0x8e, 0x2e, 0x2e, 0x21, 0x87, 0x77, 0xac, 0xe0, 0x3e, 0xce,
                0x6e, 0xa7, 0x9e, 0xb0,
            ]
        );
    }

    #[test]
    fn test_sample_bytes() {
        let mut transcript = StoneProverTranscript::new(&[0x01, 0x02]);
        assert_eq!(
            transcript.sample(8),
            vec![89, 27, 84, 161, 127, 200, 195, 181]
        );
    }

    #[test]
    fn test_sample_field_element() {
        let mut transcript = StoneProverTranscript::new(&[0x01, 0x02]);
        assert_eq!(
            transcript.sample_field_element(),
            FE::from_hex_unchecked(
                "20b962ed1a29c942e11dc63c00b51de816bcd8bf9acd221f3fa55e5201d69be"
            )
        );
    }

    #[test]
    fn test_sample_u64_element() {
        let mut transcript = StoneProverTranscript::new(&[0x01, 0x02]);
        assert_eq!(transcript.sample_u64(1024), 949);
    }

    #[test]
    fn test_sample_u64_after_appending_and_sampling_bytes() {
        let mut transcript = StoneProverTranscript::new(&[0x01, 0x02]);
        transcript.append_bytes(&[0x01, 0x02]);
        assert_eq!(transcript.sample(4), vec![0x06, 0xe5, 0x36, 0xf5]);
        assert_eq!(transcript.sample_u64(16), 5);
    }

    #[test]
    fn test_transcript_compatibility_with_stone_prover_1() {
        // This corresponds to the following run.
        // Air: `Fibonacci2ColsShifted`
        // `trace_length`: 4
        // `blowup_factor`: 4
        // `fri_number_of_queries`: 1
        let mut transcript = StoneProverTranscript::new(&[0xca, 0xfe, 0xca, 0xfe]);
        // Send hash of trace commitment
        transcript.append_bytes(
            &decode_hex("0eb9dcc0fb1854572a01236753ce05139d392aa3aeafe72abff150fe21175594")
                .unwrap(),
        );
        // Sample challenge to collapse the constraints for the composition polynomial
        assert_eq!(
            transcript.sample_field_element(),
            FE::from_hex_unchecked(
                "86105fff7b04ed4068ecccb8dbf1ed223bd45cd26c3532d6c80a818dbd4fa7"
            )
        );
        // Send hash of composition poly commitment H
        transcript.append_bytes(
            &decode_hex("7cdd8d5fe3bd62254a417e2e260e0fed4fccdb6c9005e828446f645879394f38")
                .unwrap(),
        );
        // Sample challenge Z to compute t_j(z), H(z)
        assert_eq!(
            transcript.sample_field_element(),
            FE::from_hex_unchecked(
                "317629e783794b52cd27ac3a5e418c057fec9dd42f2b537cdb3f24c95b3e550"
            )
        );
        // Append t_j(z), H(z)
        transcript.append_field_element(&FE::from_hex_unchecked(
            "70d8181785336cc7e0a0a1078a79ee6541ca0803ed3ff716de5a13c41684037",
        ));
        transcript.append_field_element(&FE::from_hex_unchecked(
            "29808fc8b7480a69295e4b61600480ae574ca55f8d118100940501b789c1630",
        ));
        transcript.append_field_element(&FE::from_hex_unchecked(
            "7d8110f21d1543324cc5e472ab82037eaad785707f8cae3d64c5b9034f0abd2",
        ));
        transcript.append_field_element(&FE::from_hex_unchecked(
            "1b58470130218c122f71399bf1e04cf75a6e8556c4751629d5ce8c02cc4e62d",
        ));
        transcript.append_field_element(&FE::from_hex_unchecked(
            "1c0b7c2275e36d62dfb48c791be122169dcc00c616c63f8efb2c2a504687e85",
        ));
        // Sample challenge Gamma to collapse the terms of the deep composition polynomial (batch open).
        // Powers of this challenge are used if more than two terms.
        assert_eq!(
            transcript.sample_field_element(),
            FE::from_hex_unchecked(
                "a0c79c1c77ded19520873d9c2440451974d23302e451d13e8124cf82fc15dd"
            )
        );
        // FRI: Sample challenge Zeta to split the polynomial in half
        assert_eq!(
            transcript.sample_field_element(),
            FE::from_hex_unchecked(
                "5c6b5a66c9fda19f583f0b10edbaade98d0e458288e62c2fa40e3da2b293cef"
            )
        );
        // FRI: Send hash of commitment at Layer 1
        transcript.append_bytes(
            &decode_hex("49c5672520e20eccc72aa28d6fa0d7ef446f1ede38d7c64fbb95d0f34a281803")
                .unwrap(),
        );
        // FRI: Sample challenge to split the polynomial in half
        assert_eq!(
            transcript.sample_field_element(),
            FE::from_hex_unchecked(
                "4243ca9a618e2127590af8e1b38c63a156863fe95e4211cc1ade9b50667bbfa"
            )
        );
        // Send field element at final layer of FRI
        transcript.append_field_element(&FE::from_hex_unchecked(
            "702ddae5809ad82a82556eed2d202202d770962b7d4d82581e183df3efa2da6",
        ));
        // Send proof of work
        transcript.append_bytes(&[0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x30, 0x4d]); // Eight bytes
                                                                                    // Sample query indices
        assert_eq!(transcript.sample_u64(8), 0);
    }

    #[test]
    fn test_transcript_compatibility_with_stone_prover_2() {
        // This corresponds to the following run.
        // Air: `Fibonacci2ColsShifted`
        // `trace_length`: 4
        // `blowup_factor`: 64
        // `fri_number_of_queries`: 2
        let mut transcript = StoneProverTranscript::new(&[0xfa, 0xfa, 0xfa, 0xee]);
        // Send hash of trace commitment
        transcript.append_bytes(
            &decode_hex("99d8d4342895c4e35a084f8ea993036be06f51e7fa965734ed9c7d41104f0848")
                .unwrap(),
        );
        // Sample challenge to collapse the constraints for the composition polynomial
        assert_eq!(
            transcript.sample_field_element(),
            FE::from_hex_unchecked(
                "3fc675742e0692558bb95f36bd34bdfe050697ed0d849e5369808685e548441"
            )
        );
        // Send hash of composition poly commitment H(z)
        transcript.append_bytes(
            &decode_hex("2f4b599828a3f1ac458202ce06ec223bc9f4ad9ac758030109d40eebcf5776fd")
                .unwrap(),
        );
        // Sample challenge Z to compute t_j(z), H(z)
        assert_eq!(
            transcript.sample_field_element(),
            FE::from_hex_unchecked(
                "7298af9e2574933e62e51b107b8ef52f253d20644fc7250e9af118b02bc8a71"
            )
        );
        // Append t_j(z), H(z)
        transcript.append_field_element(&FE::from_hex_unchecked(
            "6791c8cdbd981f7db9786d702b21b87f4128a6941f35683d8b10faafcab83d5",
        ));
        transcript.append_field_element(&FE::from_hex_unchecked(
            "3cd6d8a23d01db66ea4911d6d7b09595b674f0507278fbf1f15cd85aa4ba72d",
        ));
        transcript.append_field_element(&FE::from_hex_unchecked(
            "3123deded538b40c1faa7988310f315860a43e320ae70f8f86eaeadf3828a10",
        ));
        transcript.append_field_element(&FE::from_hex_unchecked(
            "4d2edcc28870d79cbbb87181ffcb5942f7fa1c7b5f5bd5794c43452700e00d7",
        ));
        transcript.append_field_element(&FE::from_hex_unchecked(
            "5c244407085950973147074ee245bd1c7ed6d8a019df997aab1928a4a9a1e19",
        ));
        // Sample challenge Gamma to collapse the terms of the deep composition polynomial (batch open).
        // Powers of this challenge are used if more than two terms.
        assert_eq!(
            transcript.sample_field_element(),
            FE::from_hex_unchecked(
                "12f2b9edda6bb334bdf340d99eb0e6815e57aabffb48359117f71e7d0159d93"
            )
        );
        // FRI: Sample challenge Zeta to split the polynomial in half
        assert_eq!(
            transcript.sample_field_element(),
            FE::from_hex_unchecked(
                "7549307d78354156552667acf19a0ae978d4ec4954d210e23d9979672987dc"
            )
        );
        // FRI: Send hash of commitment at Layer 1
        transcript.append_bytes(
            &decode_hex("97decf0ad3cd590e7e5a4f85b3d4fa8c02c6d4b5343388c4536127dc8ef0fbf2")
                .unwrap(),
        );
        // FRI: Sample challenge to split the polynomial in half
        assert_eq!(
            transcript.sample_field_element(),
            FE::from_hex_unchecked(
                "4b79e806108567fd0f670ded2be5468009aaefeb993b346579c4f295fa3ddd0"
            )
        );
        // Send field element at final layer of FRI
        transcript.append_field_element(&FE::from_hex_unchecked(
            "7b8aa43aef4d3f2d476608251cffc9fa1c655bedecbcac49e4cafb012c7edf4",
        ));
        // Send proof of work
        transcript.append_bytes(&[0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x3b, 0xb8]);
        assert_eq!(transcript.sample_u64(128), 28);
        assert_eq!(transcript.sample_u64(128), 31);
    }
}
