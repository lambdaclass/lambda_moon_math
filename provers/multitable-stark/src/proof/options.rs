use super::errors::InsecureOptionError;
use lambdaworks_math::field::traits::IsPrimeField;

#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::wasm_bindgen;

pub enum SecurityLevel {
    Conjecturable80Bits,
    Conjecturable100Bits,
    Conjecturable128Bits,
    Provable80Bits,
    Provable100Bits,
    Provable128Bits,
}

/// The options for the proof
///
/// - `blowup_factor`: the blowup factor for the trace
/// - `fri_number_of_queries`: the number of queries for the FRI layer
/// - `coset_offset`: the offset for the coset
/// - `grinding_factor`: the number of leading zeros that we want for the Hash(hash || nonce)
#[cfg_attr(feature = "wasm", wasm_bindgen)]
#[derive(Clone, Debug)]
pub struct ProofOptions {
    pub blowup_factor: u8,
    pub fri_number_of_queries: usize,
    pub coset_offset: u64,
    pub grinding_factor: u8,
}

impl ProofOptions {
    // TODO: Make it work for extended fields
    const EXTENSION_DEGREE: usize = 1;
    // Estimated maximum domain size. 2^40 = 1 TB
    const NUM_BITS_MAX_DOMAIN_SIZE: usize = 40;

    /// See section 5.10.1 of https://eprint.iacr.org/2021/582.pdf
    pub fn new_secure(security_level: SecurityLevel, coset_offset: u64) -> Self {
        match security_level {
            SecurityLevel::Conjecturable80Bits => ProofOptions {
                blowup_factor: 4,
                fri_number_of_queries: 31,
                coset_offset,
                grinding_factor: 20,
            },
            SecurityLevel::Conjecturable100Bits => ProofOptions {
                blowup_factor: 4,
                fri_number_of_queries: 41,
                coset_offset,
                grinding_factor: 20,
            },
            SecurityLevel::Conjecturable128Bits => ProofOptions {
                blowup_factor: 4,
                fri_number_of_queries: 55,
                coset_offset,
                grinding_factor: 20,
            },
            SecurityLevel::Provable80Bits => ProofOptions {
                blowup_factor: 4,
                fri_number_of_queries: 80,
                coset_offset,
                grinding_factor: 20,
            },
            SecurityLevel::Provable100Bits => ProofOptions {
                blowup_factor: 4,
                fri_number_of_queries: 104,
                coset_offset,
                grinding_factor: 20,
            },
            SecurityLevel::Provable128Bits => ProofOptions {
                blowup_factor: 4,
                fri_number_of_queries: 140,
                coset_offset,
                grinding_factor: 20,
            },
        }
    }

    /// Checks security of proof options given 128 bits of security
    pub fn new_with_checked_security<F: IsPrimeField>(
        blowup_factor: u8,
        fri_number_of_queries: usize,
        coset_offset: u64,
        grinding_factor: u8,
        security_target: u8,
    ) -> Result<Self, InsecureOptionError> {
        Self::check_field_security::<F>(security_target)?;

        let num_bits_blowup_factor = blowup_factor.trailing_zeros() as usize;

        if security_target as usize
            >= grinding_factor as usize + num_bits_blowup_factor * fri_number_of_queries - 1
        {
            return Err(InsecureOptionError::LowSecurityBits);
        }

        Ok(ProofOptions {
            blowup_factor,
            fri_number_of_queries,
            coset_offset,
            grinding_factor,
        })
    }

    /// Checks provable security of proof options given 128 bits of security
    /// This is an approximation. It's stricter than the formula in the paper.
    /// See https://eprint.iacr.org/2021/582.pdf
    pub fn new_with_checked_provable_security<F: IsPrimeField>(
        blowup_factor: u8,
        fri_number_of_queries: usize,
        coset_offset: u64,
        grinding_factor: u8,
        security_target: u8,
    ) -> Result<Self, InsecureOptionError> {
        Self::check_field_security::<F>(security_target)?;

        let num_bits_blowup_factor = blowup_factor.leading_zeros() as usize;

        if (security_target as usize)
            < grinding_factor as usize + num_bits_blowup_factor * fri_number_of_queries / 2
        {
            return Err(InsecureOptionError::LowSecurityBits);
        }

        Ok(ProofOptions {
            blowup_factor,
            fri_number_of_queries,
            coset_offset,
            grinding_factor,
        })
    }

    fn check_field_security<F: IsPrimeField>(
        security_target: u8,
    ) -> Result<(), InsecureOptionError> {
        if F::field_bit_size() * Self::EXTENSION_DEGREE
            <= security_target as usize + Self::NUM_BITS_MAX_DOMAIN_SIZE
        {
            return Err(InsecureOptionError::FieldSize);
        }

        Ok(())
    }

    /// Default proof options used for testing purposes.
    /// These options should never be used in production.
    pub fn default_test_options() -> Self {
        Self {
            blowup_factor: 4,
            fri_number_of_queries: 3,
            coset_offset: 3,
            grinding_factor: 1,
        }
    }
}

#[cfg(test)]
mod tests {
    use lambdaworks_math::field::fields::{
        fft_friendly::stark_252_prime_field::Stark252PrimeField, u64_prime_field::F17,
    };

    use crate::proof::{errors::InsecureOptionError, options::SecurityLevel};

    use super::ProofOptions;

    #[test]
    fn u64_prime_field_is_not_large_enough_to_be_secure() {
        let ProofOptions {
            blowup_factor,
            fri_number_of_queries,
            coset_offset,
            grinding_factor,
        } = ProofOptions::new_secure(SecurityLevel::Conjecturable128Bits, 1);

        let u64_options = ProofOptions::new_with_checked_security::<F17>(
            blowup_factor,
            fri_number_of_queries,
            coset_offset,
            grinding_factor,
            128,
        );

        assert!(matches!(u64_options, Err(InsecureOptionError::FieldSize)));
    }

    #[test]
    fn generated_stark_proof_options_for_128_bits_are_secure() {
        let ProofOptions {
            blowup_factor,
            fri_number_of_queries,
            coset_offset,
            grinding_factor,
        } = ProofOptions::new_secure(SecurityLevel::Conjecturable128Bits, 1);

        let secure_options = ProofOptions::new_with_checked_security::<Stark252PrimeField>(
            blowup_factor,
            fri_number_of_queries,
            coset_offset,
            grinding_factor,
            128,
        );

        assert!(secure_options.is_ok());
    }

    #[test]
    fn generated_proof_options_for_128_bits_with_one_fri_query_less_are_insecure() {
        let ProofOptions {
            blowup_factor,
            fri_number_of_queries,
            coset_offset,
            grinding_factor,
        } = ProofOptions::new_secure(SecurityLevel::Conjecturable128Bits, 1);

        let insecure_options = ProofOptions::new_with_checked_security::<Stark252PrimeField>(
            blowup_factor,
            fri_number_of_queries - 1,
            coset_offset,
            grinding_factor,
            128,
        );

        assert!(matches!(
            insecure_options,
            Err(InsecureOptionError::LowSecurityBits)
        ));
    }

    #[test]
    fn generated_stark_proof_options_for_100_bits_are_secure_for_100_target_bits() {
        let ProofOptions {
            blowup_factor,
            fri_number_of_queries,
            coset_offset,
            grinding_factor,
        } = ProofOptions::new_secure(SecurityLevel::Conjecturable100Bits, 1);

        let secure_options = ProofOptions::new_with_checked_security::<Stark252PrimeField>(
            blowup_factor,
            fri_number_of_queries,
            coset_offset,
            grinding_factor,
            100,
        );

        assert!(secure_options.is_ok());
    }

    #[test]
    fn generated_stark_proof_options_for_80_bits_are_secure_for_80_target_bits() {
        let ProofOptions {
            blowup_factor,
            fri_number_of_queries,
            coset_offset,
            grinding_factor,
        } = ProofOptions::new_secure(SecurityLevel::Conjecturable80Bits, 1);

        let secure_options = ProofOptions::new_with_checked_security::<Stark252PrimeField>(
            blowup_factor,
            fri_number_of_queries,
            coset_offset,
            grinding_factor,
            80,
        );

        assert!(secure_options.is_ok());
    }
}
