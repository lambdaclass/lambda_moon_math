use lambdaworks_math::{
    field::{
        element::FieldElement,
        fields::montgomery_backed_prime_fields::{IsModulus, U64PrimeField},
    },
    unsigned_integer::element::U64,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MontgomeryConfigU64Stark101Field;
impl IsModulus<U64> for MontgomeryConfigU64Stark101Field {
    //Babybear Prime p = 2^64 - 2^32 + 1
    const MODULUS: U64 = U64::from_u64(3221225473);
}

pub type U64Stark101PrimeField = U64PrimeField<MontgomeryConfigU64Stark101Field>;
