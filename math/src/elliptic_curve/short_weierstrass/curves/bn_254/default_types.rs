use crate::{
    field::{
        element::FieldElement,
        fields::montgomery_backed_prime_fields::{IsModulus, MontgomeryBackendPrimeField},
        traits::IsFFTField,
    },
    unsigned_integer::element::{UnsignedInteger, U256},
};

#[derive(Clone, Debug)]
pub struct FrConfig;

/// We define Fr where r is the number of points that the eliptic curve BN254 has.
/// I.e. r is the order of the group BN254 Curve.
/// r = 21888242871839275222246405745257275088548364400416034343698204186575808495617.
impl IsModulus<U256> for FrConfig {
    const MODULUS: U256 = U256::from_hex_unchecked(
        "30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001",
    );
}

/// FrField using MontgomeryBackend for Bn254
pub type FrField = MontgomeryBackendPrimeField<FrConfig, 4>;
/// FrElement using MontgomeryBackend for Bn254
pub type FrElement = FieldElement<FrField>;

/// TWO_ADICITY is 28 because there is a subgroup of Fr of order 2^28.
/// Note that 2^28 divides r - 1.
/// (r - 1) / 2^28 mod r = 81540058820840996586704275553141814055101440848469862132140264610111.
/// We calculated the TWO_ADIC_PRIMITVE_ROOT_OF_UNITY in the following way:
/// g = 5 is a primitive root of order r - 1, that is, g^{r-1} = 1 and g^i != 1 for i < r-1.
/// Then g^{(r-1) / 2^28} is a primitive root of order 2^28.
impl IsFFTField for FrField {
    const TWO_ADICITY: u64 = 28;
    const TWO_ADIC_PRIMITVE_ROOT_OF_UNITY: Self::BaseType = UnsignedInteger::from_hex_unchecked(
        "2A3C09F0A58A7E8500E0A7EB8EF62ABC402D111E41112ED49BD61B6E725B19F0",
    );
}
