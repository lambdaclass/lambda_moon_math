use crate::traits::{AsBytes, ByteConversion};
use crate::{
    errors::CreationError,
    field::{
        element::FieldElement,
        errors::FieldError,
        traits::{IsField, IsPrimeField},
    },
};
use core::fmt::{self, Display};

/// Represents a 31 bit integer value
/// Invariants:
///      31st bit is clear
///      n < MODULUS
#[derive(Debug, Clone, Copy, Hash, PartialOrd, Ord, PartialEq, Eq)]
pub struct Mersenne31Field;

impl Mersenne31Field {
    fn weak_reduce(n: u32) -> u32 {
        // To reduce 'n' to 31 bits we clear its MSB, then add it back in its reduced form.
        let msb = n & (1 << 31);
        let msb_reduced = msb >> 31;
        let res = msb ^ n;

        // assert msb_reduced fits within 31 bits
        debug_assert!((res >> 31) == 0 && (msb_reduced >> 1) == 0);
        res + msb_reduced
    }

    fn as_representative(n: &u32) -> u32 {
        if *n == MERSENNE_31_PRIME_FIELD_ORDER {
            0
        } else {
            *n
        }
    }

    #[inline]
    pub fn sum<I: Iterator<Item = <Self as IsField>::BaseType>>(
        iter: I,
    ) -> <Self as IsField>::BaseType {
        // Delayed reduction
        Self::from_u64(iter.map(|x| (x as u64)).sum::<u64>())
    }

    /// Computes a * 2^k, with 0 < k < 31
    pub fn mul_power_two(a: u32, k: u32) -> u32 {
        let msb = (a & (u32::MAX << (31 - k))) >> (31 - k); // The k + 1 msf shifted right .
        let lsb = (a & (u32::MAX >> (k + 1))) << k; // The 31 - k lsb shifted left.
        Self::weak_reduce(msb + lsb)
    }

    pub fn pow_2(a: &u32, order: u32) -> u32 {
        let mut res = *a;
        (0..order).for_each(|_| res = Self::square(&res));
        res
    }

    /// TODO: See if we can optimize this function.
    /// Computes 2a^2 - 1
    pub fn two_square_minus_one(a: &u32) -> u32 {
        if *a == 0 {
            MERSENNE_31_PRIME_FIELD_ORDER - 1
        } else {
            Self::from_u64(((u64::from(*a) * u64::from(*a)) << 1) - 1)
        }
    }
}

pub const MERSENNE_31_PRIME_FIELD_ORDER: u32 = (1 << 31) - 1;

//NOTE: This implementation was inspired by and borrows from the work done by the Plonky3 team
// https://github.com/Plonky3/Plonky3/blob/main/mersenne-31/src/lib.rs
// Thank you for pushing this technology forward.
impl IsField for Mersenne31Field {
    type BaseType = u32;

    /// Returns the sum of `a` and `b`.
    fn add(a: &u32, b: &u32) -> u32 {
        // We are using that if a and b are field elements of Mersenne31, then
        // a + b has at most 32 bits, so we can use the weak_reduce function to take mudulus p.
        Self::weak_reduce(a + b)
    }

    /// Returns the multiplication of `a` and `b`.
    // Note: for powers of 2 we can perform bit shifting this would involve overriding the trait implementation
    fn mul(a: &u32, b: &u32) -> u32 {
        Self::from_u64(u64::from(*a) * u64::from(*b))
    }

    fn sub(a: &u32, b: &u32) -> u32 {
        Self::weak_reduce(a + MERSENNE_31_PRIME_FIELD_ORDER - b)
    }

    /// Returns the additive inverse of `a`.
    fn neg(a: &u32) -> u32 {
        // NOTE: MODULUS known to have 31 bit clear
        MERSENNE_31_PRIME_FIELD_ORDER - a
    }

    /// Returns the multiplicative inverse of `a`.
    fn inv(x: &u32) -> Result<u32, FieldError> {
        if *x == Self::zero() || *x == MERSENNE_31_PRIME_FIELD_ORDER {
            return Err(FieldError::InvZeroError);
        }
        let p101 = Self::mul(&Self::pow_2(x, 2), x);
        let p1111 = Self::mul(&Self::square(&p101), &p101);
        let p11111111 = Self::mul(&Self::pow_2(&p1111, 4u32), &p1111);
        let p111111110000 = Self::pow_2(&p11111111, 4u32);
        let p111111111111 = Self::mul(&p111111110000, &p1111);
        let p1111111111111111 = Self::mul(&Self::pow_2(&p111111110000, 4u32), &p11111111);
        let p1111111111111111111111111111 =
            Self::mul(&Self::pow_2(&p1111111111111111, 12u32), &p111111111111);
        let p1111111111111111111111111111101 =
            Self::mul(&Self::pow_2(&p1111111111111111111111111111, 3u32), &p101);
        Ok(p1111111111111111111111111111101)
    }

    /// Returns the division of `a` and `b`.
    fn div(a: &u32, b: &u32) -> u32 {
        let b_inv = Self::inv(b).expect("InvZeroError");
        Self::mul(a, &b_inv)
    }

    /// Returns a boolean indicating whether `a` and `b` are equal or not.
    fn eq(a: &u32, b: &u32) -> bool {
        Self::as_representative(a) == Self::representative(b)
    }

    /// Returns the additive neutral element.
    fn zero() -> u32 {
        0u32
    }

    /// Returns the multiplicative neutral element.
    fn one() -> u32 {
        1u32
    }

    /// Returns the element `x * 1` where 1 is the multiplicative neutral element.
    fn from_u64(x: u64) -> u32 {
        (((((x >> 31) + x + 1) >> 31) + x) & (MERSENNE_31_PRIME_FIELD_ORDER as u64)) as u32
    }

    /// Takes as input an element of BaseType and returns the internal representation
    /// of that element in the field.
    fn from_base_type(x: u32) -> u32 {
        Self::weak_reduce(x)
    }
    fn double(a: &u32) -> u32 {
        Self::weak_reduce(a << 1)
    }
}

impl IsPrimeField for Mersenne31Field {
    type RepresentativeType = u32;

    // Since our invariant guarantees that `value` fits in 31 bits, there is only one possible value
    // `value` that is not canonical, namely 2^31 - 1 = p = 0.
    fn representative(x: &u32) -> u32 {
        debug_assert!((x >> 31) == 0);
        Self::as_representative(x)
    }

    fn field_bit_size() -> usize {
        ((MERSENNE_31_PRIME_FIELD_ORDER - 1).ilog2() + 1) as usize
    }

    fn from_hex(hex_string: &str) -> Result<Self::BaseType, CreationError> {
        let mut hex_string = hex_string;
        // Remove 0x if it's on the string
        let mut char_iterator = hex_string.chars();
        if hex_string.len() > 2
            && char_iterator.next().unwrap() == '0'
            && char_iterator.next().unwrap() == 'x'
        {
            hex_string = &hex_string[2..];
        }
        u32::from_str_radix(hex_string, 16).map_err(|_| CreationError::InvalidHexString)
    }

    #[cfg(feature = "std")]
    fn to_hex(x: &u32) -> String {
        format!("{:X}", x)
    }
}

impl FieldElement<Mersenne31Field> {
    #[cfg(feature = "alloc")]
    pub fn to_bytes_le(&self) -> alloc::vec::Vec<u8> {
        self.representative().to_le_bytes().to_vec()
    }

    #[cfg(feature = "alloc")]
    pub fn to_bytes_be(&self) -> alloc::vec::Vec<u8> {
        self.representative().to_be_bytes().to_vec()
    }
}

impl Display for FieldElement<Mersenne31Field> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:x}", self.representative())
    }
}

impl AsBytes for FieldElement<Mersenne31Field> {
    fn as_bytes(&self) -> alloc::vec::Vec<u8> {
        self.to_bytes_be()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    type F = Mersenne31Field;
    type FE = FieldElement<F>;

    #[test]
    fn mul_power_two_is_correct() {
        let a = 3u32;
        let k = 2;
        let expected_result = FE::from(&a) * FE::from(2).pow(k as u16);
        let result = F::mul_power_two(a, k);
        assert_eq!(FE::from(&result), expected_result)
    }

    #[test]
    fn mul_power_two_is_correct_2() {
        let a = 229287u32;
        let k = 4;
        let expected_result = FE::from(&a) * FE::from(2).pow(k as u16);
        let result = F::mul_power_two(a, k);
        assert_eq!(FE::from(&result), expected_result)
    }

    #[test]
    fn pow_2_is_correct() {
        let a = 3u32;
        let order = 12;
        let result = F::pow_2(&a, order);
        let expected_result = FE::pow(&FE::from(&a), 4096u32);
        assert_eq!(FE::from(&result), expected_result)
    }

    #[test]
    fn from_hex_for_b_is_11() {
        assert_eq!(F::from_hex("B").unwrap(), 11);
    }

    #[test]
    fn from_hex_for_b_is_11_v2() {
        assert_eq!(FE::from_hex("B").unwrap(), FE::from(11));
    }

    #[test]
    fn sum_delayed_reduction() {
        let up_to = u32::pow(2, 16);
        let pow = u64::pow(2, 60);

        let iter = (0..up_to).map(F::weak_reduce).map(|e| F::pow(&e, pow));

        assert_eq!(F::from_u64(2142542785), F::sum(iter));
    }

    #[test]
    fn from_hex_for_0x1_a_is_26() {
        assert_eq!(F::from_hex("0x1a").unwrap(), 26);
    }

    #[test]
    fn bit_size_of_field_is_31() {
        assert_eq!(
            <F as crate::field::traits::IsPrimeField>::field_bit_size(),
            31
        );
    }

    #[test]
    fn one_plus_1_is_2() {
        assert_eq!(FE::one() + FE::one(), FE::from(&2u32));
    }

    #[test]
    fn neg_1_plus_1_is_0() {
        assert_eq!(-FE::one() + FE::one(), FE::zero());
    }

    #[test]
    fn neg_1_plus_2_is_1() {
        assert_eq!(-FE::one() + FE::from(&2u32), FE::one());
    }

    #[test]
    fn max_order_plus_1_is_0() {
        assert_eq!(
            FE::from(&(MERSENNE_31_PRIME_FIELD_ORDER - 1)) + FE::from(1),
            FE::from(0)
        );
    }

    #[test]
    fn comparing_13_and_13_are_equal() {
        assert_eq!(FE::from(&13u32), FE::from(13));
    }

    #[test]
    fn comparing_13_and_8_they_are_not_equal() {
        assert_ne!(FE::from(&13u32), FE::from(8));
    }

    #[test]
    fn one_sub_1_is_0() {
        assert_eq!(FE::one() - FE::one(), FE::zero());
    }

    #[test]
    fn zero_sub_1_is_order_minus_1() {
        assert_eq!(
            FE::zero() - FE::one(),
            FE::from(&(MERSENNE_31_PRIME_FIELD_ORDER - 1))
        );
    }

    #[test]
    fn neg_1_sub_neg_1_is_0() {
        assert_eq!(-FE::one() - (-FE::one()), FE::zero());
    }

    #[test]
    fn neg_1_sub_0_is_neg_1() {
        assert_eq!(-FE::one() - FE::zero(), -FE::one());
    }

    #[test]
    fn mul_neutral_element() {
        assert_eq!(FE::one() * FE::from(&2u32), FE::from(&2u32));
    }

    #[test]
    fn mul_2_3_is_6() {
        assert_eq!(FE::from(&2u32) * FE::from(&3u32), FE::from(&6u32));
    }

    #[test]
    fn mul_order_neg_1() {
        assert_eq!(
            FE::from(MERSENNE_31_PRIME_FIELD_ORDER as u64 - 1)
                * FE::from(MERSENNE_31_PRIME_FIELD_ORDER as u64 - 1),
            FE::one()
        );
    }

    #[test]
    fn pow_p_neg_1() {
        assert_eq!(
            FE::pow(&FE::from(&2u32), MERSENNE_31_PRIME_FIELD_ORDER - 1),
            FE::one()
        )
    }

    #[test]
    fn inv_0_error() {
        let result = FE::inv(&FE::zero());
        assert!(matches!(result, Err(FieldError::InvZeroError)));
    }

    #[test]
    fn inv_2() {
        let result = FE::inv(&FE::from(&2u32)).unwrap();
        // sage: 1 / F(2) = 1073741824
        assert_eq!(result, FE::from(1073741824));
    }

    #[test]
    fn pow_2_3() {
        assert_eq!(FE::pow(&FE::from(&2u32), 3u64), FE::from(8));
    }

    #[test]
    fn div_1() {
        assert_eq!(FE::from(&2u32) / FE::from(&1u32), FE::from(&2u32));
    }

    #[test]
    fn div_4_2() {
        assert_eq!(FE::from(&4u32) / FE::from(&2u32), FE::from(&2u32));
    }

    #[test]
    fn div_4_3() {
        // sage: F(4) / F(3) = 1431655766
        assert_eq!(FE::from(&4u32) / FE::from(&3u32), FE::from(1431655766));
    }

    #[test]
    fn two_plus_its_additive_inv_is_0() {
        assert_eq!(FE::from(&2u32) + (-FE::from(&2u32)), FE::zero());
    }

    #[test]
    fn from_u64_test() {
        assert_eq!(FE::from(1u64), FE::one());
    }

    #[test]
    fn creating_a_field_element_from_its_representative_returns_the_same_element_1() {
        let change: u32 = MERSENNE_31_PRIME_FIELD_ORDER + 1;
        let f1 = FE::from(&change);
        let f2 = FE::from(&FE::representative(&f1));
        assert_eq!(f1, f2);
    }

    #[test]
    fn creating_a_field_element_from_its_representative_returns_the_same_element_2() {
        let change: u32 = MERSENNE_31_PRIME_FIELD_ORDER + 8;
        let f1 = FE::from(&change);
        let f2 = FE::from(&FE::representative(&f1));
        assert_eq!(f1, f2);
    }

    #[test]
    fn from_base_type_test() {
        assert_eq!(FE::from(&1u32), FE::one());
    }

    #[cfg(feature = "std")]
    #[test]
    fn to_hex_test() {
        let num = FE::from_hex("B").unwrap();
        assert_eq!(FE::to_hex(&num), "B");
    }

    #[test]
    fn double_equals_add_itself() {
        let a = FE::from(1234);
        assert_eq!(a + a, a.double())
    }

    #[test]
    fn two_square_minus_one_is_correct() {
        let a = FE::from(2147483650);
        assert_eq!(
            FE::from(&F::two_square_minus_one(a.value())),
            a.square().double() - FE::one()
        )
    }

    #[test]
    fn two_square_zero_minus_one_is_minus_one() {
        let a = FE::from(0);
        assert_eq!(
            FE::from(&F::two_square_minus_one(a.value())),
            a.square().double() - FE::one()
        )
    }

    #[test]
    fn two_square_p_minus_one_is_minus_one() {
        let a = FE::from(&MERSENNE_31_PRIME_FIELD_ORDER);
        assert_eq!(
            FE::from(&F::two_square_minus_one(a.value())),
            a.square().double() - FE::one()
        )
    }

    #[test]
    fn mul_by_inv() {
        let x = 3476715743_u32;
        assert_eq!(FE::from(&x).inv().unwrap() * FE::from(&x), FE::one());
    }
}
