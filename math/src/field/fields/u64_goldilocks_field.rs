use core::fmt::{self, Display};

#[cfg(feature = "lambdaworks-serde-binary")]
use crate::traits::ByteConversion;
use crate::{
    errors::CreationError,
    field::{
        element::FieldElement,
        errors::FieldError,
        extensions::quadratic::{HasQuadraticNonResidue, QuadraticExtensionField},
        traits::{IsField, IsPrimeField},
    },
};

/// Goldilocks Prime Field F_p where p = 2^64 - 2^32 + 1;
#[derive(Debug, Clone, Copy, Hash, PartialOrd, Ord, PartialEq, Eq)]
pub struct Goldilocks64Field;

impl Goldilocks64Field {
    pub const ORDER: u64 = 0xFFFF_FFFF_0000_0001;
    // Two's complement of `ORDER` i.e. `2^64 - ORDER = 2^32 - 1`
    pub const NEG_ORDER: u64 = Self::ORDER.wrapping_neg();
}

#[cfg(feature = "lambdaworks-serde-binary")]
impl ByteConversion for u64 {
    #[cfg(feature = "alloc")]
    fn to_bytes_be(&self) -> alloc::vec::Vec<u8> {
        self.to_be_bytes().to_vec()
    }

    #[cfg(feature = "alloc")]
    fn to_bytes_le(&self) -> alloc::vec::Vec<u8> {
        self.to_le_bytes().to_vec()
    }

    fn from_bytes_be(_bytes: &[u8]) -> Result<Self, crate::errors::ByteConversionError>
    where
        Self: Sized,
    {
        Ok(u64::from_be_bytes(_bytes.try_into().expect(
            "Tried to create a u64 from the byte array but its size was not 8",
        )))
    }

    fn from_bytes_le(_bytes: &[u8]) -> Result<Self, crate::errors::ByteConversionError>
    where
        Self: Sized,
    {
        Ok(u64::from_le_bytes(_bytes.try_into().expect(
            "Tried to create a u64 from the byte array but its size was not 8",
        )))
    }
}

//NOTE: This implementation was inspired by and borrows from the work done by the Plonky3 team
//https://github.com/Plonky3/Plonky3/blob/main/goldilocks/src/lib.rs
// Thank you for pushing this technology forward.
impl IsField for Goldilocks64Field {
    type BaseType = u64;

    fn add(a: &u64, b: &u64) -> u64 {
        let (sum, over) = a.overflowing_add(*b);
        let (mut sum, over) = sum.overflowing_add(u64::from(over) * Self::NEG_ORDER);
        if over {
            sum += Self::NEG_ORDER
        }
        Self::representative(&sum)
    }

    fn mul(a: &u64, b: &u64) -> u64 {
        Self::representative(&reduce_128(u128::from(*a) * u128::from(*b)))
    }

    fn sub(a: &u64, b: &u64) -> u64 {
        let (diff, under) = a.overflowing_sub(*b);
        let (mut diff, under) = diff.overflowing_sub(u64::from(under) * Self::NEG_ORDER);
        if under {
            diff -= Self::NEG_ORDER;
        }
        Self::representative(&diff)
    }

    fn neg(a: &u64) -> u64 {
        Self::sub(&Self::ORDER, &Self::representative(a))
    }

    /// Returns the multiplicative inverse of `a`.
    fn inv(a: &u64) -> Result<u64, FieldError> {
        if *a == Self::zero() || *a == Self::ORDER {
            return Err(FieldError::InvZeroError);
        }

        // a^11
        let t2 = Self::mul(&Self::square(a), a);

        // a^111
        let t3 = Self::mul(&Self::square(&t2), a);

        // compute base^111111 (6 ones) by repeatedly squaring t3 3 times and multiplying by t3
        let t6 = exp_acc::<3>(&t3, &t3);
        let t60 = Self::square(&t6);
        let t7 = Self::mul(&t60, a);

        // compute base^111111111111 (12 ones)
        // repeatedly square t6 6 times and multiply by t6
        let t12 = exp_acc::<5>(&t60, &t6);

        // compute base^111111111111111111111111 (24 ones)
        // repeatedly square t12 12 times and multiply by t12
        let t24 = exp_acc::<12>(&t12, &t12);

        // compute base^1111111111111111111111111111111 (31 ones)
        // repeatedly square t24 6 times and multiply by t6 first. then square t30 and multiply by base
        let t31 = exp_acc::<7>(&t24, &t7);

        // compute base^111111111111111111111111111111101111111111111111111111111111111
        // repeatedly square t31 32 times and multiply by t31
        let t63 = exp_acc::<32>(&t31, &t31);

        Ok(Self::mul(&Self::square(&t63), a))
    }

    /// Returns the division of `a` and `b`.
    fn div(a: &u64, b: &u64) -> u64 {
        let b_inv = Self::inv(b).unwrap();
        Self::mul(a, &b_inv)
    }

    /// Returns a boolean indicating whether `a` and `b` are equal or not.
    fn eq(a: &u64, b: &u64) -> bool {
        Self::representative(a) == Self::representative(b)
    }

    /// Returns the additive neutral element.
    fn zero() -> u64 {
        0u64
    }

    /// Returns the multiplicative neutral element.
    fn one() -> u64 {
        1u64
    }

    /// Returns the element `x * 1` where 1 is the multiplicative neutral element.
    fn from_u64(x: u64) -> u64 {
        Self::representative(&x)
    }

    /// Takes as input an element of BaseType and returns the internal representation
    /// of that element in the field.
    fn from_base_type(x: u64) -> u64 {
        Self::representative(&x)
    }
}

impl IsPrimeField for Goldilocks64Field {
    type RepresentativeType = u64;

    fn representative(x: &u64) -> u64 {
        let mut u = *x;
        if u >= Self::ORDER {
            u -= Self::ORDER;
        }
        u
    }

    fn field_bit_size() -> usize {
        ((self::Goldilocks64Field::ORDER - 1).ilog2() + 1) as usize
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
        u64::from_str_radix(hex_string, 16).map_err(|_| CreationError::InvalidHexString)
    }

    #[cfg(feature = "std")]
    fn to_hex(x: &u64) -> String {
        format!("{:X}", x)
    }
}

#[inline(always)]
fn reduce_128(x: u128) -> u64 {
    //possibly split apart into separate function to ensure inline
    let (x_lo, x_hi) = (x as u64, (x >> 64) as u64);
    let x_hi_hi = x_hi >> 32;
    let x_hi_lo = x_hi & Goldilocks64Field::NEG_ORDER;

    let (mut t0, borrow) = x_lo.overflowing_sub(x_hi_hi);
    if borrow {
        t0 -= Goldilocks64Field::NEG_ORDER // Cannot underflow
    }

    let t1 = x_hi_lo * Goldilocks64Field::NEG_ORDER;
    let (res_wrapped, carry) = t0.overflowing_add(t1);
    // Below cannot overflow unless the assumption if x + y < 2**64 + ORDER is incorrect.
    res_wrapped + Goldilocks64Field::NEG_ORDER * u64::from(carry)
}

#[inline(always)]
fn exp_acc<const N: usize>(base: &u64, tail: &u64) -> u64 {
    Goldilocks64Field::mul(&exp_power_of_2::<N>(base), tail)
}

#[must_use]
fn exp_power_of_2<const POWER_LOG: usize>(base: &u64) -> u64 {
    let mut res = *base;
    for _ in 0..POWER_LOG {
        res = Goldilocks64Field::square(&res);
    }
    res
}

pub type Goldilocks64ExtensionField = QuadraticExtensionField<Goldilocks64Field, Goldilocks64Field>;

impl HasQuadraticNonResidue<Goldilocks64Field> for Goldilocks64Field {
    // Verifiable in Sage with
    // `R.<x> = GF(p)[]; assert (x^2 - 7).is_irreducible()`
    fn residue() -> FieldElement<Goldilocks64Field> {
        FieldElement::from(Goldilocks64Field::from_u64(7u64))
    }
}

impl Display for FieldElement<Goldilocks64Field> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:x}", self.representative())?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    type F = Goldilocks64Field;

    // Over the Goldilocks field, the following set of equations hold
    // p               = 0
    // 2^64 - 2^32 + 1 = 0
    // 2^64            = 2^32 - 1
    #[test]
    fn from_hex_for_b_is_11() {
        assert_eq!(F::from_hex("B").unwrap(), 11);
    }

    #[test]
    fn from_hex_for_0x1_a_is_26() {
        assert_eq!(F::from_hex("0x1a").unwrap(), 26);
    }

    #[test]
    fn bit_size_of_field_is_64() {
        assert_eq!(
            <F as crate::field::traits::IsPrimeField>::field_bit_size(),
            64
        );
    }

    #[test]
    fn one_plus_one_is_two() {
        let a = F::one();
        let b = F::one();
        let c = F::add(&a, &b);
        assert_eq!(c, 2u64);
    }

    #[test]
    fn neg_one_plus_one_is_zero() {
        let a = F::neg(&F::one());
        let b = F::one();
        let c = F::add(&a, &b);
        assert_eq!(c, F::zero());
    }

    #[test]
    fn neg_one_plus_two_is_one() {
        let a = F::neg(&F::one());
        let b = F::from_base_type(2u64);
        let c = F::add(&a, &b);
        assert_eq!(c, F::one());
    }

    #[test]
    fn max_order_plus_one_is_zero() {
        let a = F::from_base_type(F::ORDER - 1);
        let b = F::one();
        let c = F::add(&a, &b);
        assert_eq!(c, F::zero());
    }

    #[test]
    fn comparing_13_and_13_are_equal() {
        let a = F::from_base_type(13);
        let b = F::from_base_type(13);
        assert_eq!(a, b);
    }

    #[test]
    fn comparing_13_and_8_they_are_not_equal() {
        let a = F::from_base_type(13);
        let b = F::from_base_type(8);
        assert_ne!(a, b);
    }

    #[test]
    fn one_sub_one_is_zero() {
        let a = F::one();
        let b = F::one();
        let c = F::sub(&a, &b);
        assert_eq!(c, F::zero());
    }

    #[test]
    fn zero_sub_one_is_order_minus_1() {
        let a = F::zero();
        let b = F::one();
        let c = F::sub(&a, &b);
        assert_eq!(c, F::ORDER - 1);
    }

    #[test]
    fn neg_one_sub_neg_one_is_zero() {
        let a = F::neg(&F::one());
        let b = F::neg(&F::one());
        let c = F::sub(&a, &b);
        assert_eq!(c, F::zero());
    }

    #[test]
    fn neg_one_sub_one_is_neg_one() {
        let a = F::neg(&F::one());
        let b = F::zero();
        let c = F::sub(&a, &b);
        assert_eq!(c, F::neg(&F::one()));
    }

    #[test]
    fn mul_neutral_element() {
        let a = F::from_base_type(1);
        let b = F::from_base_type(2);
        let c = F::mul(&a, &b);
        assert_eq!(c, F::from_base_type(2));
    }

    #[test]
    fn mul_two_three_is_six() {
        let a = F::from_base_type(2);
        let b = F::from_base_type(3);
        assert_eq!(a * b, F::from_base_type(6));
    }

    #[test]
    fn mul_order_neg_one() {
        let a = F::from_base_type(F::ORDER - 1);
        let b = F::from_base_type(F::ORDER - 1);
        let c = F::mul(&a, &b);
        assert_eq!(c, F::from_base_type(1));
    }

    #[test]
    fn pow_p_neg_one() {
        assert_eq!(F::pow(&F::from_base_type(2), F::ORDER - 1), F::one())
    }

    #[test]
    fn inv_zero_error() {
        let result = F::inv(&F::zero());
        assert!(matches!(result, Err(FieldError::InvZeroError)));
    }

    #[test]
    fn inv_two() {
        let result = F::inv(&F::from_base_type(2u64)).unwrap();
        // sage: 1 / F(2) = 9223372034707292161
        assert_eq!(result, 9223372034707292161);
    }

    #[test]
    fn pow_two_three() {
        assert_eq!(F::pow(&F::from_base_type(2), 3_u64), 8)
    }

    #[test]
    fn div_one() {
        assert_eq!(F::div(&F::from_base_type(2), &F::from_base_type(1)), 2)
    }

    #[test]
    fn div_4_2() {
        assert_eq!(F::div(&F::from_base_type(4), &F::from_base_type(2)), 2)
    }

    // 1431655766
    #[test]
    fn div_4_3() {
        // sage: F(4) / F(3) = 12297829379609722882
        assert_eq!(
            F::div(&F::from_base_type(4), &F::from_base_type(3)),
            12297829379609722882
        )
    }

    #[test]
    fn two_plus_its_additive_inv_is_0() {
        let two = F::from_base_type(2);

        assert_eq!(F::add(&two, &F::neg(&two)), F::zero())
    }

    #[test]
    fn from_u64_test() {
        let num = F::from_u64(1u64);
        assert_eq!(num, F::one());
    }

    #[test]
    fn from_u64_zero_test() {
        let num = F::from_u64(0);
        assert_eq!(num, F::zero());
    }

    #[test]
    fn from_u64_max_test() {
        let num = F::from_u64(u64::MAX);
        assert_eq!(num, u32::MAX as u64 - 1);
    }

    #[test]
    fn from_u64_order_test() {
        let num = F::from_u64(F::ORDER);
        assert_eq!(num, F::zero());
    }

    #[test]
    fn creating_a_field_element_from_its_representative_returns_the_same_element_1() {
        let change = 1;
        let f1 = F::from_base_type(F::ORDER + change);
        let f2 = F::from_base_type(F::representative(&f1));
        assert_eq!(f1, f2);
    }

    #[test]
    fn reduct_128() {
        let x = u128::MAX;
        let y = reduce_128(x);
        // The following equalitiy sequence holds, modulo p = 2^64 - 2^32 + 1
        // 2^128 - 1 = (2^64 - 1) * (2^64 + 1)
        //           = (2^32 - 1 - 1) * (2^32 - 1 + 1)
        //           = (2^32 - 2) * (2^32)
        //           = 2^64 - 2 * 2^32
        //           = 2^64 - 2^33
        //           = 2^32 - 1 - 2^33
        //           = - 2^32 - 1
        let expected_result = F::neg(&F::add(&F::from_base_type(2_u64.pow(32)), &F::one()));
        assert_eq!(y, expected_result);
    }

    #[test]
    fn u64_max_as_representative_less_than_u32_max_sub_1() {
        let f = F::from_base_type(u64::MAX);
        assert_eq!(F::representative(&f), u32::MAX as u64 - 1)
    }

    #[test]
    fn creating_a_field_element_from_its_representative_returns_the_same_element_2() {
        let change = 8;
        let f1 = F::from_base_type(F::ORDER + change);
        let f2 = F::from_base_type(F::representative(&f1));
        assert_eq!(f1, f2);
    }

    #[test]
    fn from_base_type_test() {
        let b = F::from_base_type(1u64);
        assert_eq!(b, F::one());
    }

    #[cfg(feature = "std")]
    #[test]
    fn to_hex_test() {
        let num = F::from_hex("B").unwrap();
        assert_eq!(F::to_hex(&num), "B");
    }
}
