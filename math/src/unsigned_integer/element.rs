use core::cmp::Ordering;
use core::convert::From;
use core::ops::{
    Add, BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign, Mul, Shl, Shr, ShrAssign,
    Sub,
};

#[cfg(feature = "proptest")]
use proptest::{
    arbitrary::Arbitrary,
    prelude::any,
    strategy::{SBoxedStrategy, Strategy},
};

use crate::errors::ByteConversionError;
use crate::errors::CreationError;
#[cfg(feature = "alloc")]
use crate::traits::AsBytes;
use crate::traits::ByteConversion;
use crate::unsigned_integer::traits::IsUnsignedInteger;

use core::fmt::{self, Debug, Display};

pub type U384 = UnsignedInteger<6>;
pub type U256 = UnsignedInteger<4>;
pub type U128 = UnsignedInteger<2>;
pub type U64 = UnsignedInteger<1>;

/// A big unsigned integer in base 2^{64} represented
/// as fixed-size array `limbs` of `u64` components.
/// The most significant bit is in the left-most position.
/// That is, the array `[a_n, ..., a_0]` represents the
/// integer 2^{64 * n} * a_n + ... + 2^{64} * a_1 + a_0.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub struct UnsignedInteger<const NUM_LIMBS: usize> {
    pub limbs: [u64; NUM_LIMBS],
}

impl<const NUM_LIMBS: usize> Default for UnsignedInteger<NUM_LIMBS> {
    fn default() -> Self {
        Self {
            limbs: [0; NUM_LIMBS],
        }
    }
}

// NOTE: manually implementing `PartialOrd` may seem unorthodox, but the
// derived implementation had terrible performance.
#[allow(clippy::non_canonical_partial_ord_impl)]
impl<const NUM_LIMBS: usize> PartialOrd for UnsignedInteger<NUM_LIMBS> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let mut i = 0;
        while i < NUM_LIMBS {
            if self.limbs[i] != other.limbs[i] {
                return Some(self.limbs[i].cmp(&other.limbs[i]));
            }
            i += 1;
        }
        Some(Ordering::Equal)
    }
}

// NOTE: because we implemented `PartialOrd`, clippy asks us to implement
// this manually too.
impl<const NUM_LIMBS: usize> Ord for UnsignedInteger<NUM_LIMBS> {
    fn cmp(&self, other: &Self) -> Ordering {
        let mut i = 0;
        while i < NUM_LIMBS {
            if self.limbs[i] != other.limbs[i] {
                return self.limbs[i].cmp(&other.limbs[i]);
            }
            i += 1;
        }
        Ordering::Equal
    }
}

impl<const NUM_LIMBS: usize> From<u128> for UnsignedInteger<NUM_LIMBS> {
    fn from(value: u128) -> Self {
        let mut limbs = [0u64; NUM_LIMBS];
        limbs[NUM_LIMBS - 1] = value as u64;
        limbs[NUM_LIMBS - 2] = (value >> 64) as u64;
        UnsignedInteger { limbs }
    }
}

impl<const NUM_LIMBS: usize> From<u64> for UnsignedInteger<NUM_LIMBS> {
    fn from(value: u64) -> Self {
        Self::from_u64(value)
    }
}

impl<const NUM_LIMBS: usize> From<u16> for UnsignedInteger<NUM_LIMBS> {
    fn from(value: u16) -> Self {
        let mut limbs = [0u64; NUM_LIMBS];
        limbs[NUM_LIMBS - 1] = value as u64;
        UnsignedInteger { limbs }
    }
}

impl<const NUM_LIMBS: usize> From<&str> for UnsignedInteger<NUM_LIMBS> {
    fn from(hex_str: &str) -> Self {
        Self::from_hex_unchecked(hex_str)
    }
}

impl<const NUM_LIMBS: usize> Display for UnsignedInteger<NUM_LIMBS> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut limbs_iterator = self.limbs.iter().skip_while(|limb| **limb == 0).peekable();

        if limbs_iterator.peek().is_none() {
            write!(f, "0x0")?;
        } else {
            write!(f, "0x")?;
            if let Some(most_significant_limb) = limbs_iterator.next() {
                write!(f, "{:x}", most_significant_limb)?;
            }

            for limb in limbs_iterator {
                write!(f, "{:016x}", limb)?;
            }
        }

        Ok(())
    }
}

// impl Add

impl<const NUM_LIMBS: usize> Add<&UnsignedInteger<NUM_LIMBS>> for &UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;

    fn add(self, other: &UnsignedInteger<NUM_LIMBS>) -> UnsignedInteger<NUM_LIMBS> {
        let (result, overflow) = UnsignedInteger::add(self, other);
        debug_assert!(!overflow, "UnsignedInteger addition overflow.");
        result
    }
}

impl<const NUM_LIMBS: usize> Add<UnsignedInteger<NUM_LIMBS>> for UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;

    fn add(self, other: UnsignedInteger<NUM_LIMBS>) -> UnsignedInteger<NUM_LIMBS> {
        &self + &other
    }
}

impl<const NUM_LIMBS: usize> Add<&UnsignedInteger<NUM_LIMBS>> for UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;

    fn add(self, other: &Self) -> Self {
        &self + other
    }
}

impl<const NUM_LIMBS: usize> Add<UnsignedInteger<NUM_LIMBS>> for &UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;

    fn add(self, other: UnsignedInteger<NUM_LIMBS>) -> UnsignedInteger<NUM_LIMBS> {
        self + &other
    }
}

// impl Sub

impl<const NUM_LIMBS: usize> Sub<&UnsignedInteger<NUM_LIMBS>> for &UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;

    fn sub(self, other: &UnsignedInteger<NUM_LIMBS>) -> UnsignedInteger<NUM_LIMBS> {
        let (result, overflow) = UnsignedInteger::sub(self, other);
        debug_assert!(!overflow, "UnsignedInteger subtraction overflow.");
        result
    }
}

impl<const NUM_LIMBS: usize> Sub<UnsignedInteger<NUM_LIMBS>> for UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;

    fn sub(self, other: UnsignedInteger<NUM_LIMBS>) -> UnsignedInteger<NUM_LIMBS> {
        &self - &other
    }
}

impl<const NUM_LIMBS: usize> Sub<&UnsignedInteger<NUM_LIMBS>> for UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;

    fn sub(self, other: &Self) -> Self {
        &self - other
    }
}

impl<const NUM_LIMBS: usize> Sub<UnsignedInteger<NUM_LIMBS>> for &UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;
    #[inline(always)]
    fn sub(self, other: UnsignedInteger<NUM_LIMBS>) -> UnsignedInteger<NUM_LIMBS> {
        self - &other
    }
}

/// Multi-precision multiplication.
/// Algorithm 14.12 of "Handbook of Applied Cryptography" (https://cacr.uwaterloo.ca/hac/)
impl<const NUM_LIMBS: usize> Mul<&UnsignedInteger<NUM_LIMBS>> for &UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;

    #[inline(always)]
    fn mul(self, other: &UnsignedInteger<NUM_LIMBS>) -> UnsignedInteger<NUM_LIMBS> {
        let (mut n, mut t) = (0, 0);
        for i in (0..NUM_LIMBS).rev() {
            if self.limbs[i] != 0u64 {
                n = NUM_LIMBS - 1 - i;
            }
            if other.limbs[i] != 0u64 {
                t = NUM_LIMBS - 1 - i;
            }
        }
        debug_assert!(
            n + t < NUM_LIMBS,
            "UnsignedInteger multiplication overflow."
        );

        // 1.
        let mut limbs = [0u64; NUM_LIMBS];
        // 2.
        let mut carry = 0u128;
        for i in 0..=t {
            // 2.2
            for j in 0..=n {
                let uv = (limbs[NUM_LIMBS - 1 - (i + j)] as u128)
                    + (self.limbs[NUM_LIMBS - 1 - j] as u128)
                        * (other.limbs[NUM_LIMBS - 1 - i] as u128)
                    + carry;
                carry = uv >> 64;
                limbs[NUM_LIMBS - 1 - (i + j)] = uv as u64;
            }
            if i + n + 1 < NUM_LIMBS {
                // 2.3
                limbs[NUM_LIMBS - 1 - (i + n + 1)] = carry as u64;
                carry = 0;
            }
        }
        assert_eq!(carry, 0, "UnsignedInteger multiplication overflow.");
        // 3.
        Self::Output { limbs }
    }
}

impl<const NUM_LIMBS: usize> Mul<UnsignedInteger<NUM_LIMBS>> for UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;
    #[inline(always)]
    fn mul(self, other: UnsignedInteger<NUM_LIMBS>) -> UnsignedInteger<NUM_LIMBS> {
        &self * &other
    }
}

impl<const NUM_LIMBS: usize> Mul<&UnsignedInteger<NUM_LIMBS>> for UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;
    #[inline(always)]
    fn mul(self, other: &Self) -> Self {
        &self * other
    }
}

impl<const NUM_LIMBS: usize> Mul<UnsignedInteger<NUM_LIMBS>> for &UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;
    #[inline(always)]
    fn mul(self, other: UnsignedInteger<NUM_LIMBS>) -> UnsignedInteger<NUM_LIMBS> {
        self * &other
    }
}

impl<const NUM_LIMBS: usize> Shl<usize> for &UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;
    #[inline(always)]
    fn shl(self, times: usize) -> UnsignedInteger<NUM_LIMBS> {
        self.const_shl(times)
    }
}

impl<const NUM_LIMBS: usize> Shl<usize> for UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;
    #[inline(always)]
    fn shl(self, times: usize) -> UnsignedInteger<NUM_LIMBS> {
        &self << times
    }
}

// impl Shr

impl<const NUM_LIMBS: usize> Shr<usize> for &UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;
    #[inline(always)]
    fn shr(self, times: usize) -> UnsignedInteger<NUM_LIMBS> {
        self.const_shr(times)
    }
}

impl<const NUM_LIMBS: usize> Shr<usize> for UnsignedInteger<NUM_LIMBS> {
    type Output = UnsignedInteger<NUM_LIMBS>;
    #[inline(always)]
    fn shr(self, times: usize) -> UnsignedInteger<NUM_LIMBS> {
        &self >> times
    }
}

impl<const NUM_LIMBS: usize> ShrAssign<usize> for UnsignedInteger<NUM_LIMBS> {
    fn shr_assign(&mut self, times: usize) {
        debug_assert!(
            times < 64 * NUM_LIMBS,
            "UnsignedInteger shift left overflows."
        );

        let (a, b) = (times / 64, times % 64);

        if b == 0 {
            self.limbs.copy_within(..NUM_LIMBS - a, a);
        } else {
            for i in (a + 1..NUM_LIMBS).rev() {
                self.limbs[i] = (self.limbs[i - a] >> b) | (self.limbs[i - a - 1] << (64 - b));
            }
            self.limbs[a] = self.limbs[0] >> b;
        }

        for limb in self.limbs.iter_mut().take(a) {
            *limb = 0;
        }
    }
}

/// Impl BitAnd
impl<const NUM_LIMBS: usize> BitAnd for UnsignedInteger<NUM_LIMBS> {
    type Output = Self;

    #[inline(always)]
    fn bitand(self, rhs: Self) -> Self::Output {
        let mut result = self;
        result &= rhs;
        result
    }
}

impl<const NUM_LIMBS: usize> BitAndAssign for UnsignedInteger<NUM_LIMBS> {
    fn bitand_assign(&mut self, rhs: Self) {
        for (a_i, b_i) in self.limbs.iter_mut().zip(rhs.limbs.iter()) {
            *a_i &= b_i;
        }
    }
}

/// Impl BitOr
impl<const NUM_LIMBS: usize> BitOr for UnsignedInteger<NUM_LIMBS> {
    type Output = Self;

    #[inline(always)]
    fn bitor(self, rhs: Self) -> Self::Output {
        let mut result = self;
        result |= rhs;
        result
    }
}

impl<const NUM_LIMBS: usize> BitOrAssign for UnsignedInteger<NUM_LIMBS> {
    #[inline(always)]
    fn bitor_assign(&mut self, rhs: Self) {
        for (a_i, b_i) in self.limbs.iter_mut().zip(rhs.limbs.iter()) {
            *a_i |= b_i;
        }
    }
}

/// Impl BitXor
impl<const NUM_LIMBS: usize> BitXor for UnsignedInteger<NUM_LIMBS> {
    type Output = Self;

    #[inline(always)]
    fn bitxor(self, rhs: Self) -> Self::Output {
        let mut result = self;
        result ^= rhs;
        result
    }
}

impl<const NUM_LIMBS: usize> BitXorAssign for UnsignedInteger<NUM_LIMBS> {
    #[inline(always)]
    fn bitxor_assign(&mut self, rhs: Self) {
        for (a_i, b_i) in self.limbs.iter_mut().zip(rhs.limbs.iter()) {
            *a_i ^= b_i;
        }
    }
}

impl<const NUM_LIMBS: usize> UnsignedInteger<NUM_LIMBS> {
    pub const fn from_limbs(limbs: [u64; NUM_LIMBS]) -> Self {
        Self { limbs }
    }

    #[inline(always)]
    pub const fn from_u64(value: u64) -> Self {
        let mut limbs = [0u64; NUM_LIMBS];
        limbs[NUM_LIMBS - 1] = value;
        UnsignedInteger { limbs }
    }

    #[inline(always)]
    pub const fn from_u128(value: u128) -> Self {
        let mut limbs = [0u64; NUM_LIMBS];
        limbs[NUM_LIMBS - 1] = value as u64;
        limbs[NUM_LIMBS - 2] = (value >> 64) as u64;
        UnsignedInteger { limbs }
    }

    #[inline(always)]
    const fn is_hex_string(string: &str) -> bool {
        let len: usize = string.len();
        let bytes = string.as_bytes();
        let mut i = 0;

        while i < len {
            match bytes[i] {
                b'0'..=b'9' => (),
                b'a'..=b'f' => (),
                b'A'..=b'F' => (),
                _ => return false,
            }
            i += 1;
        }

        true
    }

    /// Creates an `UnsignedInteger` from a hexstring. It can contain `0x` or not.
    /// Returns an `CreationError::InvalidHexString`if the value is not a hexstring.
    /// Returns a `CreationError::EmptyString` if the input string is empty.
    /// Returns a `CreationError::HexStringIsTooBig` if the the input hex string is bigger
    /// than the maximum amount of characters for this element.
    pub fn from_hex(value: &str) -> Result<Self, CreationError> {
        let mut string = value;
        let mut char_iterator = value.chars();
        if string.len() > 2
            && char_iterator.next().unwrap() == '0'
            && char_iterator.next().unwrap() == 'x'
        {
            string = &string[2..];
        }
        if string.is_empty() {
            return Err(CreationError::EmptyString);
        }
        if !Self::is_hex_string(string) {
            return Err(CreationError::InvalidHexString);
        }

        // Limbs are of 64 bits - 8 bytes
        // We have 16 nibbles per bytes
        let max_amount_of_hex_chars = NUM_LIMBS * 16;
        if string.len() > max_amount_of_hex_chars {
            return Err(CreationError::HexStringIsTooBig);
        }

        Ok(Self::from_hex_unchecked(string))
    }

    /// Creates an `UnsignedInteger` from a hexstring
    /// # Panics
    /// Panics if value is not a hexstring. It can contain `0x` or not.
    pub const fn from_hex_unchecked(value: &str) -> Self {
        let mut result = [0u64; NUM_LIMBS];
        let mut limb = 0;
        let mut limb_index = NUM_LIMBS - 1;
        let mut shift = 0;

        let value_bytes = value.as_bytes();

        // Remove "0x" if it's at the beginning of the string
        let mut i = 0;
        if value_bytes.len() > 2 && value_bytes[0] == b'0' && value_bytes[1] == b'x' {
            i = 2;
        }

        let mut j = value_bytes.len();
        while j > i {
            j -= 1;
            limb |= match value_bytes[j] {
                c @ b'0'..=b'9' => (c as u64 - b'0' as u64) << shift,
                c @ b'a'..=b'f' => (c as u64 - b'a' as u64 + 10) << shift,
                c @ b'A'..=b'F' => (c as u64 - b'A' as u64 + 10) << shift,
                _ => panic!("Malformed hex expression."),
            };
            shift += 4;
            if shift == 64 && limb_index > 0 {
                result[limb_index] = limb;
                limb = 0;
                limb_index -= 1;
                shift = 0;
            }
        }

        result[limb_index] = limb;
        UnsignedInteger { limbs: result }
    }

    /// Creates a hexstring from a `FieldElement` without `0x`.
    #[cfg(feature = "std")]
    pub fn to_hex(&self) -> String {
        let mut hex_string = String::new();
        for &limb in self.limbs.iter() {
            hex_string.push_str(&format!("{:016X}", limb));
        }
        hex_string.trim_start_matches('0').to_string()
    }

    pub const fn const_ne(a: &UnsignedInteger<NUM_LIMBS>, b: &UnsignedInteger<NUM_LIMBS>) -> bool {
        let mut i = 0;
        while i < NUM_LIMBS {
            if a.limbs[i] != b.limbs[i] {
                return true;
            }
            i += 1;
        }
        false
    }

    pub const fn const_le(a: &UnsignedInteger<NUM_LIMBS>, b: &UnsignedInteger<NUM_LIMBS>) -> bool {
        let mut i = 0;
        while i < NUM_LIMBS {
            if a.limbs[i] < b.limbs[i] {
                return true;
            } else if a.limbs[i] > b.limbs[i] {
                return false;
            }
            i += 1;
        }
        true
    }

    pub const fn const_shl(self, times: usize) -> Self {
        debug_assert!(
            times < 64 * NUM_LIMBS,
            "UnsignedInteger shift left overflows."
        );
        let mut limbs = [0u64; NUM_LIMBS];
        let (a, b) = (times / 64, times % 64);

        if b == 0 {
            let mut i = 0;
            while i < NUM_LIMBS - a {
                limbs[i] = self.limbs[a + i];
                i += 1;
            }
        } else {
            limbs[NUM_LIMBS - 1 - a] = self.limbs[NUM_LIMBS - 1] << b;
            let mut i = a + 1;
            while i < NUM_LIMBS {
                limbs[NUM_LIMBS - 1 - i] = (self.limbs[NUM_LIMBS - 1 - i + a] << b)
                    | (self.limbs[NUM_LIMBS - i + a] >> (64 - b));
                i += 1;
            }
        }
        Self { limbs }
    }

    pub const fn const_shr(self, times: usize) -> UnsignedInteger<NUM_LIMBS> {
        debug_assert!(
            times < 64 * NUM_LIMBS,
            "UnsignedInteger shift right overflows."
        );

        let mut limbs = [0u64; NUM_LIMBS];
        let (a, b) = (times / 64, times % 64);

        if b == 0 {
            let mut i = 0;
            while i < NUM_LIMBS - a {
                limbs[a + i] = self.limbs[i];
                i += 1;
            }
        } else {
            limbs[a] = self.limbs[0] >> b;
            let mut i = a + 1;
            while i < NUM_LIMBS {
                limbs[i] = (self.limbs[i - a - 1] << (64 - b)) | (self.limbs[i - a] >> b);
                i += 1;
            }
        }

        Self { limbs }
    }

    pub const fn add(
        a: &UnsignedInteger<NUM_LIMBS>,
        b: &UnsignedInteger<NUM_LIMBS>,
    ) -> (UnsignedInteger<NUM_LIMBS>, bool) {
        let mut limbs = [0u64; NUM_LIMBS];
        let mut carry = 0u64;
        let mut i = NUM_LIMBS;
        while i > 0 {
            let (x, cb) = a.limbs[i - 1].overflowing_add(b.limbs[i - 1]);
            let (x, cc) = x.overflowing_add(carry);
            limbs[i - 1] = x;
            carry = (cb | cc) as u64;
            i -= 1;
        }
        (UnsignedInteger { limbs }, carry > 0)
    }

    pub fn double(a: &UnsignedInteger<NUM_LIMBS>) -> (UnsignedInteger<NUM_LIMBS>, bool) {
        let mut cloned = *a;
        let overflow = cloned.double_in_place();
        (cloned, overflow)
    }

    pub fn double_in_place(&mut self) -> bool {
        let mut msb_of_previous_limb = 0;
        for i in (0..NUM_LIMBS).rev() {
            let limb_ref = &mut self.limbs[i];
            let msb_of_current_limb = *limb_ref >> 63;
            *limb_ref <<= 1;
            *limb_ref |= msb_of_previous_limb;
            msb_of_previous_limb = msb_of_current_limb;
        }
        msb_of_previous_limb != 0
    }

    /// Multi-precision subtraction.
    /// Adapted from Algorithm 14.9 of "Handbook of Applied Cryptography" (https://cacr.uwaterloo.ca/hac/)
    /// Returns the results and a flag that is set if the substraction underflowed
    #[inline(always)]
    pub const fn sub(
        a: &UnsignedInteger<NUM_LIMBS>,
        b: &UnsignedInteger<NUM_LIMBS>,
    ) -> (UnsignedInteger<NUM_LIMBS>, bool) {
        let mut limbs = [0u64; NUM_LIMBS];
        // 1.
        let mut carry = false;
        // 2.
        let mut i: usize = NUM_LIMBS;
        while i > 0 {
            i -= 1;
            let (x, cb) = a.limbs[i].overflowing_sub(b.limbs[i]);
            let (x, cc) = x.overflowing_sub(carry as u64);
            // Casting i128 to u64 drops the most significant bits of i128,
            // which effectively computes residue modulo 2^{64}
            // 2.1
            limbs[i] = x;
            // 2.2
            carry = cb | cc;
        }
        // 3.
        (Self { limbs }, carry)
    }

    /// Multi-precision multiplication.
    /// Adapted from Algorithm 14.12 of "Handbook of Applied Cryptography" (https://cacr.uwaterloo.ca/hac/)
    pub const fn mul(
        a: &UnsignedInteger<NUM_LIMBS>,
        b: &UnsignedInteger<NUM_LIMBS>,
    ) -> (UnsignedInteger<NUM_LIMBS>, UnsignedInteger<NUM_LIMBS>) {
        // 1.
        let mut hi = [0u64; NUM_LIMBS];
        let mut lo = [0u64; NUM_LIMBS];
        // Const functions don't support for loops so we use whiles
        // this is equivalent to:
        // for i in (0..NUM_LIMBS).rev()
        // 2.
        let mut i = NUM_LIMBS;
        while i > 0 {
            i -= 1;
            // 2.1
            let mut carry = 0u128;
            let mut j = NUM_LIMBS;
            // 2.2
            while j > 0 {
                j -= 1;
                let mut k = i + j;
                if k >= NUM_LIMBS - 1 {
                    k -= NUM_LIMBS - 1;
                    let uv = (lo[k] as u128) + (a.limbs[j] as u128) * (b.limbs[i] as u128) + carry;
                    carry = uv >> 64;
                    // Casting u128 to u64 takes modulo 2^{64}
                    lo[k] = uv as u64;
                } else {
                    let uv =
                        (hi[k + 1] as u128) + (a.limbs[j] as u128) * (b.limbs[i] as u128) + carry;
                    carry = uv >> 64;
                    // Casting u128 to u64 takes modulo 2^{64}
                    hi[k + 1] = uv as u64;
                }
            }
            // 2.3
            hi[i] = carry as u64;
        }
        // 3.
        (Self { limbs: hi }, Self { limbs: lo })
    }

    #[inline(always)]
    pub fn square(
        a: &UnsignedInteger<NUM_LIMBS>,
    ) -> (UnsignedInteger<NUM_LIMBS>, UnsignedInteger<NUM_LIMBS>) {
        // NOTE: we use explicit `while` loops in this function because profiling pointed
        // at iterators of the form `(<x>..<y>).rev()` as the main performance bottleneck.

        let mut hi = Self {
            limbs: [0u64; NUM_LIMBS],
        };
        let mut lo = Self {
            limbs: [0u64; NUM_LIMBS],
        };

        // Compute products between a[i] and a[j] when i != j.
        // The variable `index` below is the index of `lo` or
        // `hi` to update
        let mut i = NUM_LIMBS;
        while i > 1 {
            i -= 1;
            let mut c: u128 = 0;
            let mut j = i;
            while j > 0 {
                j -= 1;
                let k = i + j;
                if k >= NUM_LIMBS - 1 {
                    let index = k + 1 - NUM_LIMBS;
                    let cs = lo.limbs[index] as u128 + a.limbs[i] as u128 * a.limbs[j] as u128 + c;
                    c = cs >> 64;
                    lo.limbs[index] = cs as u64;
                } else {
                    let index = k + 1;
                    let cs = hi.limbs[index] as u128 + a.limbs[i] as u128 * a.limbs[j] as u128 + c;
                    c = cs >> 64;
                    hi.limbs[index] = cs as u64;
                }
            }
            hi.limbs[i] = c as u64;
        }

        // All these terms should appear twice each,
        // so we have to multiply what we got so far by two.
        let carry = lo.limbs[0] >> 63;
        lo = lo << 1;
        hi = hi << 1;
        hi.limbs[NUM_LIMBS - 1] |= carry;

        // Add the only remaning terms, which are the squares a[i] * a[i].
        // The variable `index` below is the index of `lo` or
        // `hi` to update
        let mut c = 0;
        let mut i = NUM_LIMBS;
        while i > 0 {
            i -= 1;
            if NUM_LIMBS - 1 <= i * 2 {
                let index = 2 * i + 1 - NUM_LIMBS;
                let cs = lo.limbs[index] as u128 + a.limbs[i] as u128 * a.limbs[i] as u128 + c;
                c = cs >> 64;
                lo.limbs[index] = cs as u64;
            } else {
                let index = 2 * i + 1;
                let cs = hi.limbs[index] as u128 + a.limbs[i] as u128 * a.limbs[i] as u128 + c;
                c = cs >> 64;
                hi.limbs[index] = cs as u64;
            }
            if NUM_LIMBS - 1 < i * 2 {
                let index = 2 * i - NUM_LIMBS;
                let cs = lo.limbs[index] as u128 + c;
                c = cs >> 64;
                lo.limbs[index] = cs as u64;
            } else {
                let index = 2 * i;
                let cs = hi.limbs[index] as u128 + c;
                c = cs >> 64;
                hi.limbs[index] = cs as u64;
            }
        }
        debug_assert_eq!(c, 0);
        (hi, lo)
    }

    #[inline(always)]
    /// Returns the number of bits needed to represent the number (0 for zero).
    /// If nonzero, this is equivalent to one plus the floored log2 of the number.
    pub const fn bits(&self) -> u32 {
        let mut i = NUM_LIMBS;
        while i > 0 {
            if self.limbs[i - 1] != 0 {
                return i as u32 * u64::BITS - self.limbs[i - 1].leading_zeros();
            }
            i -= 1;
        }
        0
    }

    /// Returns the truthy value if `self != 0` and the falsy value otherwise.
    #[inline]
    const fn ct_is_nonzero(ct: u64) -> u64 {
        Self::ct_from_lsb((ct | ct.wrapping_neg()) >> (u64::BITS - 1))
    }

    /// Returns the truthy value if `value == 1`, and the falsy value if `value == 0`.
    /// Panics for other values.
    const fn ct_from_lsb(value: u64) -> u64 {
        debug_assert!(value == 0 || value == 1);
        value.wrapping_neg()
    }

    /// Return `b` if `c` is truthy, otherwise return `a`.
    #[inline]
    const fn ct_select_limb(a: u64, b: u64, ct: u64) -> u64 {
        a ^ (ct & (a ^ b))
    }

    /// Return `b` if `c` is truthy, otherwise return `a`.
    #[inline]
    const fn ct_select(a: &Self, b: &Self, c: u64) -> Self {
        let mut limbs = [0_u64; NUM_LIMBS];

        let mut i = 0;
        while i < NUM_LIMBS {
            limbs[i] = Self::ct_select_limb(a.limbs[i], b.limbs[i], c);
            i += 1;
        }

        Self { limbs }
    }

    /// Computes `self - (rhs + borrow)`, returning the result along with the new borrow.
    #[inline(always)]
    const fn sbb_limbs(lhs: u64, rhs: u64, borrow: u64) -> (u64, u64) {
        let a = lhs as u128;
        let b = rhs as u128;
        let borrow = (borrow >> (u64::BITS - 1)) as u128;
        let ret = a.wrapping_sub(b + borrow);
        (ret as u64, (ret >> u64::BITS) as u64)
    }

    #[inline(always)]
    /// Computes `a - (b + borrow)`, returning the result along with the new borrow.
    pub fn sbb(&self, rhs: &Self, mut borrow: u64) -> (Self, u64) {
        let mut limbs = [0; NUM_LIMBS];

        for i in (0..NUM_LIMBS).rev() {
            let (w, b) = Self::sbb_limbs(self.limbs[i], rhs.limbs[i], borrow);
            limbs[i] = w;
            borrow = b;
        }

        (Self { limbs }, borrow)
    }

    #[inline(always)]
    /// Returns the number of bits needed to represent the number as little endian
    pub const fn bits_le(&self) -> usize {
        let mut i = 0;
        while i < NUM_LIMBS {
            if self.limbs[i] != 0 {
                return u64::BITS as usize * (NUM_LIMBS - i)
                    - self.limbs[i].leading_zeros() as usize;
            }
            i += 1;
        }
        0
    }

    /// Computes self / rhs, returns the quotient, remainder.
    pub fn div_rem(&self, rhs: &Self) -> (Self, Self) {
        debug_assert!(
            *rhs != UnsignedInteger::from_u64(0),
            "Attempted to divide by zero"
        );
        let mb = rhs.bits_le();
        let mut bd = (NUM_LIMBS * u64::BITS as usize) - mb;
        let mut rem = *self;
        let mut quo = Self::from_u64(0);
        let mut c = rhs.shl(bd);

        loop {
            let (mut r, borrow) = rem.sbb(&c, 0);
            debug_assert!(borrow == 0 || borrow == u64::MAX);
            rem = Self::ct_select(&r, &rem, borrow);
            r = quo.bitor(Self::from_u64(1));
            quo = Self::ct_select(&r, &quo, borrow);
            if bd == 0 {
                break;
            }
            bd -= 1;
            c = c.shr(1);
            quo = quo.shl(1);
        }

        let is_some = Self::ct_is_nonzero(mb as u64);
        quo = Self::ct_select(&Self::from_u64(0), &quo, is_some);
        (quo, rem)
    }

    /// Convert from a decimal string.
    pub fn from_dec_str(value: &str) -> Result<Self, CreationError> {
        if value.is_empty() {
            return Err(CreationError::InvalidDecString);
        }
        let mut res = Self::from_u64(0);
        for b in value.bytes().map(|b| b.wrapping_sub(b'0')) {
            if b > 9 {
                return Err(CreationError::InvalidDecString);
            }
            let (high, low) = Self::mul(&res, &Self::from(10_u64));
            if high > Self::from_u64(0) {
                return Err(CreationError::InvalidDecString);
            }
            res = low + Self::from(b as u64);
        }
        Ok(res)
    }

    #[cfg(feature = "proptest")]
    pub fn nonzero_uint() -> impl Strategy<Value = UnsignedInteger<NUM_LIMBS>> {
        any_uint::<NUM_LIMBS>().prop_filter("is_zero", |&x| x != UnsignedInteger::from_u64(0))
    }
}

impl<const NUM_LIMBS: usize> IsUnsignedInteger for UnsignedInteger<NUM_LIMBS> {}

impl<const NUM_LIMBS: usize> ByteConversion for UnsignedInteger<NUM_LIMBS> {
    #[cfg(feature = "alloc")]
    fn to_bytes_be(&self) -> alloc::vec::Vec<u8> {
        self.limbs
            .iter()
            .flat_map(|limb| limb.to_be_bytes())
            .collect()
    }

    #[cfg(feature = "alloc")]
    fn to_bytes_le(&self) -> alloc::vec::Vec<u8> {
        self.limbs
            .iter()
            .rev()
            .flat_map(|limb| limb.to_le_bytes())
            .collect()
    }

    fn from_bytes_be(bytes: &[u8]) -> Result<Self, ByteConversionError> {
        // We cut off extra bytes, this is useful when you use this function to generate the element from randomness
        // In the future with the right algorithm this shouldn't be needed

        let needed_bytes = bytes
            .get(0..NUM_LIMBS * 8)
            .ok_or(ByteConversionError::FromBEBytesError)?;

        let mut limbs: [u64; NUM_LIMBS] = [0; NUM_LIMBS];

        needed_bytes
            .chunks_exact(8)
            .enumerate()
            .try_for_each(|(i, chunk)| {
                let limb = u64::from_be_bytes(
                    chunk
                        .try_into()
                        .map_err(|_| ByteConversionError::FromBEBytesError)?,
                );
                limbs[i] = limb;
                Ok::<_, ByteConversionError>(())
            })?;

        Ok(Self { limbs })
    }

    fn from_bytes_le(bytes: &[u8]) -> Result<Self, ByteConversionError> {
        let needed_bytes = bytes
            .get(0..NUM_LIMBS * 8)
            .ok_or(ByteConversionError::FromBEBytesError)?;

        let mut limbs: [u64; NUM_LIMBS] = [0; NUM_LIMBS];

        needed_bytes
            .chunks_exact(8)
            .rev()
            .enumerate()
            .try_for_each(|(i, chunk)| {
                let limb = u64::from_le_bytes(
                    chunk
                        .try_into()
                        .map_err(|_| ByteConversionError::FromLEBytesError)?,
                );
                limbs[i] = limb;
                Ok::<_, ByteConversionError>(())
            })?;

        Ok(Self { limbs })
    }
}

impl<const NUM_LIMBS: usize> From<UnsignedInteger<NUM_LIMBS>> for u16 {
    fn from(value: UnsignedInteger<NUM_LIMBS>) -> Self {
        value.limbs[NUM_LIMBS - 1] as u16
    }
}

#[cfg(feature = "alloc")]
impl<const NUM_LIMBS: usize> AsBytes for UnsignedInteger<NUM_LIMBS> {
    fn as_bytes(&self) -> alloc::vec::Vec<u8> {
        self.limbs
            .into_iter()
            .fold(alloc::vec::Vec::new(), |mut acc, limb| {
                acc.extend_from_slice(&limb.as_bytes());
                acc
            })
    }
}

#[cfg(feature = "alloc")]
impl<const NUM_LIMBS: usize> From<UnsignedInteger<NUM_LIMBS>> for alloc::vec::Vec<u8> {
    fn from(val: UnsignedInteger<NUM_LIMBS>) -> Self {
        val.as_bytes()
    }
}

#[cfg(feature = "proptest")]
fn any_uint<const NUM_LIMBS: usize>() -> impl Strategy<Value = UnsignedInteger<NUM_LIMBS>> {
    any::<[u64; NUM_LIMBS]>().prop_map(UnsignedInteger::from_limbs)
}

#[cfg(feature = "proptest")]
impl<const NUM_LIMBS: usize> Arbitrary for UnsignedInteger<NUM_LIMBS> {
    type Parameters = ();

    fn arbitrary_with(_args: Self::Parameters) -> Self::Strategy {
        any_uint::<NUM_LIMBS>().sboxed()
    }

    type Strategy = SBoxedStrategy<Self>;
}

#[cfg(test)]
mod tests_u384 {
    use crate::traits::ByteConversion;
    use crate::unsigned_integer::element::{UnsignedInteger, U256, U384};
    #[cfg(feature = "proptest")]
    use proptest::prelude::*;
    #[cfg(feature = "proptest")]
    use std::ops::Shr;

    #[cfg(feature = "proptest")]
    const N_LIMBS: usize = 8;
    #[cfg(feature = "proptest")]
    type Uint = UnsignedInteger<N_LIMBS>;

    #[cfg(feature = "proptest")]
    proptest! {
        #[test]
        fn bitand(a in any::<Uint>(), b in any::<Uint>()) {
            let result = a & b;

            for i in 0..N_LIMBS {
                assert_eq!(result.limbs[i], a.limbs[i] & b.limbs[i]);
            }
        }

        #[test]
        fn bitand_assign(a in any::<Uint>(), b in any::<Uint>()) {
            let mut result = a;
            result &= b;

            for i in 0..N_LIMBS {
                assert_eq!(result.limbs[i], a.limbs[i] & b.limbs[i]);
            }
        }

        #[test]
        fn bitor(a in any::<Uint>(), b in any::<Uint>()) {
            let result = a | b;

            for i in 0..N_LIMBS {
                assert_eq!(result.limbs[i], a.limbs[i] | b.limbs[i]);
            }
        }

        #[test]
        fn bitor_assign(a in any::<Uint>(), b in any::<Uint>()) {
            let mut result = a;
            result |= b;

            for i in 0..N_LIMBS {
                assert_eq!(result.limbs[i], a.limbs[i] | b.limbs[i]);
            }
        }

        #[test]
        fn bitxor(a in any::<Uint>(), b in any::<Uint>()) {
            let result = a ^ b;

            for i in 0..N_LIMBS {
                assert_eq!(result.limbs[i], a.limbs[i] ^ b.limbs[i]);
            }
        }

        #[test]
        fn bitxor_assign(a in any::<Uint>(), b in any::<Uint>()) {
            let mut result = a;
            result ^= b;

            for i in 0..N_LIMBS {
                assert_eq!(result.limbs[i], a.limbs[i] ^ b.limbs[i]);
            }
        }

        #[test]
        fn div_rem(a in any::<Uint>(), b in any::<Uint>()) {
            let a = a.shr(256);
            let b = b.shr(256);
            assert_eq!((a * b).div_rem(&b), (a, Uint::from_u64(0)));
        }
    }

    #[test]
    fn construct_new_integer_from_limbs() {
        let a: U384 = UnsignedInteger {
            limbs: [0, 1, 2, 3, 4, 5],
        };
        assert_eq!(U384::from_limbs([0, 1, 2, 3, 4, 5]), a);
    }

    #[test]
    fn construct_new_integer_from_u64_1() {
        let a = U384::from_u64(1_u64);
        assert_eq!(a.limbs, [0, 0, 0, 0, 0, 1]);
    }

    #[test]
    fn construct_new_integer_from_u54_2() {
        let a = U384::from_u64(u64::MAX);
        assert_eq!(a.limbs, [0, 0, 0, 0, 0, u64::MAX]);
    }

    #[test]
    fn construct_new_integer_from_u128_1() {
        let a = U384::from_u128(u128::MAX);
        assert_eq!(a.limbs, [0, 0, 0, 0, u64::MAX, u64::MAX]);
    }

    #[test]
    fn construct_new_integer_from_u128_2() {
        let a = U384::from_u128(276371540478856090688472252609570374439);
        assert_eq!(
            a.limbs,
            [0, 0, 0, 0, 14982131230017065096, 14596400355126379303]
        );
    }

    #[test]
    fn construct_new_integer_from_hex_1() {
        let a = U384::from_hex_unchecked("1");
        assert_eq!(a.limbs, [0, 0, 0, 0, 0, 1]);
    }

    #[test]
    fn construct_new_integer_from_zero_x_1() {
        let a = U384::from_hex_unchecked("0x1");
        assert_eq!(a.limbs, [0, 0, 0, 0, 0, 1]);
    }

    #[test]
    fn construct_new_integer_from_hex_2() {
        let a = U384::from_hex_unchecked("f");
        assert_eq!(a.limbs, [0, 0, 0, 0, 0, 15]);
    }

    #[test]
    fn construct_new_integer_from_hex_3() {
        let a = U384::from_hex_unchecked("10000000000000000");
        assert_eq!(a.limbs, [0, 0, 0, 0, 1, 0]);
    }

    #[test]
    fn construct_new_integer_from_hex_4() {
        let a = U384::from_hex_unchecked("a0000000000000000");
        assert_eq!(a.limbs, [0, 0, 0, 0, 10, 0]);
    }

    #[test]
    fn construct_new_integer_from_hex_5() {
        let a = U384::from_hex_unchecked("ffffffffffffffffff");
        assert_eq!(a.limbs, [0, 0, 0, 0, 255, u64::MAX]);
    }

    #[test]
    fn construct_new_integer_from_hex_6() {
        let a = U384::from_hex_unchecked("eb235f6144d9e91f4b14");
        assert_eq!(a.limbs, [0, 0, 0, 0, 60195, 6872850209053821716]);
    }

    #[test]
    fn construct_new_integer_from_hex_7() {
        let a = U384::from_hex_unchecked("2b20aaa5cf482b239e2897a787faf4660cc95597854beb2");
        assert_eq!(
            a.limbs,
            [
                0,
                0,
                0,
                194229460750598834,
                4171047363999149894,
                6975114134393503410
            ]
        );
    }

    #[test]
    fn construct_new_integer_from_hex_checked_7() {
        let a = U384::from_hex("2b20aaa5cf482b239e2897a787faf4660cc95597854beb2").unwrap();
        assert_eq!(
            a.limbs,
            [
                0,
                0,
                0,
                194229460750598834,
                4171047363999149894,
                6975114134393503410
            ]
        );
    }

    #[test]
    fn construct_new_integer_from_hex_checked_7_with_zero_x() {
        let a = U384::from_hex("0x2b20aaa5cf482b239e2897a787faf4660cc95597854beb2").unwrap();
        assert_eq!(
            a.limbs,
            [
                0,
                0,
                0,
                194229460750598834,
                4171047363999149894,
                6975114134393503410
            ]
        );
    }

    #[test]
    fn construct_new_integer_from_non_hex_errs() {
        assert!(U384::from_hex("0xTEST").is_err());
    }

    #[test]
    fn construct_new_integer_from_empty_string_errs() {
        assert!(U384::from_hex("").is_err());
    }

    #[test]
    fn construct_new_integer_from_hex_checked_8() {
        let a = U384::from_hex("140f5177b90b4f96b61bb8ccb4f298ad2b20aaa5cf482b239e2897a787faf4660cc95597854beb235f6144d9e91f4b14").unwrap();
        assert_eq!(
            a.limbs,
            [
                1445463580056702870,
                13122285128622708909,
                3107671372009581347,
                11396525602857743462,
                921361708038744867,
                6872850209053821716
            ]
        );
    }

    #[test]
    fn construct_new_integer_from_hex_8() {
        let a = U384::from_hex_unchecked("140f5177b90b4f96b61bb8ccb4f298ad2b20aaa5cf482b239e2897a787faf4660cc95597854beb235f6144d9e91f4b14");
        assert_eq!(
            a.limbs,
            [
                1445463580056702870,
                13122285128622708909,
                3107671372009581347,
                11396525602857743462,
                921361708038744867,
                6872850209053821716
            ]
        );
    }

    #[test]
    fn from_hex_with_overflowing_hexstring_should_error() {
        let u256_from_big_string = U256::from_hex(&"f".repeat(65));
        assert!(u256_from_big_string.is_err());
        assert!(
            u256_from_big_string
                == Err(crate::unsigned_integer::element::CreationError::HexStringIsTooBig)
        );
    }

    #[test]
    fn from_hex_with_non_overflowing_hexstring_should_work() {
        assert_eq!(U256::from_hex(&"0".repeat(64)).unwrap().limbs, [0, 0, 0, 0])
    }

    #[test]
    fn construct_new_integer_from_dec_1() {
        let a = U384::from_dec_str("1").unwrap();
        assert_eq!(a.limbs, [0, 0, 0, 0, 0, 1]);
    }

    #[test]
    fn construct_new_integer_from_dec_2() {
        let a = U384::from_dec_str("15").unwrap();
        assert_eq!(a.limbs, [0, 0, 0, 0, 0, 15]);
    }

    #[test]
    fn construct_new_integer_from_dec_3() {
        let a = U384::from_dec_str("18446744073709551616").unwrap();
        assert_eq!(a.limbs, [0, 0, 0, 0, 1, 0]);
    }

    #[test]
    fn construct_new_integer_from_dec_4() {
        let a = U384::from_dec_str("184467440737095516160").unwrap();
        assert_eq!(a.limbs, [0, 0, 0, 0, 10, 0]);
    }

    #[test]
    fn construct_new_integer_from_dec_5() {
        let a = U384::from_dec_str("4722366482869645213695").unwrap();
        assert_eq!(a.limbs, [0, 0, 0, 0, 255, u64::MAX]);
    }

    #[test]
    fn construct_new_integer_from_dec_6() {
        let a = U384::from_dec_str("1110408632367155513346836").unwrap();
        assert_eq!(a.limbs, [0, 0, 0, 0, 60195, 6872850209053821716]);
    }

    #[test]
    fn construct_new_integer_from_dec_7() {
        let a =
            U384::from_dec_str("66092860629991288370279803883558073888453977263446474418").unwrap();
        assert_eq!(
            a.limbs,
            [
                0,
                0,
                0,
                194229460750598834,
                4171047363999149894,
                6975114134393503410
            ]
        );
    }

    #[test]
    fn construct_new_integer_from_dec_8() {
        let a = U384::from_dec_str("3087491467896943881295768554872271030441880044814691421073017731442549147034464936390742057449079000462340371991316").unwrap();
        assert_eq!(
            a.limbs,
            [
                1445463580056702870,
                13122285128622708909,
                3107671372009581347,
                11396525602857743462,
                921361708038744867,
                6872850209053821716
            ]
        );
    }

    #[test]
    fn construct_new_integer_from_dec_empty() {
        assert!(U384::from_dec_str("").is_err());
    }

    #[test]
    fn construct_new_integer_from_dec_invalid() {
        assert!(U384::from_dec_str("0xff").is_err());
    }

    #[test]
    fn equality_works_1() {
        let a = U384::from_hex_unchecked("1");
        let b = U384 {
            limbs: [0, 0, 0, 0, 0, 1],
        };
        assert_eq!(a, b);
    }
    #[test]
    fn equality_works_2() {
        let a = U384::from_hex_unchecked("f");
        let b = U384 {
            limbs: [0, 0, 0, 0, 0, 15],
        };
        assert_eq!(a, b);
    }

    #[test]
    fn equality_works_3() {
        let a = U384::from_hex_unchecked("10000000000000000");
        let b = U384 {
            limbs: [0, 0, 0, 0, 1, 0],
        };
        assert_eq!(a, b);
    }

    #[test]
    fn equality_works_4() {
        let a = U384::from_hex_unchecked("a0000000000000000");
        let b = U384 {
            limbs: [0, 0, 0, 0, 10, 0],
        };
        assert_eq!(a, b);
    }

    #[test]
    fn equality_works_5() {
        let a = U384::from_hex_unchecked("ffffffffffffffffff");
        let b = U384 {
            limbs: [0, 0, 0, 0, u8::MAX as u64, u64::MAX],
        };
        assert_eq!(a, b);
    }

    #[test]
    fn equality_works_6() {
        let a = U384::from_hex_unchecked("eb235f6144d9e91f4b14");
        let b = U384 {
            limbs: [0, 0, 0, 0, 60195, 6872850209053821716],
        };
        assert_eq!(a, b);
    }

    #[test]
    fn equality_works_7() {
        let a = U384::from_hex_unchecked("2b20aaa5cf482b239e2897a787faf4660cc95597854beb2");
        let b = U384 {
            limbs: [
                0,
                0,
                0,
                194229460750598834,
                4171047363999149894,
                6975114134393503410,
            ],
        };
        assert_eq!(a, b);
    }

    #[test]
    fn equality_works_8() {
        let a = U384::from_hex_unchecked("140f5177b90b4f96b61bb8ccb4f298ad2b20aaa5cf482b239e2897a787faf4660cc95597854beb235f6144d9e91f4b14");
        let b = U384 {
            limbs: [
                1445463580056702870,
                13122285128622708909,
                3107671372009581347,
                11396525602857743462,
                921361708038744867,
                6872850209053821716,
            ],
        };
        assert_eq!(a, b);
    }

    #[test]
    fn equality_works_9() {
        let a = U384::from_hex_unchecked("fffffff");
        let b = U384::from_hex_unchecked("fefffff");
        assert_ne!(a, b);
    }

    #[test]
    fn equality_works_10() {
        let a = U384::from_hex_unchecked("ffff000000000000");
        let b = U384::from_hex_unchecked("ffff000000100000");
        assert_ne!(a, b);
    }

    #[test]
    fn const_ne_works_1() {
        let a = U384::from_hex_unchecked("ffff000000000000");
        let b = U384::from_hex_unchecked("ffff000000100000");
        assert!(U384::const_ne(&a, &b));
    }

    #[test]
    fn const_ne_works_2() {
        let a = U384::from_hex_unchecked("140f5177b90b4f96b61bb8ccb4f298ad2b20aaa5cf482b239e2897a787faf4660cc95597854beb235f6144d9e91f4b14");
        let b = U384 {
            limbs: [
                1445463580056702870,
                13122285128622708909,
                3107671372009581347,
                11396525602857743462,
                921361708038744867,
                6872850209053821716,
            ],
        };
        assert!(!U384::const_ne(&a, &b));
    }

    #[test]
    fn double_two_384_bit_integers() {
        let a = U384::from_u64(2);
        let b = U384::from_u64(5);
        let c = U384::from_u64(7);
        assert_eq!(U384::double(&a).0, a + a);
        assert_eq!(U384::double(&b).0, b + b);
        assert_eq!(U384::double(&c).0, c + c);
    }

    #[test]
    fn add_two_384_bit_integers_1() {
        let a = U384::from_u64(2);
        let b = U384::from_u64(5);
        let c = U384::from_u64(7);
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_384_bit_integers_2() {
        let a = U384::from_u64(334);
        let b = U384::from_u64(666);
        let c = U384::from_u64(1000);
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_384_bit_integers_3() {
        let a = U384::from_hex_unchecked("ffffffffffffffff");
        let b = U384::from_hex_unchecked("1");
        let c = U384::from_hex_unchecked("10000000000000000");
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_384_bit_integers_4() {
        let a = U384::from_hex_unchecked("b58e1e0b66");
        let b = U384::from_hex_unchecked("55469d9619");
        let c = U384::from_hex_unchecked("10ad4bba17f");
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_384_bit_integers_5() {
        let a = U384::from_hex_unchecked("e8dff25cb6160f7705221da6f");
        let b = U384::from_hex_unchecked("ab879169b5f80dc8a7969f0b0");
        let c = U384::from_hex_unchecked("1946783c66c0e1d3facb8bcb1f");
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_384_bit_integers_6() {
        let a = U384::from_hex_unchecked("9adf291af3a64d59e14e7b440c850508014c551ed5");
        let b = U384::from_hex_unchecked("e7948474bce907f0feaf7e5d741a8cd2f6d1fb9448");
        let c = U384::from_hex_unchecked("18273ad8fb08f554adffdf9a1809f91daf81e50b31d");
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_384_bit_integers_7() {
        let a = U384::from_hex_unchecked(
            "f866aef803c92bf02e85c7fad0eccb4881c59825e499fa22f98e1a8fefed4cd9a03647cd3cc84",
        );
        let b = U384::from_hex_unchecked(
            "9b4000dccf01a010e196154a1b998408f949d734389626ba97cb3331ee87e01dd5badc58f41b2",
        );
        let c = U384::from_hex_unchecked(
            "193a6afd4d2cacc01101bdd44ec864f517b0f6f5a1d3020dd91594dc1de752cf775f1242630e36",
        );
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_384_bit_integers_8() {
        let a = U384::from_hex_unchecked("07df9c74fa9d5aafa74a87dbbf93215659d8a3e1706d4b06de9512284802580eb36ae12ea59f90db5b1799d0970a42e");
        let b = U384::from_hex_unchecked("d515e54973f0643a6a9957579c1f84020a6a91d5d5f27b75401c7538d2c9ea9cafff44a2c606877d46c49a3433cc85e");
        let c = U384::from_hex_unchecked("dcf581be6e8dbeea11e3df335bb2a558644335b7465fc67c1eb187611acc42ab636a25d16ba61858a1dc3404cad6c8c");
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_384_bit_integers_9() {
        let a = U384::from_hex_unchecked("92977527a0f8ba00d18c1b2f1900d965d4a70e5f5f54468ffb2d4d41519385f24b078a0e7d0281d5ad0c36724dc4233");
        let b = U384::from_hex_unchecked("46facf9953a9494822bf18836ffd7e55c48b30aa81e17fa1ace0b473015307e4622b8bd6fa68ef654796a183abde842");
        let c = U384::from_hex_unchecked("d99244c0f4a20348f44b33b288fe57bb99323f09e135c631a80e01b452e68dd6ad3315e5776b713af4a2d7f5f9a2a75");
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_384_bit_integers_10() {
        let a = U384::from_hex_unchecked("07df9c74fa9d5aafa74a87dbbf93215659d8a3e1706d4b06de9512284802580eb36ae12ea59f90db5b1799d0970a42e");
        let b = U384::from_hex_unchecked("d515e54973f0643a6a9957579c1f84020a6a91d5d5f27b75401c7538d2c9ea9cafff44a2c606877d46c49a3433cc85e");
        let c_expected = U384::from_hex_unchecked("dcf581be6e8dbeea11e3df335bb2a558644335b7465fc67c1eb187611acc42ab636a25d16ba61858a1dc3404cad6c8c");
        let (c, overflow) = U384::add(&a, &b);
        assert_eq!(c, c_expected);
        assert!(!overflow);
    }

    #[test]
    fn add_two_384_bit_integers_11() {
        let a = U384::from_hex_unchecked("92977527a0f8ba00d18c1b2f1900d965d4a70e5f5f54468ffb2d4d41519385f24b078a0e7d0281d5ad0c36724dc4233");
        let b = U384::from_hex_unchecked("46facf9953a9494822bf18836ffd7e55c48b30aa81e17fa1ace0b473015307e4622b8bd6fa68ef654796a183abde842");
        let c_expected = U384::from_hex_unchecked("d99244c0f4a20348f44b33b288fe57bb99323f09e135c631a80e01b452e68dd6ad3315e5776b713af4a2d7f5f9a2a75");
        let (c, overflow) = U384::add(&a, &b);
        assert_eq!(c, c_expected);
        assert!(!overflow);
    }

    #[test]
    fn add_two_384_bit_integers_12_with_overflow() {
        let a = U384::from_hex_unchecked("b07bc844363dd56467d9ebdd5929e9bb34a8e2577db77df6cf8f2ac45bd3d0bc2fc3078d265fe761af51d6aec5b59428");
        let b = U384::from_hex_unchecked("cbbc474761bb7995ff54e25fa5d30295604fe3545d0cde405e72d8c0acebb119e9158131679b6c34483a3dafb49deeea");
        let c_expected = U384::from_hex_unchecked("7c380f8b97f94efa672ece3cfefcec5094f8c5abdac45c372e02038508bf81d618d888be8dfb5395f78c145e7a538312");
        let (c, overflow) = U384::add(&a, &b);
        assert_eq!(c, c_expected);
        assert!(overflow);
    }

    #[test]
    fn double_384_bit_integer_12_with_overflow() {
        let a = U384::from_hex_unchecked("b07bc844363dd56467d9ebdd5929e9bb34a8e2577db77df6cf8f2ac45bd3d0bc2fc3078d265fe761af51d6aec5b59428");
        assert_eq!(U384::double(&a), U384::add(&a, &a));
    }

    #[test]
    fn sub_two_384_bit_integers_1() {
        let a = U384::from_u64(2);
        let b = U384::from_u64(5);
        let c = U384::from_u64(7);
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_384_bit_integers_2() {
        let a = U384::from_u64(334);
        let b = U384::from_u64(666);
        let c = U384::from_u64(1000);
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_384_bit_integers_3() {
        let a = U384::from_hex_unchecked("ffffffffffffffff");
        let b = U384::from_hex_unchecked("1");
        let c = U384::from_hex_unchecked("10000000000000000");
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_384_bit_integers_4() {
        let a = U384::from_hex_unchecked("b58e1e0b66");
        let b = U384::from_hex_unchecked("55469d9619");
        let c = U384::from_hex_unchecked("10ad4bba17f");
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_384_bit_integers_5() {
        let a = U384::from_hex_unchecked("e8dff25cb6160f7705221da6f");
        let b = U384::from_hex_unchecked("ab879169b5f80dc8a7969f0b0");
        let c = U384::from_hex_unchecked("1946783c66c0e1d3facb8bcb1f");
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_384_bit_integers_6() {
        let a = U384::from_hex_unchecked("9adf291af3a64d59e14e7b440c850508014c551ed5");
        let b = U384::from_hex_unchecked("e7948474bce907f0feaf7e5d741a8cd2f6d1fb9448");
        let c = U384::from_hex_unchecked("18273ad8fb08f554adffdf9a1809f91daf81e50b31d");
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_384_bit_integers_7() {
        let a = U384::from_hex_unchecked(
            "f866aef803c92bf02e85c7fad0eccb4881c59825e499fa22f98e1a8fefed4cd9a03647cd3cc84",
        );
        let b = U384::from_hex_unchecked(
            "9b4000dccf01a010e196154a1b998408f949d734389626ba97cb3331ee87e01dd5badc58f41b2",
        );
        let c = U384::from_hex_unchecked(
            "193a6afd4d2cacc01101bdd44ec864f517b0f6f5a1d3020dd91594dc1de752cf775f1242630e36",
        );
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_384_bit_integers_8() {
        let a = U384::from_hex_unchecked("07df9c74fa9d5aafa74a87dbbf93215659d8a3e1706d4b06de9512284802580eb36ae12ea59f90db5b1799d0970a42e");
        let b = U384::from_hex_unchecked("d515e54973f0643a6a9957579c1f84020a6a91d5d5f27b75401c7538d2c9ea9cafff44a2c606877d46c49a3433cc85e");
        let c = U384::from_hex_unchecked("dcf581be6e8dbeea11e3df335bb2a558644335b7465fc67c1eb187611acc42ab636a25d16ba61858a1dc3404cad6c8c");
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_384_bit_integers_9() {
        let a = U384::from_hex_unchecked("92977527a0f8ba00d18c1b2f1900d965d4a70e5f5f54468ffb2d4d41519385f24b078a0e7d0281d5ad0c36724dc4233");
        let b = U384::from_hex_unchecked("46facf9953a9494822bf18836ffd7e55c48b30aa81e17fa1ace0b473015307e4622b8bd6fa68ef654796a183abde842");
        let c = U384::from_hex_unchecked("d99244c0f4a20348f44b33b288fe57bb99323f09e135c631a80e01b452e68dd6ad3315e5776b713af4a2d7f5f9a2a75");
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_384_bit_integers_11_without_overflow() {
        let a = U384::from_u64(334);
        let b_expected = U384::from_u64(666);
        let c = U384::from_u64(1000);
        let (b, underflow) = U384::sub(&c, &a);
        assert!(!underflow);
        assert_eq!(b_expected, b);
    }

    #[test]
    fn sub_two_384_bit_integers_11_with_overflow() {
        let a = U384::from_u64(334);
        let b_expected = U384::from_hex_unchecked("fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffd66");
        let c = U384::from_u64(1000);
        let (b, underflow) = U384::sub(&a, &c);
        assert!(underflow);
        assert_eq!(b_expected, b);
    }

    #[test]
    fn partial_order_works() {
        assert!(U384::from_u64(10) <= U384::from_u64(10));
        assert!(U384::from_u64(1) < U384::from_u64(2));
        assert!(U384::from_u64(2) >= U384::from_u64(1));

        assert!(U384::from_u64(10) >= U384::from_u64(10));
        assert!(U384::from_u64(2) > U384::from_u64(1));
        assert!(U384::from_u64(1) <= U384::from_u64(2));

        let a = U384::from_hex_unchecked("92977527a0f8ba00d18c1b2f1900d965d4a70e5f5f54468ffb2d4d41519385f24b078a0e7d0281d5ad0c36724dc4233");
        let c = U384::from_hex_unchecked("d99244c0f4a20348f44b33b288fe57bb99323f09e135c631a80e01b452e68dd6ad3315e5776b713af4a2d7f5f9a2a75");

        assert!(&a <= &a);
        assert!(&a >= &a);
        assert!(&a >= &a);
        assert!(&a <= &a);
        assert!(&a < &(&a + U384::from_u64(1)));
        assert!(&a <= &(&a + U384::from_u64(1)));
        assert!(&a + U384::from_u64(1) > a);
        assert!((&a + U384::from_u64(1) >= a));
        assert!(&a <= &c);
        assert!(&a < &c);
        assert!(&a < &c);
        assert!(&a <= &c);
        assert!(&c > &a);
        assert!(&c >= &a);
        assert!(&c >= &a);
        assert!(&c > &a);
        assert!(a < c);
    }

    #[test]
    fn mul_two_384_bit_integers_works_1() {
        let a = U384::from_u64(3);
        let b = U384::from_u64(8);
        let c = U384::from_u64(3 * 8);
        assert_eq!(a * b, c);
    }

    #[test]
    fn mul_two_384_bit_integers_works_2() {
        let a = U384::from_hex_unchecked("6131d99f840b3b0");
        let b = U384::from_hex_unchecked("6f5c466db398f43");
        let c = U384::from_hex_unchecked("2a47a603a77f871dfbb937af7e5710");
        assert_eq!(a * b, c);
    }

    #[test]
    fn mul_two_384_bit_integers_works_3() {
        let a = U384::from_hex_unchecked("84a6add5db9e095b2e0f6b40eff8ee");
        let b = U384::from_hex_unchecked("2347db918f725461bec2d5c57");
        let c = U384::from_hex_unchecked("124805c476c9462adc0df6c88495d4253f5c38033afc18d78d920e2");
        assert_eq!(a * b, c);
    }

    #[test]
    fn mul_two_384_bit_integers_works_4() {
        let a = U384::from_hex_unchecked("04050753dd7c0b06c404633016f87040");
        let b = U384::from_hex_unchecked("dc3830be041b3b4476445fcad3dac0f6f3a53e4ba12da");
        let c = U384::from_hex_unchecked(
            "375342999dab7f52f4010c4abc2e18b55218015931a55d6053ac39e86e2a47d6b1cb95f41680",
        );
        assert_eq!(a * b, c);
    }

    #[test]
    fn mul_two_384_bit_integers_works_5() {
        let a = U384::from_hex_unchecked("7ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff8");
        let b = U384::from_hex_unchecked("2");
        let c_expected = U384::from_hex_unchecked(
            "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff0",
        );
        assert_eq!(a * b, c_expected);
    }

    #[test]
    #[should_panic]
    fn mul_two_384_bit_integers_works_6() {
        let a = U384::from_hex_unchecked("800000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
        let b = U384::from_hex_unchecked("2");
        let _c = a * b;
    }

    #[test]
    fn mul_two_384_bit_integers_works_7_hi_lo() {
        let a = U384::from_hex_unchecked("04050753dd7c0b06c404633016f87040");
        let b = U384::from_hex_unchecked("dc3830be041b3b4476445fcad3dac0f6f3a53e4ba12da");
        let hi_expected = U384::from_hex_unchecked("0");
        let lo_expected = U384::from_hex_unchecked(
            "375342999dab7f52f4010c4abc2e18b55218015931a55d6053ac39e86e2a47d6b1cb95f41680",
        );
        let (hi, lo) = U384::mul(&a, &b);
        assert_eq!(hi, hi_expected);
        assert_eq!(lo, lo_expected);
    }

    #[test]
    fn mul_two_384_bit_integers_works_8_hi_lo() {
        let a = U384::from_hex_unchecked("5e2d939b602a50911232731d04fe6f40c05f97da0602307099fb991f9b414e2d52bef130349ec18db1a0215ea6caf76");
        let b = U384::from_hex_unchecked("3f3ad1611ab58212f92a2484e9560935b9ac4615fe61cfed1a4861e193a74d20c94f9f88d8b2cc089543c3f699969d9");
        let hi_expected = U384::from_hex_unchecked(
            "1742daad9c7861dd3499e7ece65467e337937b27e20d641b225bfe00323d33ed62715654eadc092b057a5f19f2ad6c",
        );
        let lo_expected = U384::from_hex_unchecked("9969c0417b9304d9c16b046c860447d3533999e16710d2e90a44959a168816c015ffb44b987e8cbb82bd46b08d9e2106");
        let (hi, lo) = U384::mul(&a, &b);
        assert_eq!(hi, hi_expected);
        assert_eq!(lo, lo_expected);
    }

    #[test]
    fn mul_two_384_bit_integers_works_9_hi_lo() {
        let a = U384::from_hex_unchecked("800000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
        let b = U384::from_hex_unchecked("2");
        let hi_expected = U384::from_hex_unchecked("1");
        let lo_expected = U384::from_hex_unchecked("0");
        let (hi, lo) = U384::mul(&a, &b);
        assert_eq!(hi, hi_expected);
        assert_eq!(lo, lo_expected);
    }

    #[test]
    fn shift_left_on_384_bit_integer_works_1() {
        let a = U384::from_hex_unchecked("1");
        let b = U384::from_hex_unchecked("10");
        assert_eq!(a << 4, b);
    }

    #[test]
    fn shift_left_on_384_bit_integer_works_2() {
        let a = U384::from_u64(1);
        let b = U384::from_u128(1_u128 << 64);
        assert_eq!(a << 64, b);
    }

    #[test]
    fn shift_left_on_384_bit_integer_works_3() {
        let a = U384::from_hex_unchecked("10");
        let b = U384::from_hex_unchecked("1000");
        assert_eq!(&a << 8, b);
    }

    #[test]
    fn shift_left_on_384_bit_integer_works_4() {
        let a = U384::from_hex_unchecked("e45542992b6844553f3cb1c5ac33e7fa5");
        let b = U384::from_hex_unchecked("391550a64ada11154fcf2c716b0cf9fe940");
        assert_eq!(a << 6, b);
    }

    #[test]
    fn shift_left_on_384_bit_integer_works_5() {
        let a = U384::from_hex_unchecked(
            "03303f4d6c2d1caf0c24a6b0239b679a8390aa99bead76bc0093b1bc1a8101f5ce",
        );
        let b = U384::from_hex_unchecked("6607e9ad85a395e18494d604736cf35072155337d5aed7801276378350203eb9c0000000000000000000000000000000");
        assert_eq!(&a << 125, b);
    }

    #[test]
    fn shift_left_on_384_bit_integer_works_6() {
        let a = U384::from_hex_unchecked("762e8968bc392ed786ab132f0b5b0cacd385dd51de3a");
        let b = U384::from_hex_unchecked(
            "762e8968bc392ed786ab132f0b5b0cacd385dd51de3a00000000000000000000000000000000",
        );
        assert_eq!(&a << (64 * 2), b);
    }

    #[test]
    fn shift_left_on_384_bit_integer_works_7() {
        let a = U384::from_hex_unchecked("90823e0bd707f");
        let b = U384::from_hex_unchecked(
            "90823e0bd707f000000000000000000000000000000000000000000000000",
        );
        assert_eq!(&a << (64 * 3), b);
    }

    #[test]
    fn shift_right_on_384_bit_integer_works_1() {
        let a = U384::from_hex_unchecked("1");
        let b = U384::from_hex_unchecked("10");
        assert_eq!(b >> 4, a);
    }

    #[test]
    fn shift_right_on_384_bit_integer_works_2() {
        let a = U384::from_hex_unchecked("10");
        let b = U384::from_hex_unchecked("1000");
        assert_eq!(&b >> 8, a);
    }

    #[test]
    fn shift_right_on_384_bit_integer_works_3() {
        let a = U384::from_hex_unchecked("e45542992b6844553f3cb1c5ac33e7fa5");
        let b = U384::from_hex_unchecked("391550a64ada11154fcf2c716b0cf9fe940");
        assert_eq!(b >> 6, a);
    }

    #[test]
    fn shift_right_on_384_bit_integer_works_4() {
        let a = U384::from_hex_unchecked(
            "03303f4d6c2d1caf0c24a6b0239b679a8390aa99bead76bc0093b1bc1a8101f5ce",
        );
        let b = U384::from_hex_unchecked("6607e9ad85a395e18494d604736cf35072155337d5aed7801276378350203eb9c0000000000000000000000000000000");
        assert_eq!(&b >> 125, a);
    }

    #[test]
    fn shift_right_on_384_bit_integer_works_5() {
        let a = U384::from_hex_unchecked("ba6ab46f9a9a2f20e4061b67ce4d8c3da98091cf990d7b14ef47ffe27370abbdeb6a3ce9f9cbf5df1b2430114c8558eb");
        let b =
            U384::from_hex_unchecked("174d568df35345e41c80c36cf9c9b187b5301239f321af629de8fffc4e6");
        assert_eq!(a >> 151, b);
    }

    #[test]
    fn shift_right_on_384_bit_integer_works_6() {
        let a = U384::from_hex_unchecked(
            "076c075d2f65e39b9ecdde8bf6f8c94241962ce0f557b7739673200c777152eb7e772ad35",
        );
        let b = U384::from_hex_unchecked("ed80eba5ecbc7373d9bbd17edf19284832c59c1eaaf6ee7");
        assert_eq!(&a >> 99, b);
    }

    #[test]
    fn shift_right_on_384_bit_integer_works_7() {
        let a = U384::from_hex_unchecked("6a9ce35d8940a5ebd29604ce9a182ade76f03f7e9965760b84a8cfd1d3dd2e612669fe000e58b2af688fd90");
        let b = U384::from_hex_unchecked("6a9ce35d8940a5ebd29604ce9a182ade76f03f7");
        assert_eq!(&a >> (64 * 3), b);
    }

    #[test]
    fn shift_right_on_384_bit_integer_works_8() {
        let a = U384::from_hex_unchecked("5322c128ec84081b6c376c108ebd7fd36bbd44f71ee5e6ad6bcb3dd1c5265bd7db75c90b2665a0826d17600f0e9");
        let b =
            U384::from_hex_unchecked("5322c128ec84081b6c376c108ebd7fd36bbd44f71ee5e6ad6bcb3dd1c52");
        assert_eq!(&a >> (64 * 2), b);
    }

    #[test]
    #[cfg(feature = "alloc")]
    fn to_be_bytes_works() {
        let number = U384::from_u64(1);
        let expected_bytes = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        ];

        assert_eq!(number.to_bytes_be(), expected_bytes);
    }

    #[test]
    #[cfg(feature = "alloc")]
    fn to_le_bytes_works() {
        let number = U384::from_u64(1);
        let expected_bytes = [
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];

        assert_eq!(number.to_bytes_le(), expected_bytes);
    }

    #[test]
    fn from_bytes_be_works() {
        let bytes = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        ];
        let expected_number = U384::from_u64(1);

        assert_eq!(U384::from_bytes_be(&bytes).unwrap(), expected_number);
    }

    #[test]
    fn from_bytes_le_works() {
        let bytes = [
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let expected_number = U384::from_u64(1);

        assert_eq!(U384::from_bytes_le(&bytes).unwrap(), expected_number);
    }

    #[test]
    fn from_bytes_be_works_with_extra_data() {
        let bytes = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let expected_number = U384::from_u64(0);

        assert_eq!(U384::from_bytes_be(&bytes).unwrap(), expected_number);
    }

    #[test]
    #[should_panic]
    fn from_bytes_be_errs_with_less_data() {
        let bytes = [0, 0, 0, 0, 0];
        U384::from_bytes_be(&bytes).unwrap();
    }

    #[test]
    #[should_panic]
    fn from_bytes_le_errs_with_less_data() {
        let bytes = [0, 0, 0, 0, 0];
        U384::from_bytes_le(&bytes).unwrap();
    }

    #[test]
    fn from_bytes_le_works_with_extra_data() {
        let bytes = [
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
        ];
        let expected_number = U384::from_u64(1);

        assert_eq!(U384::from_bytes_le(&bytes).unwrap(), expected_number);
    }

    #[test]
    fn test_square_0() {
        let a = U384::from_hex_unchecked("362e35606447fb568704026c25da7a304bc7bd0aea36a61d77d4151395078cfa332b9d4928a60721eece725bbc81e158");
        let (hi, lo) = U384::square(&a);
        assert_eq!(lo, U384::from_hex_unchecked("11724caeb10c4bce5319097d74aed2246e2942b56b7365b5b2f8ceb3bb847db4828862043299d798577996e210bce40"));
        assert_eq!(hi, U384::from_hex_unchecked("b7786dbe41375b7ff64dbdc65152ef7d3fdbf499485e26486201cdbfb71b5673c77eb355a1274d08cbfbc1a4cdfdfad"));
    }

    #[test]
    fn test_square_1() {
        let a = U384::from_limbs([0, 0, 0, 0, 0, u64::MAX]);
        let (hi, lo) = U384::square(&a);
        assert_eq!(
            lo,
            U384::from_hex_unchecked("fffffffffffffffe0000000000000001")
        );
        assert_eq!(hi, U384::from_hex_unchecked("0"));
    }

    #[test]
    fn test_square_2() {
        let a = U384::from_limbs([0, 0, 0, 0, u64::MAX, 0]);
        let (hi, lo) = U384::square(&a);
        assert_eq!(
            lo,
            U384::from_hex_unchecked(
                "fffffffffffffffe000000000000000100000000000000000000000000000000"
            )
        );
        assert_eq!(hi, U384::from_hex_unchecked("0"));
    }

    #[test]
    fn test_square_3() {
        let a = U384::from_limbs([0, 0, 0, u64::MAX, 0, 0]);
        let (hi, lo) = U384::square(&a);
        assert_eq!(lo, U384::from_hex_unchecked("fffffffffffffffe00000000000000010000000000000000000000000000000000000000000000000000000000000000"));
        assert_eq!(hi, U384::from_hex_unchecked("0"));
    }

    #[test]
    fn test_square_4() {
        let a = U384::from_limbs([0, 0, u64::MAX, 0, 0, 0]);
        let (hi, lo) = U384::square(&a);
        assert_eq!(lo, U384::from_hex_unchecked("0"));
        assert_eq!(
            hi,
            U384::from_hex_unchecked("fffffffffffffffe0000000000000001")
        );
    }

    #[test]
    fn test_square_5() {
        let a = U384::from_limbs([0, 0, u64::MAX, u64::MAX, u64::MAX, u64::MAX]);
        let (hi, lo) = U384::square(&a);
        assert_eq!(lo, U384::from_hex_unchecked("fffffffffffffffffffffffffffffffe0000000000000000000000000000000000000000000000000000000000000001"));
        assert_eq!(
            hi,
            U384::from_hex_unchecked("ffffffffffffffffffffffffffffffff")
        );
    }

    #[test]
    fn test_square_6() {
        let a = U384::from_limbs([0, u64::MAX, u64::MAX, u64::MAX, u64::MAX, u64::MAX]);
        let (hi, lo) = U384::square(&a);
        assert_eq!(lo, U384::from_hex_unchecked("fffffffffffffffe00000000000000000000000000000000000000000000000000000000000000000000000000000001"));
        assert_eq!(
            hi,
            U384::from_hex_unchecked(
                "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
            )
        );
    }

    #[test]
    fn test_square_7() {
        let a = U384::from_limbs([u64::MAX, u64::MAX, u64::MAX, u64::MAX, u64::MAX, u64::MAX]);
        let (hi, lo) = U384::square(&a);
        assert_eq!(lo, U384::from_hex_unchecked("1"));
        assert_eq!(hi, U384::from_hex_unchecked("fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffe"));
    }
}

#[cfg(test)]
mod tests_u256 {
    use crate::unsigned_integer::element::ByteConversion;
    use crate::unsigned_integer::element::{UnsignedInteger, U256};
    #[cfg(feature = "proptest")]
    use proptest::prelude::*;
    #[cfg(feature = "proptest")]
    use std::ops::Shr;

    #[cfg(feature = "proptest")]
    const N_LIMBS: usize = 4;
    #[cfg(feature = "proptest")]
    type Uint = UnsignedInteger<N_LIMBS>;

    #[cfg(feature = "proptest")]
    proptest! {
        #[test]
        fn bitand(a in any::<Uint>(), b in any::<Uint>()) {
            let result = a & b;

            for i in 0..N_LIMBS {
                assert_eq!(result.limbs[i], a.limbs[i] & b.limbs[i]);
            }
        }

        #[test]
        fn bitand_assign(a in any::<Uint>(), b in any::<Uint>()) {
            let mut result = a;
            result &= b;

            for i in 0..N_LIMBS {
                assert_eq!(result.limbs[i], a.limbs[i] & b.limbs[i]);
            }
        }

        #[test]
        fn bitor(a in any::<Uint>(), b in any::<Uint>()) {
            let result = a | b;

            for i in 0..N_LIMBS {
                assert_eq!(result.limbs[i], a.limbs[i] | b.limbs[i]);
            }
        }

        #[test]
        fn bitor_assign(a in any::<Uint>(), b in any::<Uint>()) {
            let mut result = a;
            result |= b;

            for i in 0..N_LIMBS {
                assert_eq!(result.limbs[i], a.limbs[i] | b.limbs[i]);
            }
        }

        #[test]
        fn bitxor(a in any::<Uint>(), b in any::<Uint>()) {
            let result = a ^ b;

            for i in 0..N_LIMBS {
                assert_eq!(result.limbs[i], a.limbs[i] ^ b.limbs[i]);
            }
        }

        #[test]
        fn bitxor_assign(a in any::<Uint>(), b in any::<Uint>()) {
            let mut result = a;
            result ^= b;

            for i in 0..N_LIMBS {
                assert_eq!(result.limbs[i], a.limbs[i] ^ b.limbs[i]);
            }
        }

        #[test]
        fn div_rem(a in any::<Uint>(), b in any::<Uint>()) {
            let a = a.shr(128);
            let b = b.shr(128);
            assert_eq!((a * b).div_rem(&b), (a, Uint::from_u64(0)));
        }
    }

    #[test]
    fn construct_new_integer_from_limbs() {
        let a: U256 = UnsignedInteger {
            limbs: [0, 1, 2, 3],
        };
        assert_eq!(U256::from_limbs([0, 1, 2, 3]), a);
    }

    #[test]
    fn construct_new_integer_from_u64_1() {
        let a = U256::from_u64(1_u64);
        assert_eq!(a.limbs, [0, 0, 0, 1]);
    }

    #[test]
    fn construct_new_integer_from_u64_2() {
        let a = U256::from_u64(u64::MAX);
        assert_eq!(a.limbs, [0, 0, 0, u64::MAX]);
    }

    #[test]
    fn construct_new_integer_from_u128_1() {
        let a = U256::from_u128(u128::MAX);
        assert_eq!(a.limbs, [0, 0, u64::MAX, u64::MAX]);
    }

    #[test]
    fn construct_new_integer_from_u128_4() {
        let a = U256::from_u128(276371540478856090688472252609570374439);
        assert_eq!(a.limbs, [0, 0, 14982131230017065096, 14596400355126379303]);
    }

    #[test]
    fn construct_new_integer_from_hex_1() {
        let a = U256::from_hex_unchecked("1");
        assert_eq!(a.limbs, [0, 0, 0, 1]);
    }

    #[test]
    fn construct_new_integer_from_hex_2() {
        let a = U256::from_hex_unchecked("f");
        assert_eq!(a.limbs, [0, 0, 0, 15]);
    }

    #[test]
    fn construct_new_integer_from_hex_3() {
        let a = U256::from_hex_unchecked("10000000000000000");
        assert_eq!(a.limbs, [0, 0, 1, 0]);
    }

    #[test]
    fn construct_new_integer_from_hex_4() {
        let a = U256::from_hex_unchecked("a0000000000000000");
        assert_eq!(a.limbs, [0, 0, 10, 0]);
    }

    #[test]
    fn construct_new_integer_from_hex_5() {
        let a = U256::from_hex_unchecked("ffffffffffffffffff");
        assert_eq!(a.limbs, [0, 0, 255, u64::MAX]);
    }

    #[test]
    fn construct_new_integer_from_hex_6() {
        let a = U256::from_hex_unchecked("eb235f6144d9e91f4b14");
        assert_eq!(a.limbs, [0, 0, 60195, 6872850209053821716]);
    }

    #[test]
    fn construct_new_integer_from_hex_7() {
        let a = U256::from_hex_unchecked("2b20aaa5cf482b239e2897a787faf4660cc95597854beb2");
        assert_eq!(
            a.limbs,
            [
                0,
                194229460750598834,
                4171047363999149894,
                6975114134393503410
            ]
        );
    }

    #[test]
    fn construct_new_integer_from_hex_8() {
        let a = U256::from_hex_unchecked(
            "2B20AAA5CF482B239E2897A787FAF4660CC95597854BEB235F6144D9E91F4B14",
        );
        assert_eq!(
            a.limbs,
            [
                3107671372009581347,
                11396525602857743462,
                921361708038744867,
                6872850209053821716
            ]
        );
    }

    #[test]
    fn construct_new_integer_from_dec_1() {
        let a = U256::from_dec_str("1").unwrap();
        assert_eq!(a.limbs, [0, 0, 0, 1]);
    }

    #[test]
    fn construct_integer_from_invalid_hex_returns_error() {
        use crate::unsigned_integer::element::CreationError;
        assert_eq!(U256::from_hex("0xaO"), Err(CreationError::InvalidHexString));
        assert_eq!(U256::from_hex("0xOa"), Err(CreationError::InvalidHexString));
        assert_eq!(U256::from_hex("0xm"), Err(CreationError::InvalidHexString));
    }

    #[test]
    fn construct_new_integer_from_dec_2() {
        let a = U256::from_dec_str("15").unwrap();
        assert_eq!(a.limbs, [0, 0, 0, 15]);
    }

    #[test]
    fn construct_new_integer_from_dec_3() {
        let a = U256::from_dec_str("18446744073709551616").unwrap();
        assert_eq!(a.limbs, [0, 0, 1, 0]);
    }

    #[test]
    fn construct_new_integer_from_dec_4() {
        let a = U256::from_dec_str("184467440737095516160").unwrap();
        assert_eq!(a.limbs, [0, 0, 10, 0]);
    }

    #[test]
    fn construct_new_integer_from_dec_5() {
        let a = U256::from_dec_str("4722366482869645213695").unwrap();
        assert_eq!(a.limbs, [0, 0, 255, u64::MAX]);
    }

    #[test]
    fn construct_new_integer_from_dec_6() {
        let a = U256::from_dec_str("1110408632367155513346836").unwrap();
        assert_eq!(a.limbs, [0, 0, 60195, 6872850209053821716]);
    }

    #[test]
    fn construct_new_integer_from_dec_7() {
        let a =
            U256::from_dec_str("66092860629991288370279803883558073888453977263446474418").unwrap();
        assert_eq!(
            a.limbs,
            [
                0,
                194229460750598834,
                4171047363999149894,
                6975114134393503410
            ]
        );
    }

    #[test]
    fn construct_new_integer_from_dec_8() {
        let a = U256::from_dec_str(
            "19507169362252850253634654373914901165934018806002526957372506333098895428372",
        )
        .unwrap();
        assert_eq!(
            a.limbs,
            [
                3107671372009581347,
                11396525602857743462,
                921361708038744867,
                6872850209053821716
            ]
        );
    }

    #[test]
    fn construct_new_integer_from_dec_empty() {
        assert!(U256::from_dec_str("").is_err());
    }

    #[test]
    fn construct_new_integer_from_dec_invalid() {
        assert!(U256::from_dec_str("0xff").is_err());
    }

    #[test]
    fn equality_works_1() {
        let a = U256::from_hex_unchecked("1");
        let b = U256 {
            limbs: [0, 0, 0, 1],
        };
        assert_eq!(a, b);
    }
    #[test]
    fn equality_works_2() {
        let a = U256::from_hex_unchecked("f");
        let b = U256 {
            limbs: [0, 0, 0, 15],
        };
        assert_eq!(a, b);
    }

    #[test]
    fn equality_works_3() {
        let a = U256::from_hex_unchecked("10000000000000000");
        let b = U256 {
            limbs: [0, 0, 1, 0],
        };
        assert_eq!(a, b);
    }

    #[test]
    fn equality_works_4() {
        let a = U256::from_hex_unchecked("a0000000000000000");
        let b = U256 {
            limbs: [0, 0, 10, 0],
        };
        assert_eq!(a, b);
    }

    #[test]
    fn equality_works_5() {
        let a = U256::from_hex_unchecked("ffffffffffffffffff");
        let b = U256 {
            limbs: [0, 0, u8::MAX as u64, u64::MAX],
        };
        assert_eq!(a, b);
    }

    #[test]
    fn equality_works_6() {
        let a = U256::from_hex_unchecked("eb235f6144d9e91f4b14");
        let b = U256 {
            limbs: [0, 0, 60195, 6872850209053821716],
        };
        assert_eq!(a, b);
    }

    #[test]
    fn equality_works_7() {
        let a = U256::from_hex_unchecked("2b20aaa5cf482b239e2897a787faf4660cc95597854beb2");
        let b = U256 {
            limbs: [
                0,
                194229460750598834,
                4171047363999149894,
                6975114134393503410,
            ],
        };
        assert_eq!(a, b);
    }

    #[test]
    fn equality_works_8() {
        let a = U256::from_hex_unchecked(
            "2B20AAA5CF482B239E2897A787FAF4660CC95597854BEB235F6144D9E91F4B14",
        );
        let b = U256 {
            limbs: [
                3107671372009581347,
                11396525602857743462,
                921361708038744867,
                6872850209053821716,
            ],
        };
        assert_eq!(a, b);
    }

    #[test]
    fn equality_works_9() {
        let a = U256::from_hex_unchecked("fffffff");
        let b = U256::from_hex_unchecked("fefffff");
        assert_ne!(a, b);
    }

    #[test]
    fn equality_works_10() {
        let a = U256::from_hex_unchecked("ffff000000000000");
        let b = U256::from_hex_unchecked("ffff000000100000");
        assert_ne!(a, b);
    }

    #[test]
    fn double_256_bit_integer_1() {
        let a = U256::from_u64(2);
        let b = U256::from_u64(5);
        let c = U256::from_u64(7);
        assert_eq!(U256::double(&a).0, a + a);
        assert_eq!(U256::double(&b).0, b + b);
        assert_eq!(U256::double(&c).0, c + c);
    }

    #[test]
    fn add_two_256_bit_integers_1() {
        let a = U256::from_u64(2);
        let b = U256::from_u64(5);
        let c = U256::from_u64(7);
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_256_bit_integers_2() {
        let a = U256::from_u64(334);
        let b = U256::from_u64(666);
        let c = U256::from_u64(1000);
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_256_bit_integers_3() {
        let a = U256::from_hex_unchecked("ffffffffffffffff");
        let b = U256::from_hex_unchecked("1");
        let c = U256::from_hex_unchecked("10000000000000000");
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_256_bit_integers_4() {
        let a = U256::from_hex_unchecked("b58e1e0b66");
        let b = U256::from_hex_unchecked("55469d9619");
        let c = U256::from_hex_unchecked("10ad4bba17f");
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_256_bit_integers_5() {
        let a = U256::from_hex_unchecked("e8dff25cb6160f7705221da6f");
        let b = U256::from_hex_unchecked("ab879169b5f80dc8a7969f0b0");
        let c = U256::from_hex_unchecked("1946783c66c0e1d3facb8bcb1f");
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_256_bit_integers_6() {
        let a = U256::from_hex_unchecked("9adf291af3a64d59e14e7b440c850508014c551ed5");
        let b = U256::from_hex_unchecked("e7948474bce907f0feaf7e5d741a8cd2f6d1fb9448");
        let c = U256::from_hex_unchecked("18273ad8fb08f554adffdf9a1809f91daf81e50b31d");
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_256_bit_integers_7() {
        let a = U256::from_hex_unchecked(
            "10d3bc05496380cfe27bf5d97ddb99ac95eb5ecfbd3907eadf877a4c2dfa05f6",
        );
        let b = U256::from_hex_unchecked(
            "0866aef803c92bf02e85c7fad0eccb4881c59825e499fa22f98e1a8fefed4cd9",
        );
        let c = U256::from_hex_unchecked(
            "193a6afd4d2cacc01101bdd44ec864f517b0f6f5a1d3020dd91594dc1de752cf",
        );
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_256_bit_integers_8() {
        let a = U256::from_hex_unchecked(
            "07df9c74fa9d5aafa74a87dbbf93215659d8a3e1706d4b06de9512284802580f",
        );
        let b = U256::from_hex_unchecked(
            "d515e54973f0643a6a9957579c1f84020a6a91d5d5f27b75401c7538d2c9ea9c",
        );
        let c = U256::from_hex_unchecked(
            "dcf581be6e8dbeea11e3df335bb2a558644335b7465fc67c1eb187611acc42ab",
        );
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_256_bit_integers_9() {
        let a = U256::from_hex_unchecked(
            "92977527a0f8ba00d18c1b2f1900d965d4a70e5f5f54468ffb2d4d41519385f2",
        );
        let b = U256::from_hex_unchecked(
            "46facf9953a9494822bf18836ffd7e55c48b30aa81e17fa1ace0b473015307e4",
        );
        let c = U256::from_hex_unchecked(
            "d99244c0f4a20348f44b33b288fe57bb99323f09e135c631a80e01b452e68dd6",
        );
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_two_256_bit_integers_10() {
        let a = U256::from_hex_unchecked(
            "07df9c74fa9d5aafa74a87dbbf93215659d8a3e1706d4b06de9512284802580f",
        );
        let b = U256::from_hex_unchecked(
            "d515e54973f0643a6a9957579c1f84020a6a91d5d5f27b75401c7538d2c9ea9c",
        );
        let c_expected = U256::from_hex_unchecked(
            "dcf581be6e8dbeea11e3df335bb2a558644335b7465fc67c1eb187611acc42ab",
        );
        let (c, overflow) = U256::add(&a, &b);
        assert_eq!(c, c_expected);
        assert!(!overflow);
    }

    #[test]
    fn add_two_256_bit_integers_11() {
        let a = U256::from_hex_unchecked(
            "92977527a0f8ba00d18c1b2f1900d965d4a70e5f5f54468ffb2d4d41519385f2",
        );
        let b = U256::from_hex_unchecked(
            "46facf9953a9494822bf18836ffd7e55c48b30aa81e17fa1ace0b473015307e4",
        );
        let c_expected = U256::from_hex_unchecked(
            "d99244c0f4a20348f44b33b288fe57bb99323f09e135c631a80e01b452e68dd6",
        );
        let (c, overflow) = U256::add(&a, &b);
        assert_eq!(c, c_expected);
        assert!(!overflow);
    }

    #[test]
    fn add_two_256_bit_integers_12_with_overflow() {
        let a = U256::from_hex_unchecked(
            "b07bc844363dd56467d9ebdd5929e9bb34a8e2577db77df6cf8f2ac45bd3d0bc",
        );
        let b = U256::from_hex_unchecked(
            "cbbc474761bb7995ff54e25fa5d30295604fe3545d0cde405e72d8c0acebb119",
        );
        let c_expected = U256::from_hex_unchecked(
            "7c380f8b97f94efa672ece3cfefcec5094f8c5abdac45c372e02038508bf81d5",
        );
        let (c, overflow) = U256::add(&a, &b);
        assert_eq!(c, c_expected);
        assert!(overflow);
    }

    #[test]
    fn double_256_bit_integer_12_with_overflow() {
        let a = U256::from_hex_unchecked(
            "b07bc844363dd56467d9ebdd5929e9bb34a8e2577db77df6cf8f2ac45bd3d0bc",
        );
        let b = U256::from_hex_unchecked(
            "cbbc474761bb7995ff54e25fa5d30295604fe3545d0cde405e72d8c0acebb119",
        );
        assert_eq!(U256::double(&a), U256::add(&a, &a));
        assert_eq!(U256::double(&b), U256::add(&b, &b));
    }

    #[test]
    fn sub_two_256_bit_integers_1() {
        let a = U256::from_u64(2);
        let b = U256::from_u64(5);
        let c = U256::from_u64(7);
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_256_bit_integers_2() {
        let a = U256::from_u64(334);
        let b = U256::from_u64(666);
        let c = U256::from_u64(1000);
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_256_bit_integers_3() {
        let a = U256::from_hex_unchecked("ffffffffffffffff");
        let b = U256::from_hex_unchecked("1");
        let c = U256::from_hex_unchecked("10000000000000000");
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_256_bit_integers_4() {
        let a = U256::from_hex_unchecked("b58e1e0b66");
        let b = U256::from_hex_unchecked("55469d9619");
        let c = U256::from_hex_unchecked("10ad4bba17f");
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_256_bit_integers_5() {
        let a = U256::from_hex_unchecked("e8dff25cb6160f7705221da6f");
        let b = U256::from_hex_unchecked("ab879169b5f80dc8a7969f0b0");
        let c = U256::from_hex_unchecked("1946783c66c0e1d3facb8bcb1f");
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_256_bit_integers_6() {
        let a = U256::from_hex_unchecked("9adf291af3a64d59e14e7b440c850508014c551ed5");
        let b = U256::from_hex_unchecked("e7948474bce907f0feaf7e5d741a8cd2f6d1fb9448");
        let c = U256::from_hex_unchecked("18273ad8fb08f554adffdf9a1809f91daf81e50b31d");
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_256_bit_integers_7() {
        let a = U256::from_hex_unchecked(
            "9b4000dccf01a010e196154a1b998408f949d734389626ba97cb3331ee87e01d",
        );
        let b = U256::from_hex_unchecked(
            "5d26ae1b34c78bdf4cefb2b0b553473f887bc0f1ac03d36861c2e75e01656cbc",
        );
        let c = U256::from_hex_unchecked(
            "f866aef803c92bf02e85c7fad0eccb4881c59825e499fa22f98e1a8fefed4cd9",
        );
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_256_bit_integers_8() {
        let a = U256::from_hex_unchecked(
            "07df9c74fa9d5aafa74a87dbbf93215659d8a3e1706d4b06de9512284802580e",
        );
        let b = U256::from_hex_unchecked(
            "d515e54973f0643a6a9957579c1f84020a6a91d5d5f27b75401c7538d2c9ea9d",
        );
        let c = U256::from_hex_unchecked(
            "dcf581be6e8dbeea11e3df335bb2a558644335b7465fc67c1eb187611acc42ab",
        );
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_256_bit_integers_9() {
        let a = U256::from_hex_unchecked(
            "92977527a0f8ba00d18c1b2f1900d965d4a70e5f5f54468ffb2d4d41519385f2",
        );
        let b = U256::from_hex_unchecked(
            "46facf9953a9494822bf18836ffd7e55c48b30aa81e17fa1ace0b473015307e4",
        );
        let c = U256::from_hex_unchecked(
            "d99244c0f4a20348f44b33b288fe57bb99323f09e135c631a80e01b452e68dd6",
        );
        assert_eq!(c - a, b);
    }

    #[test]
    fn sub_two_256_bit_integers_11_without_overflow() {
        let a = U256::from_u64(334);
        let b_expected = U256::from_u64(666);
        let c = U256::from_u64(1000);
        let (b, overflow) = U256::sub(&c, &a);
        assert!(!overflow);
        assert_eq!(b_expected, b);
    }

    #[test]
    fn sub_two_256_bit_integers_11_with_overflow() {
        let a = U256::from_u64(334);
        let b_expected = U256::from_hex_unchecked(
            "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffd66",
        );
        let c = U256::from_u64(1000);
        let (b, overflow) = U256::sub(&a, &c);
        assert!(overflow);
        assert_eq!(b_expected, b);
    }

    #[test]
    fn const_le_works() {
        let a = U256::from_u64(334);
        let b = U256::from_u128(333);
        assert!(U256::const_le(&b, &a));
        assert!(U256::const_le(&a, &a));
        assert!(!U256::const_le(&a, &b));
    }

    #[test]
    fn partial_order_works() {
        assert!(U256::from_u64(10) <= U256::from_u64(10));
        assert!(U256::from_u64(1) < U256::from_u64(2));
        assert!(U256::from_u64(2) >= U256::from_u64(1));

        assert!(U256::from_u64(10) >= U256::from_u64(10));
        assert!(U256::from_u64(2) > U256::from_u64(1));
        assert!(U256::from_u64(1) <= U256::from_u64(2));

        let a = U256::from_hex_unchecked(
            "5d4a70e5f5f54468ffb2d4d41519385f24b078a0e7d0281d5ad0c36724dc4233",
        );
        let c = U256::from_hex_unchecked(
            "b99323f09e135c631a80e01b452e68dd6ad3315e5776b713af4a2d7f5f9a2a75",
        );

        assert!(&a <= &a);
        assert!(&a >= &a);
        assert!(&a >= &a);
        assert!(&a <= &a);
        assert!(&a < &(&a + U256::from_u64(1)));
        assert!(&a <= &(&a + U256::from_u64(1)));
        assert!(&a + U256::from_u64(1) > a);
        assert!((&a + U256::from_u64(1) >= a));
        assert!(&a <= &c);
        assert!(&a < &c);
        assert!(&a < &c);
        assert!(&a <= &c);
        assert!(&c > &a);
        assert!(&c >= &a);
        assert!(&c >= &a);
        assert!(&c > &a);
        assert!(a < c);
    }

    #[test]
    fn mul_two_256_bit_integers_works_1() {
        let a = U256::from_u64(3);
        let b = U256::from_u64(8);
        let c = U256::from_u64(3 * 8);
        assert_eq!(a * b, c);
    }

    #[test]
    fn mul_two_256_bit_integers_works_2() {
        let a = U256::from_hex_unchecked("6131d99f840b3b0");
        let b = U256::from_hex_unchecked("6f5c466db398f43");
        let c = U256::from_hex_unchecked("2a47a603a77f871dfbb937af7e5710");
        assert_eq!(a * b, c);
    }

    #[test]
    fn mul_two_256_bit_integers_works_3() {
        let a = U256::from_hex_unchecked("84a6add5db9e095b2e0f6b40eff8ee");
        let b = U256::from_hex_unchecked("2347db918f725461bec2d5c57");
        let c = U256::from_hex_unchecked("124805c476c9462adc0df6c88495d4253f5c38033afc18d78d920e2");
        assert_eq!(a * b, c);
    }

    #[test]
    fn mul_two_256_bit_integers_works_4() {
        let a = U256::from_hex_unchecked("15bf61fcf53a3f0ae1e8e555d");
        let b = U256::from_hex_unchecked("cbbc474761bb7995ff54e25fa5d5d0cde405e9f");
        let c_expected = U256::from_hex_unchecked(
            "114ec14db0c80d30b7dcb9c45948ef04cc149e612cb544f447b146553aff2ac3",
        );
        assert_eq!(a * b, c_expected);
    }

    #[test]
    fn mul_two_256_bit_integers_works_5_hi_lo() {
        let a = U256::from_hex_unchecked(
            "8e2d939b602a50911232731d04fe6f40c05f97da0602307099fb991f9b414e2d",
        );
        let b = U256::from_hex_unchecked(
            "7f3ad1611ab58212f92a2484e9560935b9ac4615fe61cfed1a4861e193a74d20",
        );
        let hi_expected = U256::from_hex_unchecked(
            "46A946D6A984FE6507DE6B8D1354256D7A7BAE4283404733BDC876A264BCE5EE",
        );
        let lo_expected = U256::from_hex_unchecked(
            "43F24263F10930EBE3EA0307466C19B13B9C7DBA6B3F7604B7F32FB0E3084EA0",
        );
        let (hi, lo) = U256::mul(&a, &b);
        assert_eq!(hi, hi_expected);
        assert_eq!(lo, lo_expected);
    }

    #[test]
    fn shift_left_on_256_bit_integer_works_1() {
        let a = U256::from_hex_unchecked("1");
        let b = U256::from_hex_unchecked("10");
        assert_eq!(a << 4, b);
    }

    #[test]
    fn shift_left_on_256_bit_integer_works_2() {
        let a = U256::from_u64(1);
        let b = U256::from_u128(1_u128 << 64);
        assert_eq!(a << 64, b);
    }

    #[test]
    fn shift_left_on_256_bit_integer_works_3() {
        let a = U256::from_hex_unchecked("10");
        let b = U256::from_hex_unchecked("1000");
        assert_eq!(&a << 8, b);
    }

    #[test]
    fn shift_left_on_256_bit_integer_works_4() {
        let a = U256::from_hex_unchecked("e45542992b6844553f3cb1c5ac33e7fa5");
        let b = U256::from_hex_unchecked("391550a64ada11154fcf2c716b0cf9fe940");
        assert_eq!(a << 6, b);
    }

    #[test]
    fn shift_left_on_256_bit_integer_works_5() {
        let a = U256::from_hex_unchecked("a8390aa99bead76bc0093b1bc1a8101f5ce");
        let b = U256::from_hex_unchecked(
            "72155337d5aed7801276378350203eb9c0000000000000000000000000000000",
        );
        assert_eq!(&a << 125, b);
    }

    #[test]
    fn shift_left_on_256_bit_integer_works_6() {
        let a = U256::from_hex_unchecked("2ed786ab132f0b5b0cacd385dd51de3a");
        let b = U256::from_hex_unchecked(
            "2ed786ab132f0b5b0cacd385dd51de3a00000000000000000000000000000000",
        );
        assert_eq!(&a << (64 * 2), b);
    }

    #[test]
    fn shift_left_on_256_bit_integer_works_7() {
        let a = U256::from_hex_unchecked("90823e0bd707f");
        let b = U256::from_hex_unchecked(
            "90823e0bd707f000000000000000000000000000000000000000000000000",
        );
        assert_eq!(&a << (64 * 3), b);
    }

    #[test]
    fn shift_right_on_256_bit_integer_works_1() {
        let a = U256::from_hex_unchecked("1");
        let b = U256::from_hex_unchecked("10");
        assert_eq!(b >> 4, a);
    }

    #[test]
    fn shift_right_on_256_bit_integer_works_2() {
        let a = U256::from_hex_unchecked("10");
        let b = U256::from_hex_unchecked("1000");
        assert_eq!(&b >> 8, a);
    }

    #[test]
    fn shift_right_on_256_bit_integer_works_3() {
        let a = U256::from_hex_unchecked("e45542992b6844553f3cb1c5ac33e7fa5");
        let b = U256::from_hex_unchecked("391550a64ada11154fcf2c716b0cf9fe940");
        assert_eq!(b >> 6, a);
    }

    #[test]
    fn shift_right_on_256_bit_integer_works_4() {
        let a = U256::from_hex_unchecked("390aa99bead76bc0093b1bc1a8101f5ce");
        let b = U256::from_hex_unchecked(
            "72155337d5aed7801276378350203eb9c0000000000000000000000000000000",
        );
        assert_eq!(&b >> 125, a);
    }

    #[test]
    fn shift_right_on_256_bit_integer_works_5() {
        let a = U256::from_hex_unchecked(
            "ba6ab46f9a9a2f20e4061b67ce4d8c3da98091cf990d7b14ef47ffe27370abbd",
        );
        let b = U256::from_hex_unchecked("174d568df35345e41c80c36cf9c");
        assert_eq!(a >> 151, b);
    }

    #[test]
    fn shift_right_on_256_bit_integer_works_6() {
        let a = U256::from_hex_unchecked(
            "076c075d2f65e39b9ecdde8bf6f8c94241962ce0f557b7739673200c777152eb",
        );
        let b = U256::from_hex_unchecked("ed80eba5ecbc7373d9bbd17edf19284832c59c");
        assert_eq!(&a >> 99, b);
    }

    #[test]
    fn shift_right_on_256_bit_integer_works_7() {
        let a = U256::from_hex_unchecked(
            "6a9ce35d8940a5ebd29604ce9a182ade76f03f7e9965760b84a8cfd1d3dd2e61",
        );
        let b = U256::from_hex_unchecked("6a9ce35d8940a5eb");
        assert_eq!(&a >> (64 * 3), b);
    }

    #[test]
    fn shift_right_on_256_bit_integer_works_8() {
        let a = U256::from_hex_unchecked(
            "5322c128ec84081b6c376c108ebd7fd36bbd44f71ee5e6ad6bcb3dd1c5265bd7",
        );
        let b = U256::from_hex_unchecked("5322c128ec84081b6c376c108ebd7fd3");
        assert_eq!(&a >> (64 * 2), b);
    }

    #[test]
    #[cfg(feature = "alloc")]
    fn to_be_bytes_works() {
        let number = U256::from_u64(1);
        let expected_bytes = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1,
        ];

        assert_eq!(number.to_bytes_be(), expected_bytes);
    }

    #[test]
    #[cfg(feature = "alloc")]
    fn to_le_bytes_works() {
        let number = U256::from_u64(1);
        let expected_bytes = [
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0,
        ];

        assert_eq!(number.to_bytes_le(), expected_bytes);
    }

    #[test]
    fn from_bytes_be_works() {
        let bytes = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1,
        ];
        let expected_number = U256::from_u64(1);
        assert_eq!(U256::from_bytes_be(&bytes).unwrap(), expected_number);
    }

    #[test]
    fn from_bytes_le_works() {
        let bytes = [
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0,
        ];
        let expected_number = U256::from_u64(1);
        assert_eq!(U256::from_bytes_le(&bytes).unwrap(), expected_number);
    }

    #[test]
    fn shr_inplace_works_1() {
        let mut n = UnsignedInteger::<3>::from(4u64);
        n >>= 1;

        assert_eq!(n, UnsignedInteger::<3>::from(2u64));
    }

    #[test]
    fn shr_inplace_on_256_bit_integer_works_1() {
        let a = U256::from_hex_unchecked("e45542992b6844553f3cb1c5ac33e7fa5");
        let mut b = U256::from_hex_unchecked("391550a64ada11154fcf2c716b0cf9fe940");
        b >>= 6;
        assert_eq!(a, b);
    }

    #[test]
    fn shr_inplace_on_254_bit_integer_works_2() {
        let a = U256::from_hex_unchecked("390aa99bead76bc0093b1bc1a8101f5ce");
        let mut b = U256::from_hex_unchecked(
            "72155337d5aed7801276378350203eb9c0000000000000000000000000000000",
        );
        b >>= 125;
        assert_eq!(a, b);
    }

    #[test]
    fn shr_inplace_on_256_bit_integer_works_3() {
        let a = U256::from_hex_unchecked("2ed786ab132f0b5b0cacd385dd51de3a");
        let mut b = U256::from_hex_unchecked(
            "2ed786ab132f0b5b0cacd385dd51de3a00000000000000000000000000000000",
        );
        b >>= 64 * 2;
        assert_eq!(a, b);
    }

    #[test]
    fn shr_inplace_on_256_bit_integer_works_4() {
        let a = U256::from_hex_unchecked("90823e0bd707f");
        let mut b = U256::from_hex_unchecked(
            "90823e0bd707f000000000000000000000000000000000000000000000000",
        );
        b >>= 64 * 3;
        assert_eq!(a, b);
    }

    #[test]
    fn shr_inplace_on_256_bit_integer_works_5() {
        let a = U256::from_hex_unchecked("24208f");
        let mut b = U256::from_hex_unchecked(
            "90823e0bd707f000000000000000000000000000000000000000000000000",
        );
        b >>= 222;
        assert_eq!(a, b);
    }

    #[test]
    fn multiplying_and_dividing_for_number_is_number_with_remainder_0() {
        let a = U256::from_u128(12678920202929299999999999282828);
        let b = U256::from_u128(9000000000000);
        assert_eq!((a * b).div_rem(&b), (a, U256::from_u64(0)));
    }

    #[test]
    fn unsigned_int_8_div_rem_3_is_2_2() {
        let a: UnsignedInteger<4> = U256::from_u64(8);
        let b = U256::from_u64(3);
        assert_eq!(a.div_rem(&b), (U256::from_u64(2), U256::from_u64(2)));
    }

    #[test]
    fn unsigned_int_500721_div_rem_5_is_100144_1() {
        let a = U256::from_u64(500721);
        let b = U256::from_u64(5);
        assert_eq!(a.div_rem(&b), (U256::from_u64(100144), U256::from_u64(1)));
    }

    #[test]
    fn div_rem_works_with_big_numbers() {
        let a = U256::from_u128(4758402376589578934275873583589345);
        let b = U256::from_u128(43950384634609);
        assert_eq!(
            a.div_rem(&b),
            (
                U256::from_u128(108267593472721187331),
                U256::from_u128(12368508650766)
            )
        );
    }

    #[cfg(feature = "std")]
    #[test]
    fn to_hex_test() {
        let a = U256::from_hex_unchecked("390aa99bead76bc0093b1bc1a8101f5ce");
        assert_eq!(U256::to_hex(&a), "390AA99BEAD76BC0093B1BC1A8101F5CE")
    }
}
