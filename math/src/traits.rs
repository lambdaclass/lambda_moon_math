use crate::{
    errors::ByteConversionError,
    field::{element::FieldElement, traits::IsField},
};

/// A trait for converting an element to and from its byte representation and
/// for getting an element from its byte representation in big-endian or
/// little-endian order.
pub trait ByteConversion {
    /// Returns the byte representation of the element in big-endian order.
    fn to_bytes_be(&self) -> Vec<u8>;

    /// Returns the byte representation of the element in little-endian order.
    fn to_bytes_le(&self) -> Vec<u8>;

    /// Returns the element from its byte representation in big-endian order.
    fn from_bytes_be(bytes: &[u8]) -> Result<Self, ByteConversionError>
    where
        Self: std::marker::Sized;

    /// Returns the element from its byte representation in little-endian order.
    fn from_bytes_le(bytes: &[u8]) -> Result<Self, ByteConversionError>
    where
        Self: std::marker::Sized;
}

/// Serialize function without options
/// Used to serialize data to feed it into Fiat Shamir
/// in some protocols, when ByteConversion is not avaible
/// For example, when using Curve Points
pub trait Serializable {
    /// Returns the byte representation of the element in big-endian order.
    fn serialize(&self) -> Vec<u8>;
}

pub trait IsRandomFieldElementGenerator<F: IsField> {
    fn generate(&self) -> FieldElement<F>;
}
