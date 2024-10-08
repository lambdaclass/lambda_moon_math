#[derive(Debug, PartialEq, Eq)]
pub enum ByteConversionError {
    FromBEBytesError,
    FromLEBytesError,
    InvalidValue,
    PointNotInSubgroup,
    ValueNotCompressed,
    InvalidLength,
}

#[derive(Debug, PartialEq, Eq)]
pub enum CreationError {
    InvalidHexString,
    InvalidDecString,
    HexStringIsTooBig,
    EmptyString,
}

#[derive(Debug, PartialEq, Eq)]
pub enum DeserializationError {
    InvalidAmountOfBytes,
    FieldFromBytesError,
    PointerSizeError,
    InvalidValue,
}

#[derive(Debug, PartialEq, Eq)]
pub enum PairingError {
    PointNotInSubgroup,
}

impl From<ByteConversionError> for DeserializationError {
    fn from(error: ByteConversionError) -> Self {
        match error {
            ByteConversionError::FromBEBytesError => DeserializationError::FieldFromBytesError,
            ByteConversionError::FromLEBytesError => DeserializationError::FieldFromBytesError,
            _ => DeserializationError::InvalidValue,
        }
    }
}
