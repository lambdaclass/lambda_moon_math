use std::io;

use lambdaworks_math::errors::DeserializationError;

#[derive(Debug)]
pub enum SrsFromFileError {
    FileError(io::Error),
    DeserializationError(lambdaworks_math::errors::DeserializationError),
}

impl From<lambdaworks_math::errors::DeserializationError> for SrsFromFileError {
    fn from(err: DeserializationError) -> SrsFromFileError {
        match err {
            DeserializationError::InvalidAmountOfBytes => {
                SrsFromFileError::DeserializationError(DeserializationError::InvalidAmountOfBytes)
            }

            DeserializationError::FieldFromBytesError => {
                SrsFromFileError::DeserializationError(DeserializationError::FieldFromBytesError)
            }

            DeserializationError::PointerSizeError => {
                SrsFromFileError::DeserializationError(DeserializationError::PointerSizeError)
            }

            DeserializationError::InvalidValue => {
                SrsFromFileError::DeserializationError(DeserializationError::InvalidValue)
            }
        }
    }
}

impl From<std::io::Error> for SrsFromFileError {
    fn from(err: std::io::Error) -> SrsFromFileError {
        SrsFromFileError::FileError(err)
    }
}

#[derive(Debug)]
pub enum ProverVerifyKeysFromFileError {
    FileError(io::Error),
    DeserializationError(lambdaworks_math::errors::DeserializationError),
}

impl From<lambdaworks_math::errors::DeserializationError> for ProverVerifyKeysFromFileError {
    fn from(err: DeserializationError) -> ProverVerifyKeysFromFileError {
        match err {
            DeserializationError::InvalidAmountOfBytes => {
                ProverVerifyKeysFromFileError::DeserializationError(
                    DeserializationError::InvalidAmountOfBytes,
                )
            }

            DeserializationError::FieldFromBytesError => {
                ProverVerifyKeysFromFileError::DeserializationError(
                    DeserializationError::FieldFromBytesError,
                )
            }

            DeserializationError::PointerSizeError => {
                ProverVerifyKeysFromFileError::DeserializationError(
                    DeserializationError::PointerSizeError,
                )
            }

            DeserializationError::InvalidValue => {
                ProverVerifyKeysFromFileError::DeserializationError(
                    DeserializationError::InvalidValue,
                )
            }
        }
    }
}

impl From<std::io::Error> for ProverVerifyKeysFromFileError {
    fn from(err: std::io::Error) -> ProverVerifyKeysFromFileError {
        ProverVerifyKeysFromFileError::FileError(err)
    }
}
