use crate::{
    cyclic_group::IsGroup,
    errors::PairingError,
    field::{element::FieldElement, traits::IsField},
};
use core::fmt::Debug;

#[derive(Debug, PartialEq, Eq)]
pub enum EllipticCurveError {
    InvalidPoint,
}

pub trait IsEllipticCurve {
    /// BaseField is the field used for each of the coordinates of a point p
    /// belonging to the curve.
    type BaseField: IsField + Clone + Debug;

    /// The representation of the point. For example it can be projective
    /// coordinates, affine coordinates, XYZZ, depending on the curve and its
    /// possible optimizations.
    type PointRepresentation: IsGroup + FromAffine<Self::BaseField>;

    /// Returns the generator of the main subgroup.
    fn generator() -> Self::PointRepresentation;

    /// Returns an affine point.
    fn create_point_from_affine(
        x: FieldElement<Self::BaseField>,
        y: FieldElement<Self::BaseField>,
    ) -> Result<Self::PointRepresentation, EllipticCurveError> {
        Self::PointRepresentation::from_affine(x, y)
    }
}

pub trait FromAffine<F: IsField>: Sized {
    fn from_affine(x: FieldElement<F>, y: FieldElement<F>) -> Result<Self, EllipticCurveError>;
}

pub trait IsPairing {
    type G1Point: IsGroup;
    type G2Point: IsGroup;
    type OutputField: IsField;
    type BaseField: IsField;

    /// Compute the product of the pairings for a list of point pairs.
    fn compute_batch(
        pairs: &[(&Self::G1Point, &Self::G2Point)],
    ) -> Result<FieldElement<Self::OutputField>, PairingError>;

    /// Compute the ate pairing between point `p` in G1 and `q` in G2.
    fn compute(
        p: &Self::G1Point,
        q: &Self::G2Point,
    ) -> Result<FieldElement<Self::OutputField>, PairingError> {
        Self::compute_batch(&[(p, q)])
    }
}
