use lambdaworks_math::field::{element::FieldElement, traits::IsField};
use std::borrow::Borrow;

use crate::fiat_shamir::transcript::Transcript;

// For Non-Hiding
// For batching operations we use a transcript to supply random values. In the case of kzg
// - Using an option for the transcript was the simplest way to enforce domain separation (prover/verifier)
//  for the future I think each protocol should have its own domain separated transcript within its instance variables
pub trait IsPolynomialCommitmentScheme<F: IsField> {
    /// Allows for Univariate vs Multilinear PCS
    type Polynomial;
    /// Point the polynomial is evaluated at
    type Point;
    /// Commitment to a Polynomial
    type Commitment;
    /// Allows for different proof structures
    type Proof;

    fn commit(&self, p: &Self::Polynomial) -> Self::Commitment;

    fn open(
        &self,
        point: impl Borrow<Self::Point>,
        eval: &FieldElement<F>,
        poly: &Self::Polynomial,
        transcript: Option<&mut dyn Transcript>,
    ) -> Self::Proof;

    fn open_batch(
        &self,
        point: impl Borrow<Self::Point>,
        eval: &[FieldElement<F>],
        polys: &[Self::Polynomial],
        transcript: Option<&mut dyn Transcript>,
    ) -> Self::Proof;

    fn verify(
        &self,
        point: impl Borrow<Self::Point>,
        eval: &FieldElement<F>,
        p_commitment: &Self::Commitment,
        proof: &Self::Proof,
        transcript: Option<&mut dyn Transcript>,
    ) -> bool;

    fn verify_batch(
        &self,
        point: impl Borrow<Self::Point>,
        evals: &[FieldElement<F>],
        p_commitments: &[Self::Commitment],
        proof: &Self::Proof,
        transcript: Option<&mut dyn Transcript>,
    ) -> bool;
}
