#[cfg(test)]
mod integration_tests;

use ark_ff::PrimeField;
use ark_relations::r1cs::{ConstraintSystemRef, Field};
use lambdaworks_groth16_bls12_381::{common::*, r1cs::R1CS, ConstraintSystem};
use lambdaworks_math::traits::ByteConversion;

use std::ops::Deref;

/// Accepts an Arkworks circuit as a ConstraintSystem reference, and creates a
/// Lambdaworks ConstraintSystem, which can then be used with the Lambdaworks
/// Groth16 backend for setup, proving, and verification.
pub fn arkworks_cs_to_lambda_cs<ArkF: PrimeField>(
    cs: &ConstraintSystemRef<ArkF>,
) -> ConstraintSystem<FrField> {
    ConstraintSystem {
        constraints: r1cs_from_arkworks_cs(cs),
        witness: extract_witness_from_arkworks_cs(cs),
    }
}

#[inline]
fn r1cs_from_arkworks_cs<F: PrimeField>(cs: &ConstraintSystemRef<F>) -> R1CS {
    cs.inline_all_lcs();

    let r1cs_matrices = cs.to_matrices().unwrap();
    let num_pub_vars = cs.num_instance_variables();
    let total_variables = cs.num_witness_variables() + num_pub_vars;

    R1CS::from_matrices(
        ark_to_lambda_matrix(&r1cs_matrices.a, total_variables),
        ark_to_lambda_matrix(&r1cs_matrices.b, total_variables),
        ark_to_lambda_matrix(&r1cs_matrices.c, total_variables),
        num_pub_vars,
    )
}

#[inline]
fn ark_to_lambda_matrix<F: Field>(
    m: &[Vec<(F, usize)>],
    total_variables: usize,
) -> Vec<Vec<FrElement>> {
    sparse_matrix_to_dense(&arkworks_matrix_fps_to_fr_elements(m), total_variables)
}

#[inline]
fn arkworks_matrix_fps_to_fr_elements<F: Field>(
    m: &[Vec<(F, usize)>],
) -> Vec<Vec<(FrElement, usize)>> {
    m.iter()
        .map(|x| {
            x.iter()
                .map(|(x, y)| (ark_fr_to_fr_element(x), *y))
                .collect()
        })
        .collect()
}

#[inline]
fn extract_witness_from_arkworks_cs<F: PrimeField>(cs: &ConstraintSystemRef<F>) -> Vec<FrElement> {
    let binding = cs.borrow().unwrap();
    let borrowed_cs_ref = binding.deref();

    // Place public variables first, then witness assignments.
    // That's how Lambdaworks Groth16 expects the witness vector.
    let mut witness = vec![];
    witness.extend(
        borrowed_cs_ref
            .instance_assignment
            .iter()
            .map(ark_fr_to_fr_element),
    );
    witness.extend(
        borrowed_cs_ref
            .witness_assignment
            .iter()
            .map(ark_fr_to_fr_element),
    );
    witness
}

#[inline]
fn ark_fr_to_fr_element<F: Field>(ark_fq: &F) -> FrElement {
    let mut buff = Vec::<u8>::new();
    ark_fq.serialize_compressed(&mut buff).unwrap();
    FrElement::from_bytes_le(&buff).unwrap()
}

#[inline]
fn sparse_matrix_to_dense(
    m: &[Vec<(FrElement, usize)>],
    total_variables: usize,
) -> Vec<Vec<FrElement>> {
    m.iter()
        .map(|row| sparse_row_to_dense(row, total_variables))
        .collect()
}

#[inline]
fn sparse_row_to_dense(row: &[(FrElement, usize)], total_variables: usize) -> Vec<FrElement> {
    let mut dense_row = vec![FrElement::from(0); total_variables];
    row.iter().for_each(|e| {
        dense_row[e.1] = e.0.clone();
    });
    dense_row
}
