pub mod fri_commitment;
pub mod fri_decommit;
mod fri_functions;

use lambdaworks_crypto::fiat_shamir::is_transcript::IsTranscript;
use lambdaworks_math::field::traits::{IsFFTField, IsField};
use lambdaworks_math::traits::AsBytes;
use lambdaworks_math::{
    fft::cpu::bit_reversing::in_place_bit_reverse_permute, field::traits::IsSubFieldOf,
};
pub use lambdaworks_math::{
    field::{element::FieldElement, fields::u64_prime_field::U64PrimeField},
    polynomial::Polynomial,
};

use crate::config::{BatchedMerkleTree, BatchedMerkleTreeBackend};

use self::fri_commitment::FriLayer;
use self::fri_decommit::FriDecommitment;
use self::fri_functions::fold_polynomial;

pub fn commit_phase<F: IsFFTField + IsSubFieldOf<E>, E: IsField>(
    number_layers: usize,
    p_0: Polynomial<FieldElement<E>>,
    transcript: &mut impl IsTranscript<E>,
    coset_offset: &FieldElement<F>,
    domain_size: usize,
) -> (
    FieldElement<E>,
    Vec<FriLayer<E, BatchedMerkleTreeBackend<E>>>,
)
where
    FieldElement<F>: AsBytes + Sync + Send,
    FieldElement<E>: AsBytes + Sync + Send,
{
    let mut domain_size = domain_size;

    let mut fri_layer_list = Vec::with_capacity(number_layers);
    let mut current_layer: FriLayer<E, BatchedMerkleTreeBackend<E>>;
    let mut current_poly = p_0;

    let mut coset_offset = coset_offset.clone();

    for _ in 1..number_layers {
        // <<<< Receive challenge 𝜁ₖ₋₁
        let zeta = transcript.sample_field_element();
        coset_offset = coset_offset.square();
        domain_size /= 2;

        // Compute layer polynomial and domain
        current_poly = FieldElement::<F>::from(2) * fold_polynomial(&current_poly, &zeta);
        current_layer = new_fri_layer(&current_poly, &coset_offset, domain_size);
        let new_data = &current_layer.merkle_tree.root;
        fri_layer_list.push(current_layer.clone()); // TODO: remove this clone

        // >>>> Send commitment: [pₖ]
        transcript.append_bytes(new_data);
    }

    // <<<< Receive challenge: 𝜁ₙ₋₁
    let zeta = transcript.sample_field_element();

    let last_poly = FieldElement::<F>::from(2) * fold_polynomial(&current_poly, &zeta);

    let last_value = last_poly
        .coefficients()
        .first()
        .unwrap_or(&FieldElement::zero())
        .clone();

    // >>>> Send value: pₙ
    transcript.append_field_element(&last_value);

    (last_value, fri_layer_list)
}

pub fn query_phase<F: IsField>(
    fri_layers: &Vec<FriLayer<F, BatchedMerkleTreeBackend<F>>>,
    iotas: &[usize],
) -> Vec<FriDecommitment<F>>
where
    FieldElement<F>: AsBytes + Sync + Send,
{
    if !fri_layers.is_empty() {
        let query_list = iotas
            .iter()
            .map(|iota_s| {
                let mut layers_evaluations_sym = Vec::new();
                let mut layers_auth_paths_sym = Vec::new();

                let mut index = *iota_s;
                for layer in fri_layers {
                    // symmetric element
                    let evaluation_sym = layer.evaluation[index ^ 1].clone();
                    let auth_path_sym = layer.merkle_tree.get_proof_by_pos(index >> 1).unwrap();
                    layers_evaluations_sym.push(evaluation_sym);
                    layers_auth_paths_sym.push(auth_path_sym);

                    index >>= 1;
                }

                FriDecommitment {
                    layers_auth_paths: layers_auth_paths_sym,
                    layers_evaluations_sym,
                }
            })
            .collect();

        query_list
    } else {
        vec![]
    }
}

pub fn new_fri_layer<F: IsFFTField + IsSubFieldOf<E>, E: IsField>(
    poly: &Polynomial<FieldElement<E>>,
    coset_offset: &FieldElement<F>,
    domain_size: usize,
) -> crate::fri::fri_commitment::FriLayer<E, BatchedMerkleTreeBackend<E>>
where
    FieldElement<F>: AsBytes + Sync + Send,
    FieldElement<E>: AsBytes + Sync + Send,
{
    let mut evaluation =
        Polynomial::evaluate_offset_fft(poly, 1, Some(domain_size), coset_offset).unwrap(); // TODO: return error

    in_place_bit_reverse_permute(&mut evaluation);

    let mut to_commit = Vec::new();
    for chunk in evaluation.chunks(2) {
        to_commit.push(vec![chunk[0].clone(), chunk[1].clone()]);
    }

    let merkle_tree = BatchedMerkleTree::build(&to_commit).unwrap();

    FriLayer::new(
        &evaluation,
        merkle_tree,
        coset_offset.clone().to_extension(),
        domain_size,
    )
}
