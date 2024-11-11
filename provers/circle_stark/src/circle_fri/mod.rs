pub mod circle_fri_commitment;
pub mod circle_fri_functions;

use crate::config::{BatchedMerkleTree, BatchedMerkleTreeBackend};
use circle_fri_commitment::CircleFriLayer;
use circle_fri_functions::{fold_circle_to_univariate, fold_univariate};
use lambdaworks_crypto::fiat_shamir::is_transcript::IsTranscript;
use lambdaworks_math::{
    circle::{cosets::Coset, point::CirclePoint, polynomial::evaluate_cfft},
    field::{
        element::FieldElement,
        fields::mersenne31::field::Mersenne31Field,
        traits::{IsField, IsSubFieldOf},
    },
    polynomial::Polynomial,
    traits::AsBytes,
};

type F = Mersenne31Field;
type FE = FieldElement<F>;

pub fn commit_phase(
    number_layers: usize,
    p_0: &Polynomial<FieldElement<F>>,
    transcript: &mut impl IsTranscript<F>,
    initial_domain: &Coset,
) -> (FieldElement<F>, Vec<CircleFriLayer>) {
    let mut domain_size = 1 << initial_domain.log_2_size;
    let mut fri_layer_list = Vec::with_capacity(number_layers);
    let mut current_evaluations = evaluate_cfft(p_0.coefficients().to_vec());

    let mut coset_shift = initial_domain.shift.clone();
    let mut domain = initial_domain.clone();

    // Create and commit the initial layer
    let current_layer = new_circle_fri_layer(&current_evaluations, &coset_shift, domain_size);
    let root = &current_layer.merkle_tree.root;
    transcript.append_bytes(root.as_ref());
    fri_layer_list.push(current_layer);

    // Process intermediate layers
    for layer_index in 0..number_layers - 1 {
        // Receive challenge beta
        let beta = transcript.sample_field_element();

        // Update domain parameters
        domain_size /= 2;
        coset_shift = CirclePoint::double(&coset_shift);
        domain = Coset::new(domain.log_2_size - 1, coset_shift.clone());

        // Fold evaluations
        current_evaluations = if layer_index == 0 {
            fold_circle_to_univariate(&current_evaluations, &beta, &domain)
        } else {
            fold_univariate(&current_evaluations, &beta, &domain, layer_index)
        };

        // Create and commit new layer
        let current_layer = new_circle_fri_layer(&current_evaluations, &coset_shift, domain_size);
        let root = &current_layer.merkle_tree.root;
        transcript.append_bytes(root.as_ref());
        fri_layer_list.push(current_layer);
    }

    // Process final layer
    let beta = transcript.sample_field_element();

    let final_evaluations = if number_layers == 1 {
        fold_circle_to_univariate(&current_evaluations, &beta, &domain)
    } else {
        fold_univariate(&current_evaluations, &beta, &domain, number_layers - 1)
    };

    let final_value = final_evaluations[0].clone();
    // Append final value to transcript
    transcript.append_field_element(&final_value);

    (final_value, fri_layer_list)
}

fn new_circle_fri_layer(
    evaluations: &[FE],
    coset_offset: &CirclePoint<F>,
    domain_size: usize,
) -> CircleFriLayer {
    // Build Merkle tree
    let merkle_tree = BatchedMerkleTree::build(&[evaluations.to_vec()]).unwrap();

    // Create new layer
    CircleFriLayer::new(
        evaluations,
        merkle_tree,
        coset_offset.x.clone(),
        domain_size,
    )
}
