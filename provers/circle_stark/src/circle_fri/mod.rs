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

pub fn commit_phase<F:isFFTField + IsSubFieldOf<E>,E:IsField(

    number_layers: usize,
    p_0: Polynomial<FieldElement<E>>,
    transcript: &mut impl IsTranscript<E>,
    coset_offset: &FieldElement<F>,
    domain_size: usize,
) -> (
    FieldElement<E>,
    Vec<FriLayer<E, BatchedMerkleTreeBackend<E>>>,
    
)

pub fn commit_phase<E: IsField>(
    number_layers: usize,
    p_0: &Polynomial<FieldElement<F>>,
    transcript: &mut impl IsTranscript<E>,
    initial_domain: &Coset,
) -> (
    FieldElement<F>,
    Vec<CircleFriLayer<F, BatchedMerkleTreeBackend<F>>>,
)
where
    FieldElement<F>: AsBytes + Sync + Send,
    FieldElement<E>: AsBytes + Sync + Send,
    F: IsSubFieldOf<E>,
    E: IsSubFieldOf<F>,
{
    let mut domain_size = 1 << initial_domain.log_2_size;
    let mut fri_layer_list = Vec::with_capacity(number_layers);
    let mut current_layer: CircleFriLayer<F, BatchedMerkleTreeBackend<F>>;
    let mut current_evaluations = evaluate_cfft(p_0.coefficients().to_vec());

    let mut coset_offset = initial_domain.shift.clone();
    let mut domain = initial_domain.clone();

    // Create and commit the initial layer
    current_layer = new_circle_fri_layer(&current_evaluations, &coset_offset, domain_size);
    let root = &current_layer.merkle_tree.root;
    transcript.append_bytes(root.as_ref());
    fri_layer_list.push(current_layer.clone());

    // Process intermediate layers
    for layer_index in 0..number_layers - 1 {
        // Receive challenge beta
        let beta_e = transcript.sample_field_element();
        let beta_f = FieldElement::<F>::from_element(&beta_e);

        // Update domain parameters
        domain_size /= 2;
        coset_offset = CirclePoint::double(&coset_offset);
        domain = Coset::new(domain.log_2_size - 1, coset_offset.clone());

        // Fold evaluations
        current_evaluations = if layer_index == 0 {
            fold_circle_to_univariate(&current_evaluations, &beta_f, &domain)
        } else {
            fold_univariate(&current_evaluations, &beta_f, &domain, layer_index)
        };

        // Create and commit new layer
        current_layer = new_circle_fri_layer(&current_evaluations, &coset_offset, domain_size);
        let root = &current_layer.merkle_tree.root;
        transcript.append_bytes(root.as_ref());
        fri_layer_list.push(current_layer.clone());
    }

    // Process final layer
    let beta_e = transcript.sample_field_element();
    let beta_f = FieldElement::<F>::from_element(&beta_e);

    let final_evaluations = if number_layers == 1 {
        fold_circle_to_univariate(&current_evaluations, &beta_f, &domain)
    } else {
        fold_univariate(&current_evaluations, &beta_f, &domain, number_layers - 1)
    };

    let final_value = final_evaluations[0].clone();
    // Convert final_value to E type for transcript
    let final_value_e = final_value.to_extension();
    transcript.append_field_element(&final_value_e);

    (final_value, fri_layer_list)
}

fn new_circle_fri_layer(
    evaluations: &[FieldElement<F>],
    coset_offset: &CirclePoint<F>,
    domain_size: usize,
) -> CircleFriLayer<F, BatchedMerkleTreeBackend<F>> {
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
