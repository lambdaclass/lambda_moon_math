use std::collections::{BTreeSet, HashMap, HashSet};

use lambdaworks_crypto::merkle_tree::proof::Proof;
use lambdaworks_math::{
    field::{
        element::FieldElement,
        fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
        traits::{IsField, IsSubFieldOf},
    },
    traits::AsBytes,
};

use crate::{
    config::Commitment,
    domain::Domain,
    fri::fri_decommit::FriDecommitment,
    table::Table,
    traits::AIR,
    transcript::StoneProverTranscript,
    verifier::{IsStarkVerifier, Verifier},
};

use super::options::ProofOptions;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PolynomialOpenings<F: IsField> {
    pub proof: Proof<Commitment>,
    pub proof_sym: Proof<Commitment>,
    pub evaluations: Vec<FieldElement<F>>,
    pub evaluations_sym: Vec<FieldElement<F>>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct DeepPolynomialOpening<F: IsSubFieldOf<E>, E: IsField> {
    pub composition_poly: PolynomialOpenings<E>,
    pub main_trace_polys: PolynomialOpenings<F>,
    pub aux_trace_polys: Option<PolynomialOpenings<E>>,
}

pub type DeepPolynomialOpenings<F, E> = Vec<DeepPolynomialOpening<F, E>>;

#[derive(Debug, serde::Serialize, serde::Deserialize)]
pub struct StarkProof<F: IsSubFieldOf<E>, E: IsField> {
    // Length of the execution trace
    pub trace_length: usize,
    // Commitments of the trace columns
    // [tⱼ]
    pub lde_trace_main_merkle_root: Commitment,
    // Commitments of auxiliary trace columns
    // [tⱼ]
    pub lde_trace_aux_merkle_root: Option<Commitment>,
    // tⱼ(zgᵏ)
    pub trace_ood_evaluations: Table<E>,
    // Commitments to Hᵢ
    pub composition_poly_root: Commitment,
    // Hᵢ(z^N)
    pub composition_poly_parts_ood_evaluation: Vec<FieldElement<E>>,
    // [pₖ]
    pub fri_layers_merkle_roots: Vec<Commitment>,
    // pₙ
    pub fri_last_value: FieldElement<E>,
    // Open(pₖ(Dₖ), −𝜐ₛ^(2ᵏ))
    pub query_list: Vec<FriDecommitment<E>>,
    // Open(H₁(D_LDE, 𝜐ᵢ), Open(H₂(D_LDE, 𝜐ᵢ), Open(tⱼ(D_LDE), 𝜐ᵢ)
    // Open(H₁(D_LDE, -𝜐ᵢ), Open(H₂(D_LDE, -𝜐ᵢ), Open(tⱼ(D_LDE), -𝜐ᵢ)
    pub deep_poly_openings: DeepPolynomialOpenings<F, E>,
    // nonce obtained from grinding
    pub nonce: Option<u64>,
}

/// Serializer compatible with Stone prover
/// (https://github.com/starkware-libs/stone-prover/)
pub struct StoneCompatibleSerializer;

impl StoneCompatibleSerializer {
    pub fn serialize_proof<A>(
        proof: &StarkProof<Stark252PrimeField, Stark252PrimeField>,
        public_inputs: &A::PublicInputs,
        options: &ProofOptions,
    ) -> Vec<u8>
    where
        A: AIR<Field = Stark252PrimeField, FieldExtension = Stark252PrimeField>,
        A::PublicInputs: AsBytes,
    {
        let mut output = Vec::new();

        Self::append_trace_commitment(proof, &mut output);
        Self::append_composition_polynomial_commitment(proof, &mut output);
        Self::append_out_of_domain_evaluations(proof, &mut output);
        Self::append_fri_commit_phase_commitments(proof, &mut output);
        Self::append_proof_of_work_nonce(proof, &mut output);

        let fri_query_indexes = Self::get_fri_query_indexes::<A>(proof, public_inputs, options);
        Self::append_fri_query_phase_first_layer(proof, &fri_query_indexes, &mut output);
        Self::append_fri_query_phase_inner_layers(proof, &fri_query_indexes, &mut output);

        output
    }

    /// Appends the root bytes of the Merkle tree for the main trace, and if there is a RAP round,
    /// it also appends the root bytes of the Merkle tree for the extended columns.
    fn append_trace_commitment(
        proof: &StarkProof<Stark252PrimeField, Stark252PrimeField>,
        output: &mut Vec<u8>,
    ) {
        output.extend_from_slice(&proof.lde_trace_main_merkle_root);

        if let Some(lde_trace_aux_merkle_root) = proof.lde_trace_aux_merkle_root {
            output.extend_from_slice(&lde_trace_aux_merkle_root);
        }
    }

    /// Appends the root bytes of the Merkle tree for the composition polynomial.
    fn append_composition_polynomial_commitment(
        proof: &StarkProof<Stark252PrimeField, Stark252PrimeField>,
        output: &mut Vec<u8>,
    ) {
        output.extend_from_slice(&proof.composition_poly_root);
    }

    /// Appends the bytes of the evaluations of the trace `t_1, ..., t_m` and composition polynomial parts
    /// `H_1, ..., H_s` at the out of domain challenge `z`, its shifts `g^i z` and its power `z^s`, respectively.
    /// These are sorted as follows: first the evaluations of the trace in increasing order of
    /// trace column and shift number. Then all the evaluations of the parts of the composition
    /// polynomial. That is:
    ///
    /// t_1(z), ..., t_1(g^K z), t_2(z), ..., t_2(g^K z), ..., t_m(g z), ..., t_m(g^K z), H_1(z^s), ..., H_s(z^s).
    ///
    /// Here, K is the length of the frame size.
    fn append_out_of_domain_evaluations(
        proof: &StarkProof<Stark252PrimeField, Stark252PrimeField>,
        output: &mut Vec<u8>,
    ) {
        for i in 0..proof.trace_ood_evaluations.width {
            for j in 0..proof.trace_ood_evaluations.height {
                output.extend_from_slice(&proof.trace_ood_evaluations.get_row(j)[i].as_bytes());
            }
        }

        for elem in proof.composition_poly_parts_ood_evaluation.iter() {
            output.extend_from_slice(&elem.as_bytes());
        }
    }

    /// Appends the commitments to the inner layers of FRI followed by the element of the last layer.
    fn append_fri_commit_phase_commitments(
        proof: &StarkProof<Stark252PrimeField, Stark252PrimeField>,
        output: &mut Vec<u8>,
    ) {
        output.extend_from_slice(
            &proof
                .fri_layers_merkle_roots
                .iter()
                .flatten()
                .cloned()
                .collect::<Vec<_>>(),
        );

        output.extend_from_slice(&proof.fri_last_value.as_bytes());
    }

    /// Appends the proof of work nonce in case there is one. There could be none if the `grinding_factor`
    /// was set to 0 during proof generation. In that case nothing is appended.
    fn append_proof_of_work_nonce(
        proof: &StarkProof<Stark252PrimeField, Stark252PrimeField>,
        output: &mut Vec<u8>,
    ) {
        if let Some(nonce_value) = proof.nonce {
            output.extend_from_slice(&nonce_value.to_be_bytes());
        }
    }

    /// Appends the values and authentication paths of the trace and composition polynomial parts
    /// needed for the first layer of FRI. Next we describe the order in which these are appended.
    ///
    /// Each FRI query index `i` determines a pair of elements `d_i` and `-d_i` on the domain of the
    /// first layer.
    /// Let BT_i be the concatenation of the bytes of the following values
    /// t_1(d_i), t_2(d_i), ..., t_m(d_i), t_1(-d_i), t_2(-d_i), ..., t_m(-d_i),
    /// where m is the total number of columns, including RAP extended ones.
    /// Similarly, let BH_i be the concatenation of the bytes of the following elements
    /// H_1(d_i), ..., H_s(d_i), H_1(-d_i), ..., H_s(-d_i),
    /// where s is the number of parts into which the composition polynomial was broken.
    ///
    /// If i_1, ..., i_k are all the FRI query indexes sorted in increasing order and without repeated
    /// values, then this method appends the following to the output:
    ///
    /// BT_{i_1} | BT_{i_2} | ... | BT_{i_k} | TraceMergedPaths | BH_{i_1} | BH_{i_2} | ... | B_{i_k} | CompositionMergedPaths.
    ///
    /// Where TraceMergedPaths is the merged authentication paths of the trace Merkle tree for all queries
    /// and similarly, CompositionMergedPaths is the merged authentication paths of the composition polynomial
    /// Merkle tree for all queries (see the `merge_authentication_paths` method).
    ///
    /// Example:
    /// If there are 6 queries [3, 1, 5, 2, 1, 3], then this method appends the
    /// following to the output:
    /// `BT_1 | BT_2 | BT_3 | BT_5 | TraceMergedPaths | BH_1 | BH_2 | BH_3 | BH_5 | CompositionMergedPaths`
    fn append_fri_query_phase_first_layer(
        proof: &StarkProof<Stark252PrimeField, Stark252PrimeField>,
        fri_query_indexes: &[usize],
        output: &mut Vec<u8>,
    ) {
        let mut fri_first_layer_openings: Vec<_> = proof
            .deep_poly_openings
            .iter()
            .zip(fri_query_indexes.iter())
            .collect();
        // Remove repeated values
        let mut seen = HashSet::new();
        fri_first_layer_openings.retain(|&(_, index)| seen.insert(index));
        // Sort by increasing value of query
        fri_first_layer_openings.sort_by(|a, b| a.1.cmp(b.1));

        // Append BT_{i_1} | BT_{i_2} | ... | BT_{i_k}
        for (opening, _) in fri_first_layer_openings.iter() {
            for elem in opening.main_trace_polys.evaluations.iter() {
                output.extend_from_slice(&elem.as_bytes());
            }
            if let Some(aux) = &opening.aux_trace_polys {
                for elem in aux.evaluations.iter() {
                    output.extend_from_slice(&elem.as_bytes());
                }
            }

            for elem in opening.main_trace_polys.evaluations_sym.iter() {
                output.extend_from_slice(&elem.as_bytes());
            }
            if let Some(aux) = &opening.aux_trace_polys {
                for elem in aux.evaluations_sym.iter() {
                    output.extend_from_slice(&elem.as_bytes());
                }
            }
        }

        let fri_trace_query_indexes: Vec<_> = fri_query_indexes
            .iter()
            .flat_map(|query| vec![query * 2, query * 2 + 1])
            .collect();

        // Append TraceMergedPaths
        //    Main trace
        let fri_trace_paths: Vec<_> = proof
            .deep_poly_openings
            .iter()
            .flat_map(|opening| {
                vec![
                    &opening.main_trace_polys.proof,
                    &opening.main_trace_polys.proof_sym,
                ]
            })
            .collect();
        let nodes = Self::merge_authentication_paths(&fri_trace_paths, &fri_trace_query_indexes);
        for node in nodes.iter() {
            output.extend_from_slice(node);
        }

        //    Aux trace
        let mut all_openings_aux_trace_polys_are_some = true;
        let mut fri_trace_paths: Vec<&Proof<Commitment>> = Vec::new();
        for opening in proof.deep_poly_openings.iter() {
            if let Some(aux_trace_polys) = &opening.aux_trace_polys {
                fri_trace_paths.push(&aux_trace_polys.proof);
                fri_trace_paths.push(&aux_trace_polys.proof_sym);
            } else {
                all_openings_aux_trace_polys_are_some = false;
            }
        }
        if all_openings_aux_trace_polys_are_some {
            let nodes =
                Self::merge_authentication_paths(&fri_trace_paths, &fri_trace_query_indexes);
            for node in nodes.iter() {
                output.extend_from_slice(node);
            }
        }

        // Append BH_{i_1} | BH_{i_2} | ... | B_{i_k}
        for (opening, _) in fri_first_layer_openings.iter() {
            for elem in opening.composition_poly.evaluations.iter() {
                output.extend_from_slice(&elem.as_bytes());
            }
            for elem in opening.composition_poly.evaluations_sym.iter() {
                output.extend_from_slice(&elem.as_bytes());
            }
        }

        // Append CompositionMergedPaths
        let fri_composition_paths: Vec<_> = proof
            .deep_poly_openings
            .iter()
            .map(|opening| &opening.composition_poly.proof)
            .collect();
        let nodes = Self::merge_authentication_paths(&fri_composition_paths, fri_query_indexes);
        for node in nodes.iter() {
            output.extend_from_slice(node);
        }
    }

    /// Appends the values and authentication paths needed for the inner layers of FRI.
    /// Just as in the append_fri_query_phase_first_layer, for each layer, the authentication
    /// paths are merged and the redundant field elements are not sent, in order to optimize
    /// the size of the proof. When having multiple queries we can have repeated field elements
    /// for two reasons: either we are sending two times the same field element because of
    /// a repeated query, or we are sending a field element that the verifier could simply
    /// derive from values from previous layers.
    ///
    /// For each layer i there are:
    /// - X_i = { p_i(-d_j), p_i(d_j) for all queries j }, the elements the verifier needs.
    /// - Y_i = { p_i( d_j) for all queries j }, the elements that the verifier computes from
    ///         previous layers.
    /// - Z_i = X_i - Y_i, the elements that the verifier needs but cannot compute from previous layers.
    ///         sorted by increasing value of query.  
    /// - MergedPathsLayer_i: the merged authentication paths for all p_i(-d_j) and p_i(d_j).
    ///
    /// This method appends:
    ///
    /// Z_1 | MergedPathsLayer_1 | Z_2 | MergedPathsLayer_2 | ... | Z_n | MergedPathsLayer_n,
    ///
    /// where n is the total number of FRI layers.
    fn append_fri_query_phase_inner_layers(
        proof: &StarkProof<Stark252PrimeField, Stark252PrimeField>,
        fri_query_indexes: &[usize],
        output: &mut Vec<u8>,
    ) {
        let mut fri_layers_evaluations: HashMap<(u64, usize, usize), FieldElement<_>> =
            HashMap::new();
        for (decommitment, query_index) in proof.query_list.iter().zip(fri_query_indexes.iter()) {
            let mut query_layer_index = *query_index;
            for (i, element) in decommitment.layers_evaluations_sym.iter().enumerate() {
                fri_layers_evaluations.insert(
                    (
                        i as u64,
                        query_layer_index >> 1,
                        (query_layer_index + 1) % 2,
                    ),
                    *element,
                );
                query_layer_index >>= 1;
            }
        }

        let mut indexes_previous_layer = fri_query_indexes.to_owned();
        for i in 0..proof.query_list[0].layers_evaluations_sym.len() {
            // Compute set Y_i
            let reconstructed_row_col: BTreeSet<_> = indexes_previous_layer
                .iter()
                .map(|index| (index >> 1, index % 2))
                .collect();

            // Compute set X_i
            let reconstructed_row_col_sym: BTreeSet<_> = reconstructed_row_col
                .iter()
                .map(|(x, y)| (*x, 1 - y))
                .collect();

            // Compute set Z_i
            let row_col_to_send: Vec<_> = reconstructed_row_col_sym
                .difference(&reconstructed_row_col)
                .collect();

            // Append Z_i
            for element in row_col_to_send
                .iter()
                .map(|(row, col)| &fri_layers_evaluations[&(i as u64, *row, *col)])
            {
                output.extend_from_slice(&element.as_bytes());
            }

            indexes_previous_layer = indexes_previous_layer
                .iter()
                .map(|index| index >> 1)
                .collect();

            let layer_auth_paths: Vec<_> = proof
                .query_list
                .iter()
                .map(|decommitment| &decommitment.layers_auth_paths[i])
                .collect();

            // Append MergedPathsLayer_i
            let nodes =
                Self::merge_authentication_paths(&layer_auth_paths, &indexes_previous_layer);
            for node in nodes.iter() {
                output.extend_from_slice(node);
            }
        }
    }

    /// Merges `n` authentication paths for `n` leaves into a list of the minimal number of nodes
    /// needed to reach the Merkle root for all of them. The nodes of the merged authentication
    /// paths are sorted from level 0 to the hightest level of the Merkle tree, and nodes at the
    /// same level are sorted from left to right.
    ///
    /// Let's consider some examples. Suppose the Merkle tree is as follows:
    ///
    /// Root              ABCD
    ///                  /    \
    /// Level 1         AB    CD
    ///                /  \  /  \
    /// Level 0       A   B C    D
    /// Leaf index    0   1 2    3
    ///
    /// All authentication paths are `pA = [B, CD]`, `pB = [A, CD]`, `pC = [D, AB]` and `pD = [C, AB]`.
    ///
    /// (1) `merge_authentication_paths([pA, pB], [0, 1]) = [CD]`.
    /// (2) `merge_authentication_paths([pC, pA], [2, 0]) = [B, D]`.
    /// (3) `merge_authentication_paths([pA, pD], [0, 3]) = [B, C]`.
    /// (4) `merge_authentication_paths([pA, pD, pB], [0, 3, 1]) = [C]`.
    /// (5) `merge_authentication_paths([pA, pB, pC, pD], [0, 1, 2, 3]) = []`.
    ///
    /// Input:
    /// `authentication_paths`: The authentication paths to be merged.
    /// `leaf_indexes`: The leaf indexes corresponding to the authentication paths.
    ///
    /// Output:
    /// The merged authentication paths
    fn merge_authentication_paths(
        authentication_paths: &[&Proof<Commitment>],
        leaf_indexes: &[usize],
    ) -> Vec<Commitment> {
        debug_assert_eq!(leaf_indexes.len(), authentication_paths.len());
        let mut merkle_tree: HashMap<(usize, usize), Commitment> = HashMap::new();
        for (index_previous_layer, path) in leaf_indexes.iter().zip(authentication_paths.iter()) {
            let mut node_index = *index_previous_layer;
            for (tree_level, node) in path.merkle_path.iter().enumerate() {
                merkle_tree.insert((tree_level, node_index ^ 1), *node);
                node_index >>= 1;
            }
        }

        let mut result = Vec::new();
        let mut level_indexes: BTreeSet<usize> = leaf_indexes.iter().copied().collect();
        let merkle_tree_height = authentication_paths[0].merkle_path.len();
        for tree_level in 0..merkle_tree_height {
            for node_index in level_indexes.iter() {
                let sibling_index = node_index ^ 1;
                if !level_indexes.contains(&sibling_index) {
                    let node = &merkle_tree[&(tree_level, sibling_index)];
                    result.push(*node);
                }
            }
            level_indexes = level_indexes.iter().map(|index| *index >> 1).collect();
        }
        result
    }
    fn get_fri_query_indexes<A>(
        proof: &StarkProof<Stark252PrimeField, Stark252PrimeField>,
        public_inputs: &A::PublicInputs,
        proof_options: &ProofOptions,
    ) -> Vec<usize>
    where
        A: AIR<Field = Stark252PrimeField, FieldExtension = Stark252PrimeField>,
        A::PublicInputs: AsBytes,
    {
        let mut transcript = StoneProverTranscript::new(&public_inputs.as_bytes());
        let air = A::new(proof.trace_length, public_inputs, proof_options);
        let domain = Domain::<Stark252PrimeField>::new(&air);
        let challenges = Verifier::step_1_replay_rounds_and_recover_challenges(
            &air,
            proof,
            &domain,
            &mut transcript,
        );
        challenges.iotas
    }
}

// #[cfg(test)]
// mod tests {
//     use lambdaworks_math::{field::element::FieldElement, traits::AsBytes};

//     use crate::{
//         examples::fibonacci_2_cols_shifted::{self, Fibonacci2ColsShifted},
//         proof::{options::ProofOptions, stark::StoneCompatibleSerializer},
//         prover::{IsStarkProver, Prover},
//         transcript::StoneProverTranscript,
//     };

//     #[test]
//     fn test_serialization_compatible_with_stone_1() {
//         let mut trace = fibonacci_2_cols_shifted::compute_trace(FieldElement::one(), 4);

//         let claimed_index = 3;
//         let claimed_value = trace.get_row(claimed_index)[0];
//         let proof_options = ProofOptions {
//             blowup_factor: 4,
//             coset_offset: 3,
//             grinding_factor: 0,
//             fri_number_of_queries: 1,
//         };

//         let pub_inputs = fibonacci_2_cols_shifted::PublicInputs {
//             claimed_value,
//             claimed_index,
//         };

//         let proof = Prover::<Fibonacci2ColsShifted<_>>::prove(
//             &mut trace,
//             &pub_inputs,
//             &proof_options,
//             StoneProverTranscript::new(&pub_inputs.as_bytes()),
//         )
//         .unwrap();

//         let expected_bytes = [
//             14, 185, 220, 192, 251, 24, 84, 87, 42, 1, 35, 103, 83, 206, 5, 19, 157, 57, 42, 163,
//             174, 175, 231, 42, 191, 241, 80, 254, 33, 23, 85, 148, 240, 15, 97, 4, 33, 250, 73, 57,
//             153, 20, 91, 112, 71, 103, 155, 245, 134, 85, 150, 224, 103, 5, 176, 183, 152, 52, 190,
//             56, 94, 184, 211, 203, 1, 10, 170, 210, 58, 15, 137, 139, 84, 215, 101, 26, 236, 253,
//             138, 16, 34, 94, 85, 246, 117, 36, 122, 25, 65, 56, 39, 64, 182, 60, 92, 149, 0, 61,
//             238, 22, 79, 89, 52, 136, 161, 125, 245, 232, 111, 27, 91, 235, 0, 112, 73, 52, 122,
//             171, 178, 11, 249, 92, 149, 195, 95, 127, 77, 90, 0, 63, 48, 208, 46, 245, 248, 39,
//             173, 179, 161, 21, 59, 173, 210, 38, 117, 61, 159, 103, 129, 41, 200, 180, 127, 152,
//             136, 37, 52, 131, 168, 143, 7, 71, 166, 221, 33, 179, 43, 105, 109, 45, 26, 161, 194,
//             171, 13, 78, 139, 52, 158, 132, 170, 241, 155, 38, 213, 231, 199, 58, 181, 248, 101,
//             136, 5, 201, 25, 47, 164, 4, 56, 196, 188, 130, 134, 39, 128, 65, 210, 9, 124, 10, 82,
//             253, 146, 34, 57, 37, 92, 71, 2, 44, 3, 248, 124, 227, 206, 179, 238, 69, 200, 177, 31,
//             7, 171, 247, 48, 97, 185, 116, 237, 171, 117, 251, 207, 4, 66, 112, 144, 10, 255, 60,
//             207, 185, 25, 7, 110, 159, 3, 120, 156, 213, 179, 46, 1, 189, 58, 131, 21, 190, 194,
//             176, 219, 255, 172, 68, 21, 117, 44, 122, 177, 139, 62, 111, 251, 21, 15, 81, 246, 120,
//             6, 115, 221, 244, 77, 82, 191, 9, 150, 178, 205, 168, 196, 218, 39, 107, 231, 32, 56,
//             92, 78, 19, 151, 21, 18, 123, 158, 163, 7, 245, 10, 237, 4, 231, 187, 232, 154, 165,
//             106, 226, 45, 101, 155, 81, 137, 180, 78, 215, 206, 64, 112, 184, 156, 39, 46, 42, 36,
//             247, 61, 70, 15, 234, 20, 185, 1, 140, 34, 11, 178, 173, 48, 7, 105, 77, 50, 87, 59,
//             37, 216, 148, 24, 223, 199, 163, 177, 236, 104, 234, 237, 132, 97, 92, 248, 10, 244,
//             20, 3, 24, 68, 23, 101, 90, 108, 206, 210, 154, 100, 174, 118, 75, 177, 40, 49, 191,
//             143, 71, 99, 216, 209, 213, 219, 8, 194, 185, 240, 21, 232, 232, 145, 176, 192, 178,
//             75, 157, 0, 6, 123, 14, 250, 181, 8, 50, 183, 108, 249, 113, 146, 9, 22, 36, 212, 43,
//             134, 116, 6, 102, 197, 211, 105, 230, 153, 59, 4, 77, 178, 36, 68, 192, 192, 235, 241,
//             9, 91, 154, 81, 250, 235, 0, 28, 155, 77, 234, 54, 171, 233, 5, 247, 22, 38, 32, 219,
//             189, 80, 23, 171, 236, 163, 63, 168, 37, 118, 181, 197, 194, 198, 23, 146, 105, 59, 72,
//             201, 212, 65, 74, 64, 126, 239, 102, 182, 2, 157, 174, 7, 234, 6, 40, 218, 216, 248,
//             96, 149, 112, 209, 255, 173, 185, 81, 55, 105, 97, 6, 239, 20, 189, 183, 213, 184, 147,
//             226, 210, 8, 224, 248, 198, 11, 186, 7, 179, 148, 95, 226, 129, 234, 27, 46, 85, 20,
//             182, 118, 241, 2, 69, 184, 39, 231, 81, 9, 28, 60, 114, 120, 53, 252, 192, 115, 40,
//             213, 33, 160, 213, 41, 195, 61, 131, 42, 105, 77, 188, 109, 118, 53, 70, 24, 141, 94,
//             101, 222, 67, 254, 29, 157, 8, 184, 145, 194, 89, 189, 95, 253, 181, 90, 70, 207, 28,
//             53, 40, 246, 178, 129, 178, 83, 109, 24, 202, 136, 140, 211, 4, 167, 36, 253, 29, 66,
//             79, 250, 184, 63, 158, 162, 206, 83, 135, 251, 125, 215, 121, 149, 118, 82, 112, 53,
//             242, 58, 127, 196, 123, 63, 110, 192, 137, 125, 95, 72, 122, 82, 8, 121, 113, 241, 76,
//             255, 36, 96, 225, 5, 158, 234, 196, 65, 220, 82, 70, 244, 90, 62, 85, 201, 169, 9, 104,
//             29, 215, 76, 179, 28, 235, 123, 86, 98, 25, 254, 123, 7, 55, 26, 190, 66, 97, 62, 111,
//             150, 110, 27, 233, 81, 131, 192, 21, 72, 16, 130, 46, 43, 43, 213, 145, 96, 22, 224,
//             211, 5, 253, 114, 94, 164, 190, 87, 143, 69, 128, 1, 253, 103, 213, 139, 7, 35, 132,
//             52, 242, 55, 248, 72, 214, 102, 108, 57, 205, 46, 20, 1, 83, 198, 32, 167, 96, 242,
//             117, 87, 201,
//         ];

//         let serialized_proof = StoneCompatibleSerializer::serialize_proof::<Fibonacci2ColsShifted<_>>(
//             &proof,
//             &pub_inputs,
//             &proof_options,
//         );
//         assert_eq!(serialized_proof, expected_bytes);
//     }

//     #[test]
//     fn test_serialization_compatible_with_stone_case_2() {
//         let mut trace = fibonacci_2_cols_shifted::compute_trace(FieldElement::one(), 4);

//         let claimed_index = 2;
//         let claimed_value = trace.get_row(claimed_index)[0];
//         let proof_options = ProofOptions {
//             blowup_factor: 2,
//             coset_offset: 3,
//             grinding_factor: 0,
//             fri_number_of_queries: 10,
//         };

//         let pub_inputs = fibonacci_2_cols_shifted::PublicInputs {
//             claimed_value,
//             claimed_index,
//         };

//         let proof = Prover::<Fibonacci2ColsShifted<_>>::prove(
//             &mut trace,
//             &pub_inputs,
//             &proof_options,
//             StoneProverTranscript::new(&pub_inputs.as_bytes()),
//         )
//         .unwrap();
//         let expected_bytes = [
//             9, 161, 59, 243, 85, 60, 44, 155, 163, 203, 128, 147, 203, 253, 93, 16, 137, 42, 94,
//             225, 173, 254, 120, 1, 43, 167, 254, 15, 49, 148, 2, 50, 47, 159, 254, 83, 209, 26, 91,
//             113, 237, 157, 107, 252, 35, 130, 117, 22, 155, 181, 217, 11, 145, 208, 53, 201, 14,
//             148, 19, 247, 19, 105, 239, 108, 2, 17, 212, 103, 137, 76, 186, 20, 31, 118, 21, 139,
//             122, 39, 100, 197, 112, 11, 188, 227, 236, 15, 127, 186, 231, 187, 158, 137, 41, 180,
//             233, 36, 3, 145, 50, 37, 100, 198, 0, 152, 68, 30, 64, 111, 78, 211, 119, 49, 142, 180,
//             178, 74, 59, 150, 209, 45, 211, 159, 6, 216, 205, 161, 255, 142, 0, 104, 163, 169, 38,
//             160, 79, 95, 195, 96, 46, 20, 45, 189, 161, 181, 95, 133, 26, 124, 224, 44, 153, 119,
//             121, 29, 187, 126, 125, 161, 4, 45, 2, 1, 113, 106, 28, 174, 255, 138, 4, 34, 227, 191,
//             33, 203, 60, 20, 34, 36, 33, 8, 44, 53, 250, 177, 127, 59, 157, 229, 179, 87, 165, 58,
//             0, 0, 83, 98, 7, 41, 90, 187, 198, 80, 159, 250, 57, 252, 211, 64, 233, 110, 223, 155,
//             56, 189, 215, 57, 80, 161, 169, 246, 65, 133, 129, 132, 129, 233, 154, 204, 187, 178,
//             244, 76, 12, 9, 30, 113, 105, 206, 46, 192, 68, 96, 27, 72, 94, 126, 101, 253, 63, 94,
//             10, 89, 116, 120, 31, 123, 5, 224, 161, 148, 232, 99, 202, 108, 45, 218, 145, 93, 103,
//             64, 177, 105, 163, 115, 34, 11, 250, 31, 46, 213, 139, 205, 219, 194, 199, 175, 220,
//             79, 3, 24, 68, 23, 101, 90, 193, 206, 210, 154, 100, 174, 118, 75, 177, 40, 49, 191,
//             143, 71, 99, 216, 209, 213, 219, 8, 194, 185, 240, 21, 237, 232, 4, 164, 102, 35, 24,
//             8, 49, 150, 59, 231, 151, 5, 177, 113, 137, 188, 74, 159, 86, 235, 21, 197, 58, 192,
//             200, 141, 36, 22, 232, 32, 229, 188, 4, 231, 187, 232, 154, 165, 64, 98, 45, 101, 155,
//             81, 137, 180, 78, 215, 206, 64, 112, 184, 156, 39, 46, 42, 36, 247, 61, 70, 15, 234,
//             18, 57, 3, 91, 153, 220, 231, 247, 223, 122, 196, 24, 104, 250, 78, 142, 118, 67, 181,
//             96, 169, 20, 234, 58, 197, 63, 55, 114, 219, 233, 23, 223, 27, 69, 6, 115, 221, 244,
//             77, 82, 191, 9, 150, 178, 205, 168, 196, 218, 39, 107, 231, 32, 56, 92, 78, 19, 151,
//             21, 18, 123, 158, 163, 7, 245, 10, 237, 4, 231, 187, 232, 154, 165, 106, 226, 45, 101,
//             155, 81, 137, 180, 78, 215, 206, 64, 112, 184, 156, 39, 46, 42, 36, 247, 61, 70, 15,
//             234, 20, 185, 1, 140, 34, 11, 178, 173, 48, 7, 105, 77, 50, 87, 59, 37, 216, 148, 24,
//             223, 199, 163, 177, 236, 104, 234, 237, 132, 97, 92, 248, 10, 244, 20, 3, 24, 68, 23,
//             101, 90, 108, 206, 210, 154, 100, 174, 118, 75, 177, 40, 49, 191, 143, 71, 99, 216,
//             209, 213, 219, 8, 194, 185, 240, 21, 232, 232, 7, 51, 41, 162, 227, 93, 103, 1, 119,
//             157, 40, 0, 161, 176, 143, 13, 53, 64, 172, 166, 74, 71, 4, 177, 30, 109, 56, 217, 15,
//             172, 148, 119, 0, 1, 40, 225, 211, 37, 67, 227, 195, 75, 15, 184, 75, 35, 121, 152, 11,
//             206, 130, 181, 93, 52, 26, 222, 54, 225, 164, 206, 119, 98, 49, 127, 1, 245, 239, 229,
//             226, 164, 147, 237, 23, 92, 189, 192, 202, 171, 211, 97, 221, 103, 41, 20, 123, 42, 73,
//             255, 19, 182, 16, 44, 170, 91, 163, 241, 3, 122, 35, 184, 126, 224, 183, 84, 233, 162,
//             161, 139, 249, 241, 173, 181, 44, 40, 254, 122, 243, 31, 209, 50, 95, 136, 54, 66, 182,
//             182, 120, 86, 1, 200, 25, 46, 236, 199, 14, 120, 215, 166, 119, 142, 184, 187, 93, 210,
//             70, 99, 209, 109, 106, 7, 38, 13, 86, 113, 77, 40, 239, 101, 112, 47, 2, 80, 119, 11,
//             163, 93, 205, 89, 30, 136, 45, 112, 44, 157, 216, 11, 39, 148, 121, 29, 189, 236, 115,
//             33, 204, 138, 68, 69, 217, 20, 115, 169, 5, 14, 205, 72, 77, 54, 231, 218, 153, 95,
//             162, 175, 218, 232, 63, 190, 166, 244, 88, 215, 208, 135, 139, 66, 119, 107, 105, 209,
//             86, 146, 86, 139, 2, 52, 60, 90, 10, 156, 32, 31, 52, 138, 33, 75, 142, 77, 0, 167,
//             160, 116, 5, 177, 241, 191, 160, 205, 157, 11, 224, 168, 248, 210, 225, 35, 1, 8, 125,
//             71, 71, 5, 38, 16, 199, 27, 91, 50, 8, 54, 219, 166, 194, 250, 10, 202, 190, 145, 195,
//             160, 218, 188, 29, 73, 184, 247, 42, 44, 3, 221, 45, 93, 101, 219, 38, 37, 3, 232, 2,
//             98, 27, 149, 66, 245, 31, 23, 114, 149, 174, 202, 119, 96, 233, 114, 114, 46, 147, 166,
//             89, 244, 1, 166, 156, 186, 246, 241, 109, 41, 76, 181, 248, 23, 253, 193, 236, 87, 42,
//             52, 162, 91, 167, 141, 227, 164, 162, 247, 95, 64, 90, 59, 117, 210, 1, 143, 44, 156,
//             144, 94, 237, 240, 120, 79, 31, 189, 37, 133, 249, 195, 95, 23, 87, 168, 26, 6, 60,
//             175, 235, 164, 121, 7, 220, 167, 78, 247, 3, 157, 144, 158, 44, 242, 203, 233, 164, 71,
//             98, 98, 68, 87, 50, 103, 230, 182, 72, 240, 222, 223, 129, 196, 249, 204, 56, 77, 80,
//             64, 39, 35, 0, 32, 49, 240, 229, 228, 251, 68, 176, 221, 45, 123, 205, 240, 137, 20,
//             168, 248, 87, 111, 142, 170, 189, 190, 226, 99, 108, 192, 61, 135, 89, 138, 5, 100, 74,
//             55, 89, 173, 9, 154, 231, 111, 119, 138, 82, 126, 3, 197, 143, 28, 74, 78, 198, 99,
//             129, 126, 84, 31, 119, 224, 42, 247, 70, 75, 6, 249, 103, 53, 199, 171, 214, 151, 83,
//             116, 110, 0, 226, 78, 69, 116, 76, 146, 140, 180, 251, 2, 154, 84, 34, 123, 74, 210,
//             202, 193, 129, 242,
//         ];

//         let serialized_proof = StoneCompatibleSerializer::serialize_proof::<Fibonacci2ColsShifted<_>>(
//             &proof,
//             &pub_inputs,
//             &proof_options,
//         );
//         assert_eq!(serialized_proof, expected_bytes);
//     }

//     #[test]
//     fn test_serialization_compatible_with_stone_case_3() {
//         let mut trace = fibonacci_2_cols_shifted::compute_trace(FieldElement::from(12345), 512);

//         let claimed_index = 420;
//         let claimed_value = trace.get_row(claimed_index)[0];
//         let proof_options = ProofOptions {
//             blowup_factor: 64,
//             coset_offset: 3,
//             grinding_factor: 0,
//             fri_number_of_queries: 1,
//         };

//         let pub_inputs = fibonacci_2_cols_shifted::PublicInputs {
//             claimed_value,
//             claimed_index,
//         };

//         let proof = Prover::<Fibonacci2ColsShifted<_>>::prove(
//             &mut trace,
//             &pub_inputs,
//             &proof_options,
//             StoneProverTranscript::new(&pub_inputs.as_bytes()),
//         )
//         .unwrap();

//         let expected_bytes = [
//             109, 49, 221, 0, 3, 137, 116, 189, 229, 254, 12, 94, 58, 118, 95, 141, 220, 130, 42,
//             93, 243, 37, 79, 202, 133, 161, 149, 10, 224, 32, 140, 190, 239, 19, 119, 217, 232, 38,
//             44, 103, 224, 84, 37, 164, 175, 0, 176, 5, 228, 209, 131, 135, 35, 160, 245, 180, 101,
//             20, 17, 193, 139, 244, 214, 182, 4, 75, 226, 131, 251, 225, 219, 95, 239, 73, 57, 144,
//             157, 77, 141, 185, 96, 12, 72, 162, 220, 59, 28, 165, 125, 180, 59, 196, 125, 175, 147,
//             35, 3, 166, 113, 39, 138, 98, 152, 242, 130, 179, 98, 207, 75, 213, 4, 166, 18, 174,
//             48, 180, 163, 178, 171, 151, 243, 160, 172, 30, 19, 10, 81, 246, 7, 162, 17, 184, 194,
//             192, 238, 157, 46, 254, 115, 59, 129, 110, 64, 254, 132, 223, 32, 210, 127, 58, 127,
//             190, 163, 99, 231, 54, 160, 7, 227, 245, 1, 28, 36, 75, 6, 234, 197, 16, 123, 234, 123,
//             154, 110, 86, 44, 23, 105, 219, 66, 35, 196, 86, 208, 208, 157, 109, 255, 213, 31, 138,
//             123, 204, 0, 132, 81, 145, 33, 88, 101, 30, 95, 38, 58, 112, 158, 89, 220, 144, 206,
//             77, 119, 14, 188, 69, 181, 70, 203, 71, 219, 116, 215, 142, 82, 164, 61, 94, 225, 126,
//             73, 242, 43, 24, 127, 130, 247, 244, 127, 165, 142, 63, 223, 61, 62, 47, 58, 50, 239,
//             208, 158, 217, 33, 145, 153, 52, 14, 146, 60, 152, 42, 0, 41, 233, 231, 56, 238, 53,
//             148, 74, 50, 241, 205, 55, 228, 202, 41, 162, 225, 98, 181, 163, 113, 121, 186, 191,
//             251, 237, 32, 214, 95, 103, 181, 76, 217, 46, 234, 0, 133, 191, 106, 153, 5, 56, 85,
//             104, 188, 40, 43, 58, 21, 49, 220, 52, 59, 55, 57, 161, 205, 139, 63, 152, 192, 83,
//             136, 28, 1, 55, 103, 147, 6, 199, 241, 92, 63, 22, 56, 220, 43, 42, 4, 48, 68, 250, 77,
//             111, 227, 149, 14, 210, 133, 198, 182, 114, 56, 21, 64, 202, 15, 57, 86, 25, 159, 1,
//             35, 149, 231, 230, 219, 194, 124, 165, 228, 46, 208, 208, 164, 239, 100, 33, 111, 18,
//             227, 115, 181, 168, 180, 24, 131, 150, 226, 64, 87, 166, 222, 137, 57, 189, 96, 66, 46,
//             148, 142, 4, 142, 243, 124, 107, 115, 149, 111, 232, 87, 214, 39, 169, 76, 71, 189, 27,
//             56, 42, 148, 156, 18, 98, 128, 83, 135, 142, 177, 240, 195, 181, 221, 74, 3, 246, 155,
//             238, 56, 177, 216, 10, 65, 176, 176, 212, 50, 153, 84, 128, 133, 72, 64, 210, 71, 10,
//             48, 15, 164, 237, 45, 117, 5, 187, 184, 110, 119, 11, 176, 88, 244, 33, 128, 177, 152,
//             84, 180, 245, 96, 179, 0, 229, 209, 249, 139, 216, 125, 21, 8, 97, 40, 101, 126, 12,
//             86, 7, 58, 198, 234, 37, 156, 164, 28, 147, 117, 228, 144, 110, 220, 76, 177, 162, 6,
//             197, 149, 26, 240, 74, 208, 137, 170, 135, 203, 150, 205, 215, 51, 209, 52, 226, 185,
//             209, 170, 24, 64, 229, 90, 4, 210, 205, 115, 74, 164, 106, 7, 92, 83, 16, 123, 27, 14,
//             166, 87, 98, 205, 51, 110, 212, 153, 51, 231, 14, 51, 252, 2, 143, 196, 251, 26, 139,
//             172, 22, 53, 159, 102, 132, 0, 72, 76, 107, 117, 128, 25, 180, 198, 189, 189, 5, 45,
//             174, 88, 212, 140, 156, 71, 133, 166, 85, 44, 208, 109, 241, 218, 237, 126, 135, 159,
//             245, 4, 252, 211, 29, 148, 213, 86, 146, 213, 64, 168, 202, 144, 142, 118, 62, 66, 222,
//             252, 200, 165, 10, 226, 2, 159, 80, 239, 120, 74, 233, 33, 183, 135, 240, 9, 224, 18,
//             8, 242, 173, 64, 189, 254, 101, 25, 116, 224, 85, 71, 75, 28, 107, 156, 197, 233, 58,
//             215, 153, 14, 110, 80, 2, 254, 31, 10, 152, 9, 182, 44, 18, 128, 22, 64, 20, 175, 135,
//             117, 227, 111, 88, 148, 211, 188, 21, 233, 222, 56, 64, 219, 219, 57, 0, 242, 47, 240,
//             142, 172, 143, 46, 26, 205, 28, 28, 204, 175, 7, 16, 55, 157, 18, 59, 85, 96, 12, 161,
//             195, 113, 218, 83, 47, 191, 147, 91, 134, 244, 105, 179, 250, 211, 242, 252, 220, 59,
//             165, 146, 120, 202, 78, 196, 59, 169, 9, 10, 215, 21, 233, 236, 62, 126, 42, 108, 101,
//             150, 41, 198, 240, 102, 175, 43, 7, 117, 175, 8, 173, 182, 110, 124, 184, 31, 68, 251,
//             58, 156, 75, 22, 152, 171, 231, 122, 174, 65, 152, 228, 111, 145, 59, 97, 216, 8, 27,
//             170, 190, 95, 60, 179, 225, 115, 7, 87, 179, 156, 39, 151, 104, 248, 150, 153, 45, 198,
//             251, 253, 90, 22, 78, 15, 184, 68, 216, 106, 19, 24, 5, 136, 206, 229, 175, 254, 173,
//             177, 240, 221, 3, 156, 244, 92, 5, 237, 51, 202, 54, 101, 247, 95, 25, 234, 164, 217,
//             81, 5, 98, 193, 159, 81, 27, 161, 210, 195, 10, 61, 163, 168, 196, 190, 166, 31, 153,
//             103, 44, 57, 82, 245, 233, 21, 163, 5, 255, 63, 98, 250, 230, 134, 239, 130, 241, 40,
//             12, 66, 150, 73, 61, 231, 25, 155, 136, 16, 177, 131, 254, 17, 199, 46, 55, 83, 209,
//             203, 148, 20, 253, 89, 194, 121, 212, 156, 46, 42, 211, 13, 209, 251, 66, 136, 118, 51,
//             218, 166, 53, 97, 13, 162, 17, 58, 241, 225, 49, 231, 110, 253, 23, 186, 100, 143, 106,
//             209, 98, 176, 88, 131, 128, 22, 171, 147, 210, 79, 79, 2, 185, 240, 232, 248, 244, 216,
//             102, 170, 223, 42, 221, 209, 223, 88, 168, 52, 161, 250, 144, 138, 51, 57, 14, 225, 94,
//             44, 92, 159, 185, 161, 161, 23, 39, 102, 155, 37, 60, 158, 87, 223, 77, 117, 116, 86,
//             16, 207, 223, 9, 147, 47, 107, 172, 83, 246, 7, 133, 107, 110, 59, 167, 172, 203, 229,
//             114, 13, 183, 168, 172, 135, 192, 136, 80, 119, 67, 14, 254, 168, 56, 252, 111, 46, 13,
//             45, 166, 94, 99, 87, 193, 50, 17, 196, 0, 113, 135, 60, 104, 213, 3, 147, 151, 69, 123,
//             184, 247, 86, 208, 191, 63, 49, 124, 39, 153, 208, 240, 78, 91, 251, 222, 197, 237,
//             101, 221, 226, 189, 61, 1, 223, 152, 242, 135, 190, 30, 33, 118, 28, 134, 96, 255, 2,
//             3, 82, 102, 242, 105, 129, 86, 250, 249, 81, 246, 78, 207, 234, 47, 49, 44, 6, 1, 209,
//             121, 88, 181, 75, 98, 7, 171, 77, 227, 63, 203, 30, 108, 141, 33, 153, 2, 105, 125,
//             163, 197, 212, 90, 23, 87, 50, 97, 60, 102, 158, 230, 76, 9, 115, 86, 216, 215, 16, 41,
//             20, 47, 34, 20, 221, 144, 161, 102, 251, 212, 29, 24, 77, 76, 54, 95, 227, 133, 112,
//             134, 113, 197, 181, 151, 16, 103, 221, 115, 146, 226, 114, 240, 147, 205, 155, 155, 86,
//             216, 168, 102, 56, 73, 115, 164, 174, 76, 140, 62, 82, 221, 172, 69, 127, 175, 68, 15,
//             147, 11, 43, 162, 106, 224, 3, 137, 66, 240, 52, 244, 138, 231, 45, 192, 40, 124, 32,
//             166, 117, 7, 148, 179, 148, 142, 82, 165, 240, 211, 183, 159, 248, 144, 67, 85, 229,
//             10, 202, 39, 181, 181, 102, 190, 31, 198, 67, 124, 154, 3, 253, 122, 78, 56, 107, 195,
//             216, 243, 27, 85, 12, 18, 216, 134, 123, 107, 113, 210, 28, 125, 136, 75, 135, 213,
//             194, 210, 102, 183, 254, 69, 92, 157, 221, 119, 71, 105, 218, 10, 238, 112, 57, 203,
//             249, 35, 19, 39, 97, 37, 188, 192, 242, 83, 199, 255, 0, 117, 70, 193, 4, 85, 195, 164,
//             14, 84, 194, 117, 30, 169, 247, 101, 164, 228, 170, 210, 136, 236, 252, 88, 3, 80, 239,
//             109, 10, 250, 224, 169, 142, 110, 119, 160, 78, 188, 154, 155, 105, 88, 197, 81, 165,
//             90, 33, 139, 91, 31, 183, 236, 204, 100, 17, 183, 214, 33, 104, 243, 8, 25, 251, 149,
//             95, 247, 98, 248, 23, 7, 37, 128, 202, 37, 67, 107, 16, 18, 228, 252, 175, 55, 9, 97,
//             142, 228, 163, 25, 196, 247, 92, 61, 161, 20, 218, 252, 80, 215, 22, 83, 30, 246, 167,
//             75, 109, 44, 118, 81, 74, 26, 27, 180, 35, 115, 206, 175, 229, 90, 242, 120, 74, 164,
//             96, 54, 167, 67, 41, 219, 129, 224, 28, 52, 119, 29, 4, 11, 72, 86, 28, 7, 74, 76, 52,
//             179, 185, 57, 71, 134, 103, 63, 57, 168, 244, 194, 194, 205, 158, 233, 203, 155, 250,
//             118, 18, 146, 84, 173, 87, 201, 49, 146, 135, 113, 254, 63, 199, 36, 227, 189, 241,
//             197, 204, 119, 50, 10, 253, 201, 207, 9, 149, 79, 218, 123, 149, 26, 141, 53, 177, 179,
//             68, 183, 13, 158, 212, 231, 236, 212, 188, 192, 201, 129, 161, 121, 63, 225, 161, 104,
//             154, 203, 221, 53, 171, 235, 154, 137, 254, 247, 95, 215, 23, 109, 83, 148, 46, 160,
//             77, 164, 166, 156, 157, 27, 38, 111, 30, 127, 243, 163, 104, 251, 95, 15, 122, 132, 65,
//             88, 201, 15, 185, 146, 151, 169, 148, 184, 180, 44, 198, 244, 243, 170, 15, 217, 170,
//             157, 126, 163, 201, 71, 221, 97, 138, 4, 2, 178, 39, 118, 28, 107, 55, 24, 8, 38, 127,
//             160, 68, 204, 40, 139, 172, 121, 229, 232, 158, 197, 74, 241, 49, 224, 75, 127, 223,
//             150, 31, 221, 154, 209, 209, 152, 148, 206, 239, 63, 228, 54, 141, 73, 239, 241, 84,
//             43, 141, 223, 155, 160, 240, 208, 129, 125, 152, 209, 191, 84, 11, 5, 187, 0, 170, 82,
//             169, 142, 76, 191, 30, 221, 84, 80, 85, 30, 56, 127, 77, 37, 102, 212, 247, 168, 64,
//             201, 96, 85, 253, 58, 56, 106, 91, 99, 208, 16, 69, 15, 159, 57, 149, 124, 215, 138, 2,
//             143, 145, 0, 84, 233, 123, 126, 237, 240, 60, 107, 148, 1, 101, 199, 24, 180, 211, 147,
//             152, 98, 32, 76, 154, 112, 35, 187, 21, 72, 186, 22, 128, 171, 48, 179, 120, 132, 15,
//             118, 107, 103, 161, 76, 83, 216, 232, 22, 68, 203, 109, 26, 146, 160, 183, 39, 225, 43,
//             187, 121, 209, 176, 223, 62, 117, 154, 70, 218, 179, 56, 226, 186, 133, 203, 244, 25,
//             206, 121, 195, 33, 107, 43, 220, 183, 192, 194, 70, 157, 122, 236, 45, 93, 120, 252,
//             248, 17, 245, 187, 196, 57, 151, 50, 153, 151, 109, 52, 120, 229, 244, 193, 34, 219,
//             251, 244, 167, 245, 195, 171, 248, 100, 9, 126, 141, 121, 16, 126, 134, 105, 107, 129,
//             79, 237, 140, 144, 220, 219, 122, 51, 231, 78, 137, 36, 8, 218, 220, 34, 149, 198, 142,
//             240, 173, 27, 6, 207, 72, 131, 65, 155, 150, 106, 49, 193, 239, 160, 80, 92, 42, 149,
//             182, 43, 176, 230, 27, 245, 49, 83, 131, 61, 194, 116, 24, 215, 79, 48, 230, 78, 165,
//             79, 146, 16, 70, 134, 210, 213, 208, 210, 14, 200, 147, 148, 54, 108, 154, 155, 249,
//             71, 250, 199, 14, 249, 151, 234, 17, 170, 85, 201, 59, 40, 40, 251, 74, 54, 198, 250,
//             57, 221, 230, 132, 25, 201, 197, 205, 137, 200, 153, 255, 44, 79, 241, 83, 195, 206,
//             144, 12, 31, 251, 87, 99, 31, 254, 146, 82, 214, 51, 10, 80, 26, 30, 184, 224, 65, 137,
//             100, 61, 128, 193, 225, 85, 253, 192, 236, 29, 213, 90, 213, 197, 142, 47, 245, 243,
//             93, 192, 159, 235, 203, 206, 90, 231, 62, 52, 59, 254, 134, 213, 179, 164, 196, 239,
//             142, 28, 242, 185, 18, 100, 211, 39, 57, 206, 55, 107, 6, 223, 53, 127, 85, 196, 175,
//             202, 226, 83, 23, 58, 7, 131, 55, 49, 241, 47, 146, 75, 95, 131, 61, 30, 0, 201, 216,
//             226, 73, 15, 46, 156, 255, 3, 186, 188, 131, 118, 74, 228, 112, 156, 31, 41, 182, 25,
//             184, 126, 254, 35, 63, 132, 216, 83, 121, 200, 232, 213, 131, 208, 66, 34, 145, 114,
//             109, 109, 51, 174, 164, 89, 152, 45, 94, 205, 231, 136, 125, 201, 154, 104, 14, 178,
//             68, 126, 116, 246, 165, 127, 2, 175, 198, 98, 187, 95, 66, 220, 190, 132, 223, 75, 145,
//             173, 26, 84, 194, 177, 171, 144, 33, 48, 172, 94, 207, 22, 186, 146, 161, 29, 23, 174,
//             1, 138, 169, 178, 43, 126, 52, 58, 236, 47, 114, 147, 69, 86, 246, 196, 50, 100, 45,
//             94, 179, 93, 17, 115, 166, 83, 229, 220, 206, 83, 221, 235, 89, 247, 239, 177, 38, 235,
//             141, 64, 144, 240, 255, 127, 99, 132, 227, 183, 22, 174, 175, 71, 81, 236, 34, 191, 50,
//             107, 101, 113, 96, 88, 175, 23, 191, 219, 184, 69, 72, 164, 187, 168, 41, 103, 150, 67,
//             9, 153, 217, 25, 156, 235, 69, 144, 44, 180, 171, 70, 193, 209, 127, 10, 60, 203, 125,
//             154, 192, 65, 176, 80, 70, 51, 8, 154, 150, 6, 58, 210, 220, 6, 123, 123, 166, 141, 8,
//             167, 19, 93, 144, 21, 247, 104, 86, 255, 126, 188, 137, 59, 137, 119, 94, 153, 27, 246,
//             192, 11, 87, 153, 232, 34, 126, 167, 245, 15, 69, 188, 36, 41, 183, 230, 216, 179, 195,
//             21, 234, 109, 240, 4, 64, 72, 250, 43, 54, 72, 16, 69, 162, 235, 97, 59, 211, 123, 247,
//             59, 215, 110, 177, 7, 159, 210, 14, 25, 36, 123, 9, 232, 160, 211, 130, 213, 22, 254,
//             43, 151, 90, 14, 118, 78, 108, 248, 31, 118, 90, 243, 203, 163, 20, 109, 35, 30, 110,
//             65, 155, 182, 98, 64, 161, 15, 144, 89, 66, 63, 234, 111, 90, 184, 20, 175, 114, 103,
//             254, 203, 103, 95, 7, 128, 147, 58, 19, 73, 102, 145, 125, 98, 12, 87, 121, 57, 117,
//             114, 53, 170, 232, 50, 35, 228, 145, 21, 3, 152, 118, 250, 78, 24, 6, 140, 15, 125,
//             118, 247, 183, 198, 103, 98, 20, 112, 8, 31, 60, 15, 33, 3, 235, 119, 184, 102, 137,
//             145, 228, 147, 115, 95, 34, 212, 207, 125, 251, 240, 16, 186, 247, 188, 9, 29, 179,
//             176, 194, 250, 115, 71, 59, 197, 244, 199, 142, 19, 25, 14, 254, 109, 182, 160, 10,
//             123, 202, 55, 120, 80, 108, 3, 115, 66, 65, 77, 190, 248, 35, 116, 208, 23, 74, 164,
//             158, 133, 51, 234, 23, 152, 122, 84, 62, 203, 162, 53, 9, 152, 117, 121, 72, 45, 85,
//             191, 236, 150, 54, 238, 150, 171, 151, 166, 214, 67, 91, 104, 42, 119, 51, 168, 168,
//             158, 115, 193, 41, 125, 217, 78, 176, 73, 36, 132, 83, 7, 240, 82, 52, 147, 37, 38,
//             145, 191, 75, 133, 216, 143, 0, 55, 233, 230, 25, 113, 243, 251, 230, 43, 78, 183, 170,
//             250, 145, 168, 234, 180, 68, 55, 233, 53, 166, 80, 174, 43, 137, 157, 46, 107, 90, 154,
//             44, 229, 27, 83, 77, 226, 203, 24, 203, 105, 254, 242, 175, 230, 96, 45, 189, 42, 17,
//             118, 163, 89, 172, 138, 250, 245, 105, 82, 166, 238, 46, 95, 13, 255, 123, 14, 58, 40,
//             103, 179, 52, 94, 41, 32, 249, 39, 36, 83, 198, 144, 205, 114, 103, 219, 159, 78, 10,
//             77, 249, 30, 245, 207, 9, 160, 137, 200, 126, 97, 63, 194, 134, 137, 62, 53, 210, 224,
//             81, 225, 101, 120, 35, 158, 44, 65, 156, 23, 103, 57, 81, 82, 137, 94, 34, 251, 210,
//             238, 81, 25, 227, 16, 142, 230, 220, 255, 225, 98, 206, 238, 189, 54, 231, 80, 148,
//             212, 244, 13, 129, 55, 109, 164, 214, 15, 38, 246, 130, 51, 187, 83, 116, 150, 164, 35,
//             224, 162, 131, 189, 173, 117, 225, 236, 23, 102, 111, 146, 69, 146, 116, 125, 33, 214,
//             232, 112, 241, 169, 247, 80, 150, 27, 38, 179, 32, 107, 245, 170, 93, 32, 204, 243, 13,
//             130, 190, 57, 125, 14, 64, 16, 151, 77, 238, 178, 24, 88, 240, 41, 191, 116, 166, 126,
//             224, 42, 65, 166, 208, 35, 95, 24, 3, 216, 234, 60, 156, 110, 225, 78, 43, 147, 24,
//             103, 154, 220, 185, 7, 5, 4, 61, 222, 22, 83, 143, 156, 15, 78, 144, 46, 252, 193, 225,
//             101, 140, 154, 168, 107, 113, 213, 104, 147, 160, 232, 77, 143, 57, 170, 196, 109, 152,
//             31, 232, 70, 103, 179, 72, 49, 19, 12, 227, 98, 118, 107, 187, 30, 170, 40, 54, 88,
//             202, 222, 242, 54, 96, 28, 57, 147, 131, 105, 103, 68, 147, 6, 189, 201, 222, 223, 135,
//             111, 85, 12, 239, 17, 7, 137, 188, 172, 230, 129, 36, 165, 140, 136, 55, 123, 116, 23,
//             48, 122, 247, 235, 81, 63, 43, 180, 192, 175, 1, 243, 80, 116, 158, 228, 228, 115, 131,
//             124, 1, 7, 236, 97, 3, 168, 220, 43, 95, 118, 146, 51, 156, 84, 0, 24, 131, 237, 83,
//             117, 89, 109, 216, 173, 196, 148, 232, 170, 111, 188, 147, 58, 133, 152, 207, 48, 210,
//             195, 237, 0, 139, 177, 188, 113, 41, 133, 146, 249, 190, 184, 228, 4, 22, 147, 89, 66,
//             50, 77, 59, 179, 215, 87, 234, 166, 211, 186, 154, 125, 151, 100, 206, 176, 84, 65,
//             232, 108, 35, 39, 8, 234, 195, 214, 6, 138, 132, 172, 164, 215, 16, 200, 213, 203, 99,
//             9, 69, 6, 74, 34, 144, 49, 7, 236, 156, 73, 142, 156, 168, 40, 7, 126, 134, 43, 28,
//             148, 255, 94, 170, 23, 20, 228, 114, 165, 152, 148, 184, 159, 58, 214, 66, 118, 106,
//             165, 69, 100, 105, 106, 185, 187, 119, 205, 52, 233, 185, 167, 111, 117, 164, 60, 14,
//             201, 70, 24, 216, 149, 52, 57, 229, 198, 225, 12, 129, 52, 49, 228, 146, 86, 65, 251,
//             215, 61, 91, 104, 1, 141, 148, 5, 230, 212, 21, 96, 107, 220, 252, 150, 128, 37, 34,
//             97, 68, 189, 141, 59, 224, 93, 150, 199, 249, 170, 64, 78, 61, 231, 169, 212, 171, 98,
//             64, 2, 247, 190, 105, 31, 113, 201, 165, 102, 216, 100, 128, 233, 190, 150, 14, 71, 40,
//             142, 37, 45, 213, 56, 212, 229, 59, 42, 88, 175, 9, 231, 220, 37, 28, 230, 171, 14,
//             122, 78, 96, 111, 179, 18, 239, 193, 88, 125, 155, 99, 178, 211, 190, 100, 122, 126,
//             203, 219, 83, 29, 235, 242, 129, 50, 201, 15, 175, 95, 45, 244, 217, 186, 242, 98, 71,
//             197, 169, 155, 63, 107, 59, 155, 88, 5, 224, 20, 238, 113, 147, 11, 93, 228, 112, 180,
//             69, 229, 154, 188, 118, 210, 52, 6, 238, 83, 70, 184, 34, 136, 7, 200, 187, 40, 144,
//             76, 106, 56, 179, 198, 253, 47, 46, 16, 90, 126, 240, 26, 116, 102, 96, 30, 91, 235,
//             120, 183, 80, 62, 102, 111, 191, 241, 184, 129, 55, 79, 142, 204, 46, 70, 196, 8, 228,
//             67, 56, 202, 252, 91, 72, 74, 89, 195, 203, 81, 172, 198, 48, 187, 224, 77, 223, 216,
//             208, 31, 25, 95, 107, 141, 104, 101, 92, 190, 31, 223, 164, 221, 155, 180, 33, 202,
//             248, 217, 182, 253, 56, 132, 220, 212, 10, 18, 53, 192, 155, 232, 22, 68, 74, 62, 180,
//             253, 77, 21, 199, 56, 60, 101, 201, 215, 218, 221, 133, 229, 199, 82, 230, 190, 74, 69,
//             236, 226, 12, 203, 3, 216, 125, 254, 176, 176, 91, 39, 67, 182, 25, 47, 45, 58, 13, 38,
//             229, 130, 162, 255, 117, 80, 141, 188, 23, 87, 96, 101, 35, 128, 124, 112, 124, 29,
//             121, 114, 89, 253, 244, 161, 22, 231, 229, 178, 30, 99, 34, 220, 28, 172, 61, 229, 40,
//             19, 195, 32, 235, 203, 174, 78, 124, 131, 68, 31, 139, 142, 156, 38, 255, 44, 166, 171,
//             9, 196, 173, 49, 156, 16, 208, 244, 227, 30, 83, 150, 82, 61, 109, 200, 237, 163, 171,
//             188, 121, 118, 149, 61, 245, 224, 62, 179, 139, 217, 62, 245, 3, 110, 216, 142, 207,
//             17, 125, 51, 198, 124, 60, 46, 240, 26, 119, 215, 14, 68, 142, 150, 210, 197, 110, 107,
//             204, 128, 43, 214, 90, 187, 142, 60, 51, 56, 227, 239, 17, 13, 80, 199, 101, 8, 50,
//             122, 1, 190, 158, 212, 115, 61, 197, 51, 50, 235, 184, 200, 4, 250, 244, 203, 201, 222,
//             206, 252, 138, 34, 8, 16, 148, 198, 112, 157, 103, 243, 74, 209, 244, 186, 198, 244,
//             21, 229, 168, 246, 211, 151, 87, 7, 71, 253, 147, 109, 108, 78, 62, 148, 129, 220, 1,
//             104, 76, 69, 95, 81, 99, 13, 65, 108, 42, 176, 94, 189, 120, 205, 48, 112, 207, 188,
//             186, 93, 250, 13, 168, 247, 53, 112, 56, 227, 140, 246, 130, 229, 36, 200, 253, 212,
//             189, 217, 128, 107, 163, 78, 78, 90, 137, 44, 20, 0, 194, 70, 89, 185, 230, 32, 161,
//             37, 134, 142, 94, 155, 99, 18, 117, 154, 146, 129, 178, 53, 75, 113, 214, 116, 30, 218,
//             53, 180, 61, 237, 204, 145, 98, 132, 255, 167, 27, 45, 6, 60, 23, 50, 254, 108, 76, 42,
//             93, 131, 90, 209, 55, 0, 139, 71, 242, 18, 42, 225, 247, 183, 240, 147, 227, 254, 164,
//             19, 184, 41, 3, 112, 41, 191, 21, 140, 43, 230, 124, 207, 193, 233, 123, 45, 79, 195,
//             231, 59, 2, 181, 216, 54, 205, 186, 103, 203, 67, 10, 132, 122, 62, 203, 249, 223, 174,
//             234, 104, 130, 42, 32, 32, 94, 54, 135, 161, 186, 14, 34, 125, 102, 219, 251, 209, 221,
//             221, 125, 232, 18, 223, 208, 175, 206, 133, 160, 11, 228, 119, 142, 197, 232, 21, 10,
//             67, 68, 132, 11, 182, 80, 158, 67, 254, 252, 90, 15, 113, 112, 15, 60, 238, 90, 147,
//             207, 140, 90, 162, 211, 240, 193, 215, 182, 74, 101, 220, 176, 7, 168, 217, 104, 189,
//             59, 228, 142, 204, 180, 13, 203, 209, 83, 41, 148, 179, 159, 160, 250, 91, 79, 207,
//             233, 40, 4, 253, 168, 98, 209, 63, 58, 150, 164, 222, 133, 127, 213, 98, 158, 58, 14,
//             207, 208, 218, 25, 71, 224, 191, 62, 178, 234, 214, 111, 105, 0, 17, 167, 226, 93, 161,
//             245, 194, 161, 164, 156, 61, 174, 86, 76, 143, 78, 180, 242, 198, 189, 40, 27, 11, 119,
//             165, 160, 115, 105, 8, 68, 133, 214, 163, 223, 64, 137, 185, 146, 82, 106, 38, 76, 144,
//             61, 148, 57, 245, 29, 254, 246, 193, 9, 33, 84, 17, 84, 226, 201, 68, 103, 97, 220,
//             221, 25, 237, 18, 242, 72, 56, 226, 26, 140, 128, 82, 224, 238, 85, 137, 243, 204, 53,
//             148, 134, 103, 204, 93, 19, 253, 244, 200, 232, 100, 184, 242, 204, 155, 106, 154, 29,
//             46, 168, 152, 246, 179, 85, 253, 60, 243, 159, 62, 130, 41, 55, 131, 110, 146, 68, 137,
//             187, 75, 56, 9, 39, 10, 136, 120, 165, 205, 202, 76, 248, 123, 127, 159, 91, 51, 205,
//             165, 76, 241, 68, 25, 38, 138, 166, 228, 98, 242, 234, 155, 147, 252, 240, 208, 145,
//             23, 51, 115, 165, 22, 139, 227, 18, 81, 177, 10, 109, 86, 47, 179, 4, 198, 87, 207,
//             141, 110, 129, 53, 221, 82, 182, 220, 12, 127, 17, 191, 138, 130, 95, 80, 239, 36, 90,
//             185, 3, 128, 221, 252, 150, 126, 234, 79, 35, 89, 15, 175, 197, 74, 129, 167, 195, 243,
//             123, 23, 179, 96, 62, 217, 158, 90, 49, 137, 144, 199, 9, 129, 91, 66, 28, 29, 24, 245,
//             37, 16, 130, 47, 127, 16, 244, 220, 210, 248, 45, 98, 59, 90, 108, 153, 39, 37, 14, 22,
//             72, 227, 132, 1, 72, 2, 53, 42, 150, 201, 241, 221, 21, 41, 207, 135, 191, 50, 217, 11,
//             113, 73, 85, 139, 221, 233, 222, 250, 176, 230, 238, 67, 61, 26, 185, 162, 113, 139,
//             10, 189, 31, 151, 154, 199, 83, 47, 229, 68, 33, 248, 228, 32, 193, 252, 30, 134, 203,
//             16, 50, 157, 223, 66, 40, 207, 113, 1, 182, 195, 6, 66, 35, 187, 28, 20, 53, 118, 168,
//             184, 81, 202, 34, 48, 155, 135, 157, 205, 227, 227, 187, 4, 195, 139, 2, 15, 15, 174,
//             222, 109, 249, 146, 214, 77, 48, 63, 235, 45, 156, 211, 148, 188, 83, 45, 141, 250,
//             253, 71, 38, 83, 145, 199, 161, 210, 213, 169, 37, 253, 175, 81, 52, 133, 76, 109, 54,
//             107, 153, 205, 56, 74, 252, 58, 214, 195, 124, 19, 207, 237, 106, 205, 70, 165, 79, 28,
//             182, 202, 45, 254, 55, 118, 16, 17, 3, 21, 179, 64, 215, 69, 97, 216, 41, 183, 144,
//             191, 194, 113, 105, 206, 6, 180, 78, 139, 209, 44, 153, 89, 39, 248, 42, 128, 178, 19,
//             125, 246, 177, 221, 39, 89, 240, 200, 158, 25, 51, 162, 105, 42, 11, 254, 156, 62, 255,
//             183, 47, 228, 161, 87, 75, 132, 11, 107, 45, 45, 160, 169, 115, 73, 0, 14, 163,
//         ];

//         let serialized_proof = StoneCompatibleSerializer::serialize_proof::<Fibonacci2ColsShifted<_>>(
//             &proof,
//             &pub_inputs,
//             &proof_options,
//         );
//         assert_eq!(serialized_proof, expected_bytes);
//     }

//     #[test]
//     fn test_serialization_compatible_with_stone_4() {
//         let mut trace = fibonacci_2_cols_shifted::compute_trace(FieldElement::one(), 4);

//         let claimed_index = 2;
//         let claimed_value = trace.get_row(claimed_index)[0];
//         let proof_options = ProofOptions {
//             blowup_factor: 2,
//             coset_offset: 3,
//             grinding_factor: 0,
//             fri_number_of_queries: 2,
//         };

//         let pub_inputs = fibonacci_2_cols_shifted::PublicInputs {
//             claimed_value,
//             claimed_index,
//         };

//         let proof = Prover::<Fibonacci2ColsShifted<_>>::prove(
//             &mut trace,
//             &pub_inputs,
//             &proof_options,
//             StoneProverTranscript::new(&pub_inputs.as_bytes()),
//         )
//         .unwrap();

//         let expected_bytes = [
//             9, 161, 59, 243, 85, 60, 44, 155, 163, 203, 128, 147, 203, 253, 93, 16, 137, 42, 94,
//             225, 173, 254, 120, 1, 43, 167, 254, 15, 49, 148, 2, 50, 47, 159, 254, 83, 209, 26, 91,
//             113, 237, 157, 107, 252, 35, 130, 117, 22, 155, 181, 217, 11, 145, 208, 53, 201, 14,
//             148, 19, 247, 19, 105, 239, 108, 2, 17, 212, 103, 137, 76, 186, 20, 31, 118, 21, 139,
//             122, 39, 100, 197, 112, 11, 188, 227, 236, 15, 127, 186, 231, 187, 158, 137, 41, 180,
//             233, 36, 3, 145, 50, 37, 100, 198, 0, 152, 68, 30, 64, 111, 78, 211, 119, 49, 142, 180,
//             178, 74, 59, 150, 209, 45, 211, 159, 6, 216, 205, 161, 255, 142, 0, 104, 163, 169, 38,
//             160, 79, 95, 195, 96, 46, 20, 45, 189, 161, 181, 95, 133, 26, 124, 224, 44, 153, 119,
//             121, 29, 187, 126, 125, 161, 4, 45, 2, 1, 113, 106, 28, 174, 255, 138, 4, 34, 227, 191,
//             33, 203, 60, 20, 34, 36, 33, 8, 44, 53, 250, 177, 127, 59, 157, 229, 179, 87, 165, 58,
//             0, 0, 83, 98, 7, 41, 90, 187, 198, 80, 159, 250, 57, 252, 211, 64, 233, 110, 223, 155,
//             56, 189, 215, 57, 80, 161, 169, 246, 65, 133, 129, 132, 129, 233, 154, 204, 187, 178,
//             244, 76, 12, 9, 30, 113, 105, 206, 46, 192, 68, 96, 27, 72, 94, 126, 101, 253, 63, 94,
//             10, 89, 116, 120, 31, 123, 5, 224, 161, 148, 232, 99, 202, 108, 45, 218, 145, 93, 103,
//             64, 177, 105, 163, 115, 34, 11, 250, 31, 46, 213, 139, 205, 219, 194, 199, 175, 220,
//             79, 7, 51, 41, 162, 227, 93, 103, 1, 119, 157, 40, 0, 161, 176, 143, 13, 53, 64, 172,
//             166, 74, 71, 4, 177, 30, 109, 56, 217, 15, 172, 148, 119, 0, 1, 40, 225, 211, 37, 67,
//             227, 195, 75, 15, 184, 75, 35, 121, 152, 11, 206, 130, 181, 93, 52, 26, 222, 54, 225,
//             164, 206, 119, 98, 49, 127, 1, 245, 239, 229, 226, 164, 147, 237, 23, 92, 189, 192,
//             202, 171, 211, 97, 221, 103, 41, 20, 123, 42, 73, 255, 19, 182, 16, 44, 170, 91, 163,
//             241, 3, 122, 35, 184, 126, 224, 183, 84, 233, 162, 161, 139, 249, 241, 173, 181, 44,
//             40, 254, 122, 243, 31, 209, 50, 95, 136, 54, 66, 182, 182, 120, 86, 1, 200, 25, 46,
//             236, 199, 14, 120, 215, 166, 119, 142, 184, 187, 93, 210, 70, 99, 209, 109, 106, 7, 38,
//             13, 86, 113, 77, 40, 239, 101, 112, 47, 2, 80, 119, 11, 163, 93, 205, 89, 30, 136, 45,
//             112, 44, 157, 216, 11, 39, 148, 121, 29, 189, 236, 115, 33, 204, 138, 68, 69, 217, 20,
//             115, 169, 5, 14, 205, 72, 77, 54, 231, 218, 153, 95, 162, 175, 218, 232, 63, 190, 166,
//             244, 88, 215, 208, 135, 139, 66, 119, 107, 105, 209, 86, 146, 86, 139, 2, 52, 60, 90,
//             10, 156, 32, 31, 52, 138, 33, 75, 142, 77, 0, 167, 160, 116, 5, 177, 241, 191, 160,
//             205, 157, 11, 224, 168, 248, 210, 225, 35, 28, 249, 169, 95, 21, 155, 125, 184, 161,
//             209, 104, 28, 40, 157, 113, 186, 88, 83, 80, 52, 130, 162, 139, 20, 152, 253, 6, 236,
//             251, 188, 248, 74, 3, 157, 144, 158, 44, 242, 203, 233, 164, 71, 98, 98, 68, 87, 50,
//             103, 230, 182, 72, 240, 222, 223, 129, 196, 249, 204, 56, 77, 80, 64, 39, 35, 0, 32,
//             49, 240, 229, 228, 251, 68, 176, 221, 45, 123, 205, 240, 137, 20, 168, 248, 87, 111,
//             142, 170, 189, 190, 226, 99, 108, 192, 61, 135, 89, 138, 5, 100, 74, 55, 89, 173, 9,
//             154, 231, 111, 119, 138, 82, 126, 3, 197, 143, 28, 74, 78, 198, 99, 129, 126, 84, 31,
//             119, 224, 42, 247, 70, 75, 6, 249, 103, 53, 199, 171, 214, 151, 83, 116, 110, 0, 226,
//             78, 69, 116, 76, 146, 140, 180, 251, 2, 154, 84, 34, 123, 74, 210, 202, 193, 129, 242,
//             181, 142, 85, 140, 84, 138, 69, 121, 69, 23, 14, 219, 249, 133, 141, 242, 128, 253, 44,
//             159, 125, 93, 13, 89, 70, 107, 195, 118, 133, 114, 4, 202, 76, 185, 171, 27, 107, 95,
//             178, 68, 155, 72, 25, 53, 160, 89, 109, 77, 67, 112, 240, 99, 114, 26, 51, 240, 83,
//             134, 72, 157, 118, 238, 0, 156,
//         ];

//         let serialized_proof = StoneCompatibleSerializer::serialize_proof::<Fibonacci2ColsShifted<_>>(
//             &proof,
//             &pub_inputs,
//             &proof_options,
//         );
//         assert_eq!(serialized_proof, expected_bytes);
//     }

//     #[test]
//     fn test_serialization_compatible_with_stone_5() {
//         let mut trace = fibonacci_2_cols_shifted::compute_trace(FieldElement::one(), 128);

//         let claimed_index = 111;
//         let claimed_value = trace.get_row(claimed_index)[0];
//         let proof_options = ProofOptions {
//             blowup_factor: 4,
//             coset_offset: 3,
//             grinding_factor: 0,
//             fri_number_of_queries: 3,
//         };

//         let pub_inputs = fibonacci_2_cols_shifted::PublicInputs {
//             claimed_value,
//             claimed_index,
//         };

//         let proof = Prover::<Fibonacci2ColsShifted<_>>::prove(
//             &mut trace,
//             &pub_inputs,
//             &proof_options,
//             StoneProverTranscript::new(&pub_inputs.as_bytes()),
//         )
//         .unwrap();

//         let expected_bytes = [
//             68, 228, 98, 183, 28, 139, 62, 73, 131, 192, 34, 97, 48, 52, 113, 8, 60, 34, 63, 155,
//             94, 81, 98, 20, 135, 161, 25, 61, 234, 184, 129, 198, 144, 128, 202, 95, 113, 181, 23,
//             21, 159, 52, 240, 15, 152, 224, 44, 222, 4, 55, 135, 79, 36, 227, 217, 2, 99, 161, 149,
//             115, 30, 184, 45, 230, 4, 77, 37, 128, 110, 52, 56, 37, 193, 196, 32, 179, 243, 141,
//             31, 42, 204, 120, 141, 60, 220, 222, 222, 215, 24, 213, 46, 45, 197, 81, 12, 217, 6,
//             192, 96, 58, 36, 138, 6, 26, 193, 18, 57, 204, 116, 181, 43, 73, 201, 23, 56, 191, 204,
//             196, 103, 248, 81, 175, 7, 191, 7, 96, 94, 249, 1, 121, 108, 49, 11, 225, 107, 207,
//             252, 200, 206, 7, 175, 20, 138, 144, 147, 251, 124, 82, 97, 200, 54, 7, 85, 59, 200,
//             14, 98, 254, 17, 4, 7, 75, 53, 15, 137, 76, 197, 75, 96, 177, 216, 83, 24, 248, 153,
//             197, 35, 234, 125, 210, 179, 239, 38, 3, 147, 48, 3, 215, 224, 97, 158, 61, 3, 126, 7,
//             213, 168, 94, 76, 45, 126, 222, 108, 126, 98, 94, 181, 180, 118, 69, 73, 214, 126, 171,
//             171, 202, 3, 187, 25, 139, 137, 61, 168, 34, 228, 73, 162, 238, 201, 149, 8, 247, 182,
//             167, 58, 131, 254, 110, 116, 66, 36, 194, 73, 58, 230, 242, 105, 34, 119, 228, 51, 251,
//             56, 120, 109, 169, 9, 39, 243, 26, 57, 57, 60, 178, 75, 236, 199, 241, 184, 94, 150,
//             101, 202, 22, 99, 11, 6, 137, 207, 124, 220, 239, 95, 42, 177, 251, 103, 130, 56, 45,
//             74, 17, 203, 52, 106, 210, 111, 13, 90, 244, 147, 58, 194, 151, 31, 72, 196, 213, 43,
//             197, 113, 132, 5, 75, 120, 170, 187, 23, 187, 216, 67, 113, 205, 179, 18, 20, 92, 32,
//             204, 197, 25, 31, 253, 56, 204, 167, 77, 68, 218, 98, 186, 246, 237, 91, 67, 166, 157,
//             87, 193, 17, 28, 208, 248, 9, 158, 213, 14, 232, 27, 170, 208, 10, 28, 87, 85, 107, 16,
//             114, 130, 65, 211, 63, 185, 200, 32, 196, 50, 39, 234, 172, 62, 236, 108, 203, 61, 28,
//             143, 60, 61, 88, 54, 20, 228, 26, 62, 158, 49, 35, 64, 201, 182, 76, 166, 91, 237, 106,
//             66, 123, 144, 37, 119, 205, 156, 42, 142, 179, 6, 51, 80, 39, 144, 181, 147, 155, 195,
//             147, 115, 126, 133, 228, 85, 152, 154, 188, 114, 10, 9, 215, 176, 133, 207, 112, 52,
//             226, 238, 30, 74, 18, 1, 11, 196, 150, 186, 177, 35, 112, 40, 217, 137, 137, 207, 123,
//             104, 13, 239, 231, 201, 87, 108, 76, 75, 73, 12, 103, 0, 43, 168, 13, 40, 42, 5, 186,
//             109, 184, 112, 11, 56, 245, 168, 76, 222, 106, 188, 154, 219, 181, 69, 45, 158, 243,
//             89, 144, 86, 144, 232, 148, 25, 39, 89, 247, 95, 181, 7, 164, 13, 119, 129, 40, 59,
//             108, 21, 234, 12, 216, 35, 47, 214, 135, 136, 133, 146, 114, 106, 44, 87, 87, 239, 49,
//             161, 206, 102, 25, 33, 63, 5, 107, 72, 49, 191, 149, 185, 133, 160, 79, 228, 98, 219,
//             82, 69, 215, 230, 78, 145, 23, 161, 64, 200, 183, 50, 48, 71, 146, 173, 140, 39, 9, 7,
//             123, 153, 177, 94, 211, 89, 72, 83, 58, 26, 218, 17, 85, 196, 107, 97, 207, 15, 248,
//             122, 85, 194, 237, 25, 133, 54, 221, 54, 247, 35, 111, 0, 208, 142, 88, 50, 205, 210,
//             163, 184, 46, 195, 83, 185, 188, 104, 102, 247, 40, 141, 20, 55, 29, 25, 120, 135, 51,
//             69, 129, 224, 51, 150, 114, 3, 253, 93, 25, 49, 124, 36, 249, 212, 37, 126, 251, 218,
//             93, 50, 228, 222, 233, 229, 72, 116, 109, 3, 212, 79, 53, 195, 179, 134, 49, 176, 25,
//             3, 92, 61, 142, 129, 23, 152, 177, 153, 80, 223, 227, 47, 26, 39, 36, 22, 10, 96, 130,
//             211, 103, 128, 244, 243, 171, 100, 175, 26, 197, 156, 241, 2, 186, 125, 100, 79, 57,
//             119, 22, 36, 98, 178, 130, 49, 83, 132, 195, 21, 231, 138, 65, 235, 145, 236, 235, 73,
//             175, 183, 200, 78, 52, 95, 7, 1, 154, 106, 98, 38, 34, 125, 79, 194, 182, 75, 26, 36,
//             120, 239, 212, 90, 215, 225, 60, 57, 99, 159, 79, 145, 116, 235, 251, 142, 226, 56,
//             197, 3, 208, 28, 167, 217, 109, 97, 164, 78, 9, 221, 139, 241, 209, 198, 38, 234, 171,
//             166, 174, 244, 182, 228, 121, 71, 4, 219, 46, 146, 32, 153, 79, 0, 229, 11, 36, 224,
//             61, 63, 63, 223, 95, 190, 126, 79, 27, 233, 156, 232, 147, 3, 96, 233, 120, 240, 197,
//             168, 68, 173, 9, 87, 73, 150, 56, 230, 74, 62, 76, 205, 143, 234, 82, 186, 92, 209,
//             154, 181, 231, 146, 150, 218, 233, 192, 114, 135, 240, 252, 59, 28, 169, 254, 204, 37,
//             25, 50, 39, 199, 195, 40, 184, 197, 230, 68, 196, 171, 61, 133, 181, 140, 132, 116,
//             169, 211, 164, 3, 5, 17, 37, 149, 79, 160, 145, 107, 194, 58, 130, 226, 93, 38, 247,
//             17, 150, 146, 72, 188, 109, 75, 15, 170, 128, 97, 229, 188, 170, 188, 133, 103, 217,
//             153, 41, 64, 154, 159, 87, 167, 80, 240, 89, 123, 16, 152, 42, 23, 235, 165, 232, 71,
//             32, 101, 31, 18, 27, 79, 122, 243, 65, 83, 76, 90, 200, 108, 203, 252, 64, 53, 97, 230,
//             117, 194, 55, 249, 168, 98, 203, 12, 179, 20, 223, 151, 185, 253, 89, 29, 232, 86, 1,
//             26, 123, 245, 220, 240, 198, 90, 200, 187, 99, 80, 22, 100, 113, 163, 105, 109, 69, 87,
//             49, 150, 68, 124, 149, 68, 102, 62, 17, 139, 36, 205, 201, 119, 12, 47, 172, 148, 107,
//             234, 240, 95, 111, 193, 142, 215, 149, 216, 239, 133, 171, 180, 88, 68, 35, 128, 205,
//             214, 29, 34, 123, 166, 211, 173, 22, 129, 23, 116, 30, 79, 148, 38, 183, 196, 205, 233,
//             208, 166, 133, 158, 5, 61, 143, 89, 15, 119, 45, 160, 34, 100, 233, 242, 174, 246, 156,
//             28, 68, 157, 216, 96, 95, 144, 145, 188, 251, 88, 211, 67, 245, 224, 233, 154, 145, 75,
//             126, 207, 89, 206, 219, 207, 64, 79, 155, 172, 175, 211, 148, 237, 102, 130, 249, 13,
//             40, 229, 74, 140, 198, 170, 0, 153, 157, 83, 183, 177, 41, 219, 229, 16, 233, 69, 95,
//             239, 241, 93, 54, 219, 200, 204, 154, 81, 162, 234, 16, 164, 68, 147, 97, 213, 197,
//             180, 198, 243, 58, 92, 161, 203, 230, 106, 191, 110, 167, 115, 124, 216, 61, 251, 9,
//             16, 190, 50, 60, 230, 237, 2, 38, 27, 18, 85, 78, 225, 44, 251, 232, 48, 217, 96, 56,
//             25, 230, 224, 118, 25, 253, 198, 21, 78, 22, 153, 101, 159, 239, 209, 40, 98, 44, 36,
//             74, 251, 150, 171, 107, 83, 164, 200, 154, 113, 18, 162, 204, 179, 56, 170, 71, 253,
//             206, 68, 63, 92, 114, 207, 64, 213, 160, 36, 169, 121, 170, 226, 126, 37, 64, 53, 62,
//             16, 96, 187, 32, 14, 186, 29, 127, 178, 23, 199, 56, 133, 139, 164, 78, 79, 235, 158,
//             131, 189, 178, 8, 75, 197, 139, 88, 250, 133, 155, 95, 212, 56, 135, 122, 194, 247, 94,
//             49, 108, 65, 3, 182, 15, 43, 226, 153, 135, 142, 229, 178, 0, 93, 91, 116, 163, 228,
//             145, 112, 112, 78, 109, 23, 206, 245, 65, 10, 19, 45, 118, 96, 91, 184, 162, 217, 74,
//             135, 106, 233, 151, 97, 69, 106, 54, 22, 135, 48, 189, 165, 191, 9, 113, 238, 226, 16,
//             154, 162, 141, 15, 81, 110, 33, 61, 253, 218, 230, 112, 35, 2, 64, 253, 59, 89, 128,
//             85, 216, 157, 240, 241, 118, 248, 132, 182, 59, 137, 73, 171, 88, 152, 227, 205, 19,
//             220, 85, 60, 122, 87, 94, 150, 221, 2, 135, 13, 66, 187, 49, 197, 41, 101, 243, 246,
//             183, 197, 212, 15, 107, 193, 156, 220, 63, 123, 224, 16, 184, 134, 114, 73, 33, 26, 35,
//             110, 152, 2, 222, 134, 92, 157, 96, 123, 189, 210, 214, 78, 114, 52, 51, 33, 49, 124,
//             75, 224, 108, 130, 162, 20, 43, 193, 94, 229, 228, 174, 33, 162, 230, 2, 67, 193, 53,
//             92, 25, 154, 95, 29, 158, 67, 35, 255, 194, 83, 170, 26, 76, 53, 98, 2, 138, 27, 103,
//             6, 203, 183, 226, 48, 0, 22, 254, 0, 3, 94, 32, 155, 0, 191, 57, 203, 150, 21, 83, 133,
//             242, 130, 244, 91, 10, 83, 5, 113, 44, 10, 248, 142, 20, 73, 10, 71, 240, 157, 161, 5,
//             143, 247, 49, 252, 183, 33, 50, 137, 74, 197, 4, 42, 50, 113, 36, 206, 204, 213, 198,
//             100, 197, 80, 206, 40, 224, 114, 159, 59, 145, 231, 176, 80, 20, 252, 161, 230, 26,
//             114, 245, 38, 71, 150, 233, 168, 195, 228, 168, 78, 149, 0, 95, 55, 24, 159, 157, 76,
//             227, 4, 185, 5, 192, 102, 240, 103, 40, 9, 69, 72, 236, 44, 181, 98, 29, 58, 121, 63,
//             168, 190, 225, 27, 35, 59, 228, 230, 15, 227, 211, 90, 205, 108, 203, 228, 145, 219,
//             212, 254, 92, 100, 220, 50, 57, 62, 36, 194, 133, 80, 11, 243, 190, 227, 83, 152, 168,
//             50, 114, 146, 4, 134, 77, 236, 36, 235, 104, 113, 72, 140, 148, 53, 109, 109, 154, 203,
//             128, 21, 253, 192, 18, 44, 220, 150, 63, 118, 67, 30, 7, 0, 81, 226, 168, 92, 185, 135,
//             192, 152, 18, 213, 217, 120, 15, 145, 206, 194, 168, 56, 158, 189, 112, 143, 45, 11,
//             197, 229, 69, 245, 105, 217, 9, 97, 177, 107, 59, 144, 247, 23, 132, 50, 129, 75, 134,
//             125, 95, 78, 172, 175, 3, 193, 3, 247, 125, 11, 25, 125, 21, 198, 124, 108, 234, 120,
//             52, 104, 64, 204, 147, 109, 117, 108, 45, 3, 29, 163, 208, 221, 199, 64, 162, 214, 72,
//             56, 32, 221, 13, 220, 59, 239, 242, 232, 210, 113, 43, 75, 149, 170, 229, 221, 49, 161,
//             253, 78, 106, 113, 65, 227, 169, 99, 72, 225, 65, 62, 57, 127, 53, 167, 25, 73, 164, 0,
//             208, 84, 56, 86, 132, 13, 17, 83, 183, 27, 164, 26, 196, 193, 214, 195, 76, 176, 210,
//             135, 88, 151, 69, 253, 24, 47, 146, 0, 120, 3, 211, 113, 128, 191, 227, 235, 105, 198,
//             181, 240, 186, 255, 242, 196, 193, 216, 15, 192, 101, 165, 160, 19, 243, 52, 166, 254,
//             137, 11, 156, 192, 63, 70, 91, 251, 0, 229, 197, 209, 129, 198, 232, 97, 49, 238, 248,
//             141, 210, 9, 80, 14, 115, 251, 82, 235, 132, 209, 3, 123, 43, 25, 29, 59, 175, 204,
//             127, 144, 241, 61, 137, 123, 6, 130, 155, 200, 55, 190, 33, 194, 50, 48, 238, 239, 132,
//             118, 216, 63, 203, 178, 81, 227, 87, 184, 177, 147, 192, 254, 206, 134, 77, 2, 120, 58,
//             180, 95, 159, 37, 207, 64, 121, 101, 134, 179, 165, 105, 154, 212, 50, 195, 23, 39, 66,
//             190, 216, 32, 56, 224, 165, 191, 114, 84, 96, 155, 85, 71, 135, 46, 198, 47, 80, 151,
//             176, 94, 211, 249, 48, 134, 114, 110, 131, 32, 21, 12, 162, 245, 7, 186, 30, 199, 218,
//             204, 232, 115, 160, 85, 45, 0, 80, 227, 65, 212, 135, 143, 151, 84, 168, 237, 85, 38,
//             141, 154, 216, 217, 241, 77, 141, 113, 207, 196, 132, 156, 240, 130, 249, 118, 251, 61,
//             112, 4, 68, 121, 110, 140, 255, 49, 123, 233, 40, 222, 225, 213, 160, 81, 1, 236, 126,
//             136, 42, 145, 123, 239, 96, 215, 233, 81, 172, 231, 11, 138, 194, 116, 4, 54, 43, 117,
//             64, 159, 170, 166, 162, 143, 245, 175, 100, 116, 156, 227, 64, 0, 200, 44, 239, 120,
//             98, 68, 96, 27, 218, 61, 5, 62, 159, 225, 248, 250, 172, 77, 61, 190, 158, 84, 143, 48,
//             177, 68, 18, 225, 147, 106, 93, 5, 109, 237, 41, 154, 223, 225, 6, 45, 41, 48, 243,
//             184, 129, 161, 224, 67, 2, 34, 71, 61, 36, 96, 139, 217, 194, 107, 15, 224, 21, 144,
//             191, 0, 194, 170, 179, 77, 194, 164, 32, 255, 81, 252, 211, 137, 45, 90, 238, 151, 229,
//             25, 31, 24, 176, 150, 252, 202, 14, 176, 159, 232, 199, 115, 147, 184, 236, 254, 115,
//             18, 10, 113, 29, 15, 120, 125, 204, 34, 25, 14, 33, 103, 7, 118, 8, 120, 95, 3, 225,
//             123, 82, 87, 41, 136, 152, 150, 185, 1, 1, 192, 119, 80, 225, 235, 11, 171, 28, 189,
//             62, 157, 244, 240, 117, 104, 154, 69, 114, 58, 152, 188, 62, 185, 174, 151, 204, 194,
//             160, 36, 146, 90, 3, 181, 253, 91, 180, 196, 210, 102, 20, 237, 201, 45, 116, 22, 178,
//             36, 161, 133, 84, 39, 220, 181, 133, 235, 72, 128, 198, 186, 55, 160, 81, 36, 98, 140,
//             128, 131, 140, 199, 201, 142, 171, 106, 33, 210, 198, 66, 164, 9, 161, 46, 232, 216,
//             239, 244, 66, 88, 166, 111, 82, 201, 46, 53, 194, 117, 244, 34, 89, 156, 148, 180, 84,
//             49, 176, 38, 193, 5, 25, 132, 91, 210, 101, 13, 205, 228, 222, 84, 81, 160, 94, 34,
//             136, 192, 223, 119, 224, 16, 1, 74, 62, 42, 177, 100, 202, 183, 170, 108, 146, 201, 99,
//             175, 196, 55, 189, 170, 110, 233, 86, 223, 233, 249, 165, 68, 178, 26, 87, 49, 19, 212,
//             152, 40, 114, 175, 166, 145, 133, 202, 93, 153, 57, 19, 229, 216, 148, 234, 15, 46,
//             167, 5, 95, 51, 144, 242, 159, 28, 236, 158, 12, 206, 109, 165, 35, 105, 158, 238, 207,
//             197, 179, 150, 119, 151, 199, 21, 36, 101, 188, 116, 106, 240, 101, 159, 154, 194, 240,
//             39, 128, 152, 160, 178, 251, 56, 249, 195, 50, 113, 227, 202, 175, 5, 9, 249, 117, 148,
//             203, 104, 14, 169, 174, 197, 121, 245, 81, 140, 16, 129, 47, 255, 7, 125, 169, 239,
//             111, 235, 138, 243, 52, 63, 230, 187, 163, 234, 134, 184, 36, 136, 24, 181, 226, 243,
//             153, 8, 61, 242, 126, 123, 64, 245, 196, 11, 189, 149, 238, 56, 228, 248, 87, 19, 215,
//             198, 29, 145, 155, 118, 246, 120, 198, 170, 107, 1, 174, 81, 237, 113, 79, 100, 102,
//             237, 28, 10, 198, 210, 178, 250, 8, 138, 64, 184, 187, 13, 171, 107, 236, 127, 198, 41,
//             240, 158, 96, 243, 229, 191, 251, 102, 191, 202, 186, 90, 255, 54, 15, 172, 46, 135,
//             247, 116, 238, 184, 227, 57, 252, 227, 149, 219, 69, 92, 24, 245, 83, 49, 250, 130,
//             212, 115, 8, 166, 14, 145, 240, 119, 40, 147, 9, 247, 235, 232, 159, 65, 72, 204, 131,
//             132, 94, 6, 155, 127, 65, 84, 141, 54, 213, 93, 217, 118, 100, 175, 20, 55, 38, 12, 63,
//             109, 69, 113, 169, 95, 29, 227, 137, 105, 100, 255, 166, 229, 216, 1, 148, 154, 105,
//             227, 201, 229, 195, 134, 12, 251, 164, 76, 103, 227, 205, 19, 188, 141, 66, 24, 241,
//             59, 49, 178, 95, 11, 154, 240, 182, 83, 33, 7, 62, 98, 233, 175, 150, 136, 1, 137, 152,
//             28, 83, 237, 100, 236, 147, 131, 196, 119, 219, 143, 113, 0, 18, 232, 195, 92, 90, 191,
//             243, 64, 181, 240, 217, 13, 101, 72, 197, 237, 199, 138, 60, 21, 70, 144, 86, 178, 175,
//             145, 95, 76, 174, 43, 188, 233, 139, 161, 203, 91, 235, 178, 182, 225, 155, 23, 219,
//             42, 119, 143, 246, 211, 154, 55, 126, 69, 61, 8, 40, 189, 190, 92, 130, 85, 54, 143,
//             231, 191, 19, 243, 31, 192, 242, 150, 131, 195, 130, 111, 181, 32, 147, 190, 138, 172,
//             146, 222, 208, 74, 95, 209, 122, 185, 72, 1, 96, 71, 34, 116, 182, 82, 35, 126, 133,
//             51, 101, 17, 2, 83, 15, 149, 98, 49, 47, 172, 137, 223, 248, 89, 82, 186, 152, 89, 4,
//             60, 74, 7, 181, 205, 181, 14, 92, 40, 203, 101, 155, 35, 206, 112, 144, 102, 184, 219,
//             16, 87, 237, 188, 68, 169, 90, 189, 61, 201, 209, 192, 99, 154, 134, 56, 19, 102, 107,
//             25, 241, 237, 174, 105, 136, 199, 55, 111, 225, 58, 152, 112, 159, 125, 111, 81, 237,
//             21, 39, 54, 253, 120, 189, 220, 164, 138, 135, 25, 61, 182, 242, 116, 100, 159, 41, 27,
//             144, 78, 104, 71, 129, 64, 142, 5, 63, 97, 239, 225, 163, 229, 9, 143, 156, 101, 201,
//             237, 57, 80, 103, 1, 36, 135, 63, 61, 247, 89, 202, 132, 140, 177, 178, 98, 7, 151, 10,
//             240, 59, 152, 109, 124, 11, 28, 218, 94, 131, 163, 101, 71, 187, 17, 162, 35, 22, 22,
//             9, 237, 30, 120, 118, 15, 50, 179, 52, 50, 5, 183, 194, 137, 254, 74, 80, 158, 238,
//             236, 186, 186, 121, 197, 231, 114, 183, 27, 6, 238, 104, 30, 254, 130, 247, 149, 224,
//             129, 200, 162, 49, 206, 20, 197, 45, 206, 179, 118, 169, 128, 184, 157, 85, 212, 198,
//             192, 22, 208, 130, 116, 99, 218, 56, 14, 249, 204, 234, 50, 74, 12, 7, 155, 116, 192,
//             213, 201, 115, 2, 22, 203, 145, 45, 140, 35, 17, 156, 209, 62, 73, 171, 132, 196, 84,
//             231, 146, 76, 123, 253, 109, 178, 117, 99, 52, 6, 221, 115, 207, 36, 232, 3, 46, 133,
//             43, 190, 220, 39, 242, 179, 120, 205, 120, 221, 101, 249, 43, 226, 131, 59, 159, 214,
//             202, 202, 90, 133, 51, 135, 141, 179, 123, 114, 167, 4, 65, 149, 218, 30, 196, 221, 44,
//             249, 170, 179, 73, 171, 89, 59, 177, 194, 95, 94, 61, 62, 44, 172, 172, 94, 228, 92,
//             17, 68, 16, 157, 1, 119, 137, 91, 67, 52, 14, 235, 57, 92, 102, 224, 208, 182, 222, 43,
//             196, 69, 33, 14, 184, 75, 54, 5, 59, 127, 69, 156, 222, 116, 185, 224, 106, 233, 94,
//             189, 183, 90, 204, 198, 142, 51, 129, 134, 251, 108, 187, 61, 78, 135, 102, 59, 165,
//             49, 38, 102, 97, 38, 102, 160, 229, 238, 170, 20, 24, 214, 83, 0, 202, 171, 190, 74,
//             55, 243, 17, 36, 116, 241, 94, 26, 17, 55, 245, 71, 104, 130, 141, 251, 1, 47, 77, 74,
//             17, 113, 48, 43, 204, 233, 157, 153, 136, 42, 132, 58, 173, 224, 23, 245, 9, 220, 224,
//             49, 142, 179, 49, 0, 135, 25, 109, 220, 226, 252, 51, 213, 235, 150, 101, 30, 63, 128,
//             82, 112, 187, 90, 64, 81, 30, 32, 116, 75, 233, 128, 79, 246, 155, 95, 193, 172, 93,
//             193, 84, 146, 143, 172, 21, 172, 6, 41, 41, 213, 50, 26, 178, 196, 126, 186, 74, 45,
//             141, 211, 182, 92, 154, 186, 207, 99, 221, 21, 177, 107, 199, 7, 98, 231, 144, 252,
//             120, 113, 217, 138, 23, 234, 162, 127, 180, 130, 126, 159, 1, 17, 50, 36, 249, 181,
//             123, 236, 110, 59, 26, 38, 212, 137, 58, 2, 73, 224, 135, 178, 0, 31, 157, 255, 68, 49,
//             205, 4, 109, 154, 29, 232, 129, 92, 67, 158, 161, 213, 129, 113, 53, 231, 159, 234,
//             114, 104, 17, 79, 89, 214, 106, 55, 150, 252, 35, 39, 78, 186, 59, 6, 12, 215, 97, 45,
//             88, 146, 80, 229, 50, 153, 202, 83, 158, 58, 247, 124, 78, 52, 194, 26, 198, 250, 14,
//             192, 162, 161, 25, 146, 190, 171, 241, 7, 94, 203, 139, 69, 222, 13, 1, 57, 86, 64, 84,
//             159, 226, 195, 239, 40, 26, 4, 120, 56, 127, 123, 209, 35, 127, 95, 157, 15, 155, 6,
//             143, 124, 37, 211, 186, 113, 213, 29, 101, 209, 238, 20, 207, 127, 45, 90, 245, 44,
//             220, 94, 224, 57, 204, 96, 83, 1, 90, 132, 111, 221, 5, 210, 186, 230, 39, 151, 46, 57,
//             31, 96, 218, 59, 115, 108, 29, 78, 23, 55, 231, 88, 142, 61, 147, 44, 57, 9, 112, 84,
//             55, 136, 254, 87, 83, 214, 23, 173, 31, 156, 202, 26, 193, 84, 130, 158, 88, 208, 209,
//             118, 231, 92, 160, 51, 32, 210, 125, 5, 114, 230, 119, 152, 165, 153, 98, 76, 145, 143,
//             78, 35, 201, 207, 251, 9, 44, 167, 198, 112, 22, 54, 224, 178, 216, 201, 248, 217, 139,
//             103, 86, 83, 220, 71, 244, 164, 126, 22, 91, 122, 154, 205, 30, 4, 76, 248, 75, 200,
//             191, 201, 95, 209, 20, 107, 13, 70, 88, 212, 15, 33, 160, 178, 202, 221, 23, 159, 1,
//             115, 152, 141, 54, 105, 37, 188, 106, 216, 119, 188, 233, 128, 226, 25, 12, 101, 193,
//             171, 81, 34, 156, 229, 241, 99, 243, 146, 33, 89, 193, 48, 48, 134, 213, 134, 232, 209,
//             177, 91, 29, 82, 242, 106, 241, 216, 132, 39, 20, 166, 59, 199, 184, 187, 139, 174, 40,
//             171, 149, 158, 160, 163, 255, 210, 111, 24, 201, 96, 54, 190, 244, 214, 85, 200, 239,
//             61, 99, 124, 239, 244, 170, 247, 153, 202, 47, 20, 136, 236, 17, 58, 164, 17, 196, 171,
//             171, 7, 235, 126, 171, 148, 60, 19, 1, 205, 202, 6, 230, 164, 222, 254, 83, 237, 80,
//             32, 177, 77, 12, 67, 106, 39, 48, 156, 107, 178, 36, 72, 125, 131, 179, 165, 124, 40,
//             139, 172, 178, 1, 170, 7, 247, 141, 97, 68, 98, 180, 164, 54, 120, 128, 134, 192, 248,
//             3, 197, 136, 207, 82, 119, 185, 10, 106, 216, 84, 173, 87, 176, 0, 21, 151, 48, 220,
//             196, 109, 236, 149, 52, 82, 251, 14, 201, 97, 226, 75, 177, 52, 16, 249, 36, 158, 103,
//             210, 33, 191, 114, 98, 40, 235, 19, 219, 101, 88, 189,
//         ];

//         let serialized_proof = StoneCompatibleSerializer::serialize_proof::<Fibonacci2ColsShifted<_>>(
//             &proof,
//             &pub_inputs,
//             &proof_options,
//         );
//         assert_eq!(serialized_proof, expected_bytes);
//     }
// }
