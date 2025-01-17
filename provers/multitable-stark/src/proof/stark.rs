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
