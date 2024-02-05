use std::marker::PhantomData;
#[cfg(feature = "instruments")]
use std::time::Instant;

use lambdaworks_math::fft::cpu::bit_reversing::{in_place_bit_reverse_permute, reverse_index};
use lambdaworks_math::fft::errors::FFTError;

use lambdaworks_math::field::traits::{IsField, IsSubFieldOf};
use lambdaworks_math::traits::AsBytes;
use lambdaworks_math::{
    field::{element::FieldElement, traits::IsFFTField},
    polynomial::Polynomial,
};
use log::info;

#[cfg(feature = "parallel")]
use rayon::prelude::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

#[cfg(debug_assertions)]
use crate::debug::validate_trace;
use crate::fri;
use crate::proof::stark::{DeepPolynomialOpenings, PolynomialOpenings};
use crate::table::Table;
use crate::trace::{columns2rows, LDETraceTable};
use crate::transcript::IsStarkTranscript;

use super::config::{BatchedMerkleTree, Commitment};
use super::constraints::evaluator::ConstraintEvaluator;
use super::domain::Domain;
use super::fri::fri_decommit::FriDecommitment;
use super::grinding;
use super::proof::options::ProofOptions;
use super::proof::stark::{DeepPolynomialOpening, StarkProof};
use super::trace::TraceTable;
use super::traits::AIR;

/// A default STARK prover implementing `IsStarkProver`.
pub struct Prover<A: AIR> {
    phantom: PhantomData<A>,
}

impl<A: AIR> IsStarkProver<A> for Prover<A> {}

#[derive(Debug)]
pub enum ProvingError {
    WrongParameter(String),
}

/// A container for the intermediate results of the commitments to a trace table, main or auxiliary in case of RAP,
/// in the first round of the STARK Prove protocol.
pub struct Round1CommitmentData<F>
where
    F: IsField,
    FieldElement<F>: AsBytes + Send + Sync,
{
    /// The result of the interpolation of the columns of the trace table.
    pub(crate) trace_polys: Vec<Polynomial<FieldElement<F>>>,
    /// The Merkle trees constructed to obtain the commitment of the entire trace table.
    pub(crate) lde_trace_merkle_tree: BatchedMerkleTree<F>,
    /// The root of the Merkle tree in `lde_trace_merkle_tree`.
    pub(crate) lde_trace_merkle_root: Commitment,
}

/// A container for the results of the first round of the STARK Prove protocol.
pub struct Round1<A>
where
    A: AIR,
    FieldElement<A::FieldExtension>: AsBytes + Sync + Send,
    FieldElement<A::Field>: AsBytes + Sync + Send,
{
    /// The table of evaluations over the LDE of the main and auxiliary trace tables.
    pub(crate) lde_trace: LDETraceTable<A::Field, A::FieldExtension>,
    /// The intermediate results of the commitment to the main trace table.
    pub(crate) main: Round1CommitmentData<A::Field>,
    /// The intermediate results of the commitment to the auxiliary trace table in case of RAP.
    pub(crate) aux: Option<Round1CommitmentData<A::FieldExtension>>,
    /// The challenges of the RAP round.
    pub(crate) rap_challenges: Vec<FieldElement<A::FieldExtension>>,
}

impl<A> Round1<A>
where
    A: AIR,
    FieldElement<A::FieldExtension>: AsBytes + Sync + Send,
    FieldElement<A::Field>: AsBytes + Sync + Send,
{
    /// Returns the full list of the polynomials interpolating the trace. It includes both
    /// main and auxiliary trace polynomials. The main trace polynomials are casted to
    /// polynomials with coefficients over `Self::FieldExtension`.
    fn all_trace_polys(&self) -> Vec<Polynomial<FieldElement<A::FieldExtension>>> {
        let mut trace_polys: Vec<_> = self
            .main
            .trace_polys
            .clone()
            .into_iter()
            .map(|poly| poly.to_extension())
            .collect();

        if let Some(aux) = &self.aux {
            trace_polys.extend_from_slice(&aux.trace_polys.to_owned())
        }
        trace_polys
    }
}

/// A container for the results of the second round of the STARK Prove protocol.
pub struct Round2<F>
where
    F: IsField,
    FieldElement<F>: AsBytes + Sync + Send,
{
    /// The list of polynomials `H₀, ..., Hₙ` such that `H = ∑ᵢXⁱH(Xⁿ)`, where H is the composition polynomial.
    pub(crate) composition_poly_parts: Vec<Polynomial<FieldElement<F>>>,
    /// Evaluations of the composition polynomial parts over the LDE domain.
    pub(crate) lde_composition_poly_evaluations: Vec<Vec<FieldElement<F>>>,
    /// The Merkle tree built to compute the commitment to the composition polynomial parts.
    pub(crate) composition_poly_merkle_tree: BatchedMerkleTree<F>,
    /// The commitment to the composition polynomial parts.
    pub(crate) composition_poly_root: Commitment,
}

/// A container for the results of the third round of the STARK Prove protocol.
pub struct Round3<F: IsField> {
    /// Evaluations of the trace polynomials, main ans auxiliary, at the out-of-domain challenge.
    trace_ood_evaluations: Table<F>,
    /// Evaluations of the composition polynomial parts at the out-of-domain challenge.
    composition_poly_parts_ood_evaluation: Vec<FieldElement<F>>,
}

/// A container for the results of the fourth round of the STARK Prove protocol.
pub struct Round4<F: IsSubFieldOf<E>, E: IsField> {
    /// The final value resulting from folding the Deep composition polynomial all the way down to a constant value.
    fri_last_value: FieldElement<E>,
    /// The commitments to the fold polynomials of the inner layers of FRI.
    fri_layers_merkle_roots: Vec<Commitment>,
    /// The values and proofs of validity of the evaluations of the trace polynomials and the composition polynomials
    /// parts at the domain values corresponding to the FRI query challenges and their symmetric counterparts.
    deep_poly_openings: DeepPolynomialOpenings<F, E>,
    /// The values and proofs of validity of the evaluations of the fold polynomials of the inner
    /// layers of FRI at the values corresponding to the symmetrics of the FRI query challenges.
    query_list: Vec<FriDecommitment<E>>,
    /// The proof of work nonce.
    nonce: Option<u64>,
}

/// Returns the evaluations of the polynomial `p` over the lde domain defined by the given
/// `blowup_factor`, `domain_size` and `offset`. The number of evaluations returned is `domain_size
/// * blowup_factor`. The domain generator used is the one given by the implementation of `F` as `IsFFTField`.
pub fn evaluate_polynomial_on_lde_domain<F, E>(
    p: &Polynomial<FieldElement<E>>,
    blowup_factor: usize,
    domain_size: usize,
    offset: &FieldElement<F>,
) -> Result<Vec<FieldElement<E>>, FFTError>
where
    F: IsFFTField + IsSubFieldOf<E>,
    E: IsField,
{
    let evaluations = Polynomial::evaluate_offset_fft(p, blowup_factor, Some(domain_size), offset)?;
    let step = evaluations.len() / (domain_size * blowup_factor);
    match step {
        1 => Ok(evaluations),
        _ => Ok(evaluations.into_iter().step_by(step).collect()),
    }
}

/// The functionality of a STARK prover providing methods to run the STARK Prove protocol
/// https://lambdaclass.github.io/lambdaworks/starks/protocol.html
/// The default implementation is complete and is compatible with Stone prover
/// https://github.com/starkware-libs/stone-prover
pub trait IsStarkProver<A: AIR> {
    /// Returns the Merkle tree and the commitment to the vectors `vectors`.
    fn batch_commit_main(
        vectors: &[Vec<FieldElement<A::Field>>],
    ) -> (BatchedMerkleTree<A::Field>, Commitment)
    where
        FieldElement<A::Field>: AsBytes + Sync + Send,
        // FieldElement<A::FieldExtension>: AsBytes + Sync + Send,
        // E: IsSubFieldOf<A::FieldExtension>,
        // A::Field: IsSubFieldOf<E>,
    {
        let tree = BatchedMerkleTree::build(vectors);
        let commitment = tree.root;
        (tree, commitment)
    }

    /// Returns the Merkle tree and the commitment to the vectors `vectors`.
    fn batch_commit_extension(
        vectors: &[Vec<FieldElement<A::FieldExtension>>],
    ) -> (BatchedMerkleTree<A::FieldExtension>, Commitment)
    where
        FieldElement<A::Field>: AsBytes + Sync + Send,
        FieldElement<A::FieldExtension>: AsBytes + Sync + Send,
        // E: IsSubFieldOf<A::FieldExtension>,
        // A::Field: IsSubFieldOf<E>,
    {
        let tree = BatchedMerkleTree::build(vectors);
        let commitment = tree.root;
        (tree, commitment)
    }

    /// Given a `TraceTable`, this method interpolates its columns, computes the commitment to the
    /// table and appends it to the transcript.
    /// Output: a touple of length 4 with the following:
    /// • The polynomials interpolating the columns of `trace`.
    /// • The evaluations of the above polynomials over the domain `domain`.
    /// • The Merkle tree of evaluations of the above polynomials over the domain `domain`.
    /// • The roots of the above Merkle trees.
    #[allow(clippy::type_complexity)]
    fn interpolate_and_commit_main(
        trace: &TraceTable<A::Field, A::FieldExtension>,
        domain: &Domain<A::Field>,
        transcript: &mut impl IsStarkTranscript<A::FieldExtension>,
    ) -> (
        Vec<Polynomial<FieldElement<A::Field>>>,
        Vec<Vec<FieldElement<A::Field>>>,
        BatchedMerkleTree<A::Field>,
        Commitment,
    )
    where
        FieldElement<A::Field>: AsBytes + Send + Sync,
        // FieldElement<E>: AsBytes + Send + Sync,
        FieldElement<A::FieldExtension>: AsBytes + Send + Sync,
        A::Field: IsSubFieldOf<A::FieldExtension>,
    {
        // Interpolate columns of `trace`.
        let trace_polys = trace.compute_trace_polys_main::<A::Field>();

        // Evaluate those polynomials t_j on the large domain D_LDE.
        let lde_trace_evaluations =
            Self::compute_lde_trace_evaluations::<A::Field>(&trace_polys, domain);

        let mut lde_trace_permuted = lde_trace_evaluations.clone();
        for col in lde_trace_permuted.iter_mut() {
            in_place_bit_reverse_permute(col);
        }

        // Compute commitment.
        let lde_trace_permuted_rows = columns2rows(lde_trace_permuted);
        let (lde_trace_merkle_tree, lde_trace_merkle_root) =
            Self::batch_commit_main(&lde_trace_permuted_rows);

        // >>>> Send commitment.
        transcript.append_bytes(&lde_trace_merkle_root);

        (
            trace_polys,
            lde_trace_evaluations,
            lde_trace_merkle_tree,
            lde_trace_merkle_root,
        )
    }

    /// Given a `TraceTable`, this method interpolates its columns, computes the commitment to the
    /// table and appends it to the transcript.
    /// Output: a touple of length 4 with the following:
    /// • The polynomials interpolating the columns of `trace`.
    /// • The evaluations of the above polynomials over the domain `domain`.
    /// • The Merkle tree of evaluations of the above polynomials over the domain `domain`.
    /// • The roots of the above Merkle trees.
    #[allow(clippy::type_complexity)]
    fn interpolate_and_commit_aux(
        trace: &TraceTable<A::Field, A::FieldExtension>,
        domain: &Domain<A::Field>,
        transcript: &mut impl IsStarkTranscript<A::FieldExtension>,
    ) -> (
        Vec<Polynomial<FieldElement<A::FieldExtension>>>,
        Vec<Vec<FieldElement<A::FieldExtension>>>,
        BatchedMerkleTree<A::FieldExtension>,
        Commitment,
    )
    where
        FieldElement<A::Field>: AsBytes + Send + Sync,
        // FieldElement<E>: AsBytes + Send + Sync,
        FieldElement<A::FieldExtension>: AsBytes + Send + Sync,
        // E: IsSubFieldOf<A::FieldExtension> + IsFFTField,
        A::Field: IsSubFieldOf<A::FieldExtension> + IsFFTField,
    {
        // Interpolate columns of `trace`.
        let trace_polys = trace.compute_trace_polys_aux::<A::Field>();

        // Evaluate those polynomials t_j on the large domain D_LDE.
        let lde_trace_evaluations = Self::compute_lde_trace_evaluations(&trace_polys, domain);

        let mut lde_trace_permuted = lde_trace_evaluations.clone();
        for col in lde_trace_permuted.iter_mut() {
            in_place_bit_reverse_permute(col);
        }

        // Compute commitment.
        let lde_trace_permuted_rows = columns2rows(lde_trace_permuted);
        let (lde_trace_merkle_tree, lde_trace_merkle_root) =
            Self::batch_commit_extension(&lde_trace_permuted_rows);

        // >>>> Send commitment.
        transcript.append_bytes(&lde_trace_merkle_root);

        (
            trace_polys,
            lde_trace_evaluations,
            lde_trace_merkle_tree,
            lde_trace_merkle_root,
        )
    }

    /// Evaluate polynomials `trace_polys` over the domain `domain`.
    /// The i-th entry of the returned vector contains the evaluations of the i-th polynomial in `trace_polys`.
    fn compute_lde_trace_evaluations<E>(
        trace_polys: &[Polynomial<FieldElement<E>>],
        domain: &Domain<A::Field>,
    ) -> Vec<Vec<FieldElement<E>>>
    where
        FieldElement<E>: Send + Sync,
        FieldElement<A::Field>: Send + Sync,
        E: IsSubFieldOf<A::FieldExtension>,
        A::Field: IsSubFieldOf<E>,
        // F: IsFFTField,
    {
        #[cfg(not(feature = "parallel"))]
        let trace_polys_iter = trace_polys.iter();
        #[cfg(feature = "parallel")]
        let trace_polys_iter = trace_polys.par_iter();

        trace_polys_iter
            .map(|poly| {
                evaluate_polynomial_on_lde_domain(
                    poly,
                    domain.blowup_factor,
                    domain.interpolation_domain_size,
                    &domain.coset_offset,
                )
            })
            .collect::<Result<Vec<Vec<FieldElement<E>>>, FFTError>>()
            .unwrap()
    }

    /// Returns the result of the first round of the STARK Prove protocol.
    fn round_1_randomized_air_with_preprocessing(
        air: &A,
        trace: &mut TraceTable<A::Field, A::FieldExtension>,
        domain: &Domain<A::Field>,
        transcript: &mut impl IsStarkTranscript<A::FieldExtension>,
    ) -> Result<Round1<A>, ProvingError>
    where
        FieldElement<A::Field>: AsBytes + Send + Sync,
        A::FieldExtension: IsFFTField,
        FieldElement<A::FieldExtension>: AsBytes + Send + Sync,
    {
        let (trace_polys, evaluations, main_merkle_tree, main_merkle_root) =
            // Self::interpolate_and_commit_main::<A::Field>(trace, domain, transcript);
            Self::interpolate_and_commit_main(trace, domain, transcript);

        let main = Round1CommitmentData::<A::Field> {
            trace_polys,
            lde_trace_merkle_tree: main_merkle_tree,
            lde_trace_merkle_root: main_merkle_root,
        };

        let rap_challenges = air.build_rap_challenges(transcript);
        let (aux, aux_evaluations) = if air.has_trace_interaction() {
            air.build_auxiliary_trace(trace, &rap_challenges);
            let (aux_trace_polys, aux_trace_polys_evaluations, aux_merkle_tree, aux_merkle_root) =
                Self::interpolate_and_commit_aux(trace, domain, transcript);
            let aux_evaluations = aux_trace_polys_evaluations;
            let aux = Some(Round1CommitmentData::<A::FieldExtension> {
                trace_polys: aux_trace_polys,
                lde_trace_merkle_tree: aux_merkle_tree,
                lde_trace_merkle_root: aux_merkle_root,
            });
            (aux, aux_evaluations)
        } else {
            (None, Vec::new())
        };

        let lde_trace = LDETraceTable::from_columns(
            evaluations,
            aux_evaluations,
            A::STEP_SIZE,
            domain.blowup_factor,
        );

        Ok(Round1 {
            lde_trace,
            main,
            aux,
            rap_challenges,
        })
    }

    /// Returns the Merkle tree and the commitment to the evaluations of the parts of the
    /// composition polynomial.
    fn commit_composition_polynomial(
        lde_composition_poly_parts_evaluations: &[Vec<FieldElement<A::FieldExtension>>],
    ) -> (BatchedMerkleTree<A::FieldExtension>, Commitment)
    where
        FieldElement<A::Field>: AsBytes + Sync + Send,
        FieldElement<A::FieldExtension>: AsBytes + Sync + Send,
    {
        // TODO: Remove clones
        let mut lde_composition_poly_evaluations = Vec::new();
        for i in 0..lde_composition_poly_parts_evaluations[0].len() {
            let mut row = Vec::new();
            for evaluation in lde_composition_poly_parts_evaluations.iter() {
                row.push(evaluation[i].clone());
            }
            lde_composition_poly_evaluations.push(row);
        }

        in_place_bit_reverse_permute(&mut lde_composition_poly_evaluations);

        let mut lde_composition_poly_evaluations_merged = Vec::new();
        for chunk in lde_composition_poly_evaluations.chunks(2) {
            let (mut chunk0, chunk1) = (chunk[0].clone(), &chunk[1]);
            chunk0.extend_from_slice(chunk1);
            lde_composition_poly_evaluations_merged.push(chunk0);
        }

        Self::batch_commit_extension(&lde_composition_poly_evaluations_merged)
    }

    /// Returns the result of the second round of the STARK Prove protocol.
    fn round_2_compute_composition_polynomial(
        air: &A,
        domain: &Domain<A::Field>,
        round_1_result: &Round1<A>,
        transition_coefficients: &[FieldElement<A::FieldExtension>],
        boundary_coefficients: &[FieldElement<A::FieldExtension>],
    ) -> Round2<A::FieldExtension>
    where
        A: Send + Sync,
        FieldElement<A::Field>: AsBytes + Send + Sync,
        FieldElement<A::FieldExtension>: AsBytes + Send + Sync,
    {
        // Compute the evaluations of the composition polynomial on the LDE domain.
        let evaluator = ConstraintEvaluator::new(air, &round_1_result.rap_challenges);
        let constraint_evaluations = evaluator.evaluate(
            air,
            &round_1_result.lde_trace,
            domain,
            transition_coefficients,
            boundary_coefficients,
            &round_1_result.rap_challenges,
        );

        // Get coefficients of the composition poly H
        let composition_poly =
            Polynomial::interpolate_offset_fft(&constraint_evaluations, &domain.coset_offset)
                .unwrap();

        let number_of_parts = air.composition_poly_degree_bound() / air.trace_length();
        let composition_poly_parts = composition_poly.break_in_parts(number_of_parts);

        let lde_composition_poly_parts_evaluations: Vec<_> = composition_poly_parts
            .iter()
            .map(|part| {
                evaluate_polynomial_on_lde_domain(
                    part,
                    domain.blowup_factor,
                    domain.interpolation_domain_size,
                    &domain.coset_offset,
                )
                .unwrap()
            })
            .collect();

        let (composition_poly_merkle_tree, composition_poly_root) =
            Self::commit_composition_polynomial(&lde_composition_poly_parts_evaluations);

        Round2 {
            lde_composition_poly_evaluations: lde_composition_poly_parts_evaluations,
            composition_poly_parts,
            composition_poly_merkle_tree,
            composition_poly_root,
        }
    }

    /// Returns the result of the third round of the STARK Prove protocol.
    fn round_3_evaluate_polynomials_in_out_of_domain_element(
        air: &A,
        domain: &Domain<A::Field>,
        round_1_result: &Round1<A>,
        round_2_result: &Round2<A::FieldExtension>,
        z: &FieldElement<A::FieldExtension>,
    ) -> Round3<A::FieldExtension>
    where
        FieldElement<A::Field>: AsBytes + Sync + Send,
        FieldElement<A::FieldExtension>: AsBytes + Sync + Send,
    {
        let z_power = z.pow(round_2_result.composition_poly_parts.len());

        // Evaluate H_i in z^N for all i, where N is the number of parts the composition poly was
        // broken into.
        let composition_poly_parts_ood_evaluation: Vec<_> = round_2_result
            .composition_poly_parts
            .iter()
            .map(|part| part.evaluate(&z_power))
            .collect();

        // Returns the Out of Domain Frame for the given trace polynomials, out of domain evaluation point (called `z` in the literature),
        // frame offsets given by the AIR and primitive root used for interpolating the trace polynomials.
        // An out of domain frame is nothing more than the evaluation of the trace polynomials in the points required by the
        // verifier to check the consistency between the trace and the composition polynomial.
        //
        // In the fibonacci example, the ood frame is simply the evaluations `[t(z), t(z * g), t(z * g^2)]`, where `t` is the trace
        // polynomial and `g` is the primitive root of unity used when interpolating `t`.
        let trace_ood_evaluations =
            crate::trace::get_trace_evaluations::<A::Field, A::FieldExtension>(
                &round_1_result.main.trace_polys,
                round_1_result
                    .aux
                    .as_ref()
                    .map(|aux| &aux.trace_polys)
                    .unwrap_or(&vec![]),
                z,
                &air.context().transition_offsets,
                &domain.trace_primitive_root,
                A::STEP_SIZE,
            );

        Round3 {
            trace_ood_evaluations,
            composition_poly_parts_ood_evaluation,
        }
    }

    /// Returns the result of the fourth round of the STARK Prove protocol.
    fn round_4_compute_and_run_fri_on_the_deep_composition_polynomial(
        air: &A,
        domain: &Domain<A::Field>,
        round_1_result: &Round1<A>,
        round_2_result: &Round2<A::FieldExtension>,
        round_3_result: &Round3<A::FieldExtension>,
        z: &FieldElement<A::FieldExtension>,
        transcript: &mut impl IsStarkTranscript<A::FieldExtension>,
    ) -> Round4<A::Field, A::FieldExtension>
    where
        FieldElement<A::Field>: AsBytes + Send + Sync,
        FieldElement<A::FieldExtension>: AsBytes + Send + Sync,
    {
        let coset_offset_u64 = air.context().proof_options.coset_offset;
        let coset_offset = FieldElement::<A::Field>::from(coset_offset_u64);

        let gamma = transcript.sample_field_element();

        let n_terms_composition_poly = round_2_result.lde_composition_poly_evaluations.len();
        let num_terms_trace =
            air.context().transition_offsets.len() * A::STEP_SIZE * air.context().trace_columns;

        // <<<< Receive challenges: 𝛾, 𝛾'
        let mut deep_composition_coefficients: Vec<_> =
            core::iter::successors(Some(FieldElement::one()), |x| Some(x * &gamma))
                .take(n_terms_composition_poly + num_terms_trace)
                .collect();

        // let trace_poly_coeffients: Vec<_> = deep_composition_coefficients
        //     .drain(..num_terms_trace)
        //     .collect();

        let trace_term_coeffs: Vec<_> = deep_composition_coefficients
            .drain(..num_terms_trace)
            .collect::<Vec<_>>()
            .chunks(air.context().transition_offsets.len() * A::STEP_SIZE)
            .map(|chunk| chunk.to_vec())
            .collect();

        // <<<< Receive challenges: 𝛾ⱼ, 𝛾ⱼ'
        let gammas = deep_composition_coefficients;

        // Compute p₀ (deep composition polynomial)
        let deep_composition_poly = Self::compute_deep_composition_poly(
            &round_1_result.all_trace_polys(),
            round_2_result,
            round_3_result,
            z,
            &domain.trace_primitive_root,
            &gammas,
            &trace_term_coeffs,
        );

        let domain_size = domain.lde_roots_of_unity_coset.len();

        // FRI commit and query phases
        let (fri_last_value, fri_layers) = fri::commit_phase::<A::Field, A::FieldExtension>(
            domain.root_order as usize,
            deep_composition_poly,
            transcript,
            &coset_offset,
            domain_size,
        );

        // grinding: generate nonce and append it to the transcript
        let security_bits = air.context().proof_options.grinding_factor;
        let mut nonce = None;
        if security_bits > 0 {
            let nonce_value = grinding::generate_nonce(&transcript.state(), security_bits)
                .expect("nonce not found");
            transcript.append_bytes(&nonce_value.to_be_bytes());
            nonce = Some(nonce_value);
        }

        let number_of_queries = air.options().fri_number_of_queries;
        let iotas = Self::sample_query_indexes(number_of_queries, domain, transcript);

        let query_list = fri::query_phase(&fri_layers, &iotas);

        let fri_layers_merkle_roots: Vec<_> = fri_layers
            .iter()
            .map(|layer| layer.merkle_tree.root)
            .collect();

        let deep_poly_openings =
            Self::open_deep_composition_poly(domain, round_1_result, round_2_result, &iotas);

        Round4 {
            fri_last_value,
            fri_layers_merkle_roots,
            deep_poly_openings,
            query_list,
            nonce,
        }
    }

    fn sample_query_indexes(
        number_of_queries: usize,
        domain: &Domain<A::Field>,
        transcript: &mut impl IsStarkTranscript<A::FieldExtension>,
    ) -> Vec<usize> {
        let domain_size = domain.lde_roots_of_unity_coset.len() as u64;
        (0..number_of_queries)
            .map(|_| (transcript.sample_u64(domain_size >> 1)) as usize)
            .collect::<Vec<usize>>()
    }

    /// Returns the DEEP composition polynomial that the prover then commits to using
    /// FRI. This polynomial is a linear combination of the trace polynomial and the
    /// composition polynomial, with coefficients sampled by the verifier (i.e. using Fiat-Shamir).
    #[allow(clippy::too_many_arguments)]
    fn compute_deep_composition_poly(
        trace_polys: &[Polynomial<FieldElement<A::FieldExtension>>],
        round_2_result: &Round2<A::FieldExtension>,
        round_3_result: &Round3<A::FieldExtension>,
        z: &FieldElement<A::FieldExtension>,
        primitive_root: &FieldElement<A::Field>,
        composition_poly_gammas: &[FieldElement<A::FieldExtension>],
        trace_terms_gammas: &[Vec<FieldElement<A::FieldExtension>>],
    ) -> Polynomial<FieldElement<A::FieldExtension>>
    where
        FieldElement<A::Field>: AsBytes + Send + Sync,
        FieldElement<A::FieldExtension>: AsBytes + Send + Sync,
    {
        let z_power = z.pow(round_2_result.composition_poly_parts.len());

        // ∑ᵢ 𝛾ᵢ ( Hᵢ − Hᵢ(z^N) ) / ( X − z^N )
        let mut h_terms = Polynomial::zero();
        for (i, part) in round_2_result.composition_poly_parts.iter().enumerate() {
            // h_i_eval is the evaluation of the i-th part of the composition polynomial at z^N,
            // where N is the number of parts of the composition polynomial.
            let h_i_eval = &round_3_result.composition_poly_parts_ood_evaluation[i];
            let h_i_term = &composition_poly_gammas[i] * (part - h_i_eval);
            h_terms = h_terms + h_i_term;
        }
        assert_eq!(h_terms.evaluate(&z_power), FieldElement::zero());
        h_terms.ruffini_division_inplace(&z_power);

        // Get trace evaluations needed for the trace terms of the deep composition polynomial
        let trace_frame_evaluations = &round_3_result.trace_ood_evaluations;

        // Compute the sum of all the trace terms of the deep composition polynomial.
        // There is one term for every trace polynomial and for every row in the frame.
        // ∑ ⱼₖ [ 𝛾ₖ ( tⱼ − tⱼ(z) ) / ( X − zgᵏ )]

        let trace_evaluations_columns = &trace_frame_evaluations.columns();

        #[cfg(feature = "parallel")]
        let trace_terms = trace_polys
            .par_iter()
            .enumerate()
            .fold(Polynomial::zero, |trace_terms, (i, t_j)| {
                let gammas_i = &trace_terms_gammas[i];
                let trace_evaluations_i = &trace_evaluations_columns[i];
                Self::compute_trace_term(
                    &trace_terms,
                    t_j,
                    gammas_i,
                    trace_evaluations_i,
                    (z, primitive_root),
                )
            })
            .reduce(Polynomial::zero, |a, b| a + b);

        #[cfg(not(feature = "parallel"))]
        let trace_terms =
            trace_polys
                .iter()
                .enumerate()
                .fold(Polynomial::zero(), |trace_terms, (i, t_j)| {
                    let gammas_i = &trace_terms_gammas[i];
                    let trace_evaluations_i = &trace_evaluations_columns[i];
                    Self::compute_trace_term(
                        &trace_terms,
                        t_j,
                        gammas_i,
                        trace_evaluations_i,
                        (z, primitive_root),
                    )
                });

        h_terms + trace_terms
    }

    // FIXME: FIX THIS DOCS!
    /// Adds to `accumulator` the term corresponding to the trace polynomial `t_j` of the Deep
    /// composition polynomial. That is, returns `accumulator + \sum_i \gamma_i \frac{ t_j - t_j(zg^i) }{ X - zg^i }`,
    /// where `i` ranges from `T * j` to `T * j + T - 1`, where `T` is the number of offsets in every frame.
    fn compute_trace_term(
        accumulator: &Polynomial<FieldElement<A::FieldExtension>>,
        trace_term_poly: &Polynomial<FieldElement<A::FieldExtension>>,
        trace_terms_gammas: &[FieldElement<A::FieldExtension>],
        trace_frame_evaluations: &[FieldElement<A::FieldExtension>],
        (z, primitive_root): (&FieldElement<A::FieldExtension>, &FieldElement<A::Field>),
    ) -> Polynomial<FieldElement<A::FieldExtension>>
    where
        FieldElement<A::Field>: AsBytes + Send + Sync,
        FieldElement<A::FieldExtension>: AsBytes + Send + Sync,
    {
        let trace_int = trace_frame_evaluations
            .iter()
            .enumerate()
            .zip(trace_terms_gammas)
            .fold(
                Polynomial::zero(),
                |trace_agg, ((offset, trace_term_poly_evaluation), trace_gamma)| {
                    // @@@ this can be pre-computed
                    let z_shifted = primitive_root.pow(offset) * z;
                    let mut poly = trace_term_poly - trace_term_poly_evaluation;
                    poly.ruffini_division_inplace(&z_shifted);
                    trace_agg + poly * trace_gamma
                },
            );

        accumulator + trace_int
    }

    /// Computes values and validity proofs of the evaluations of the composition polynomial parts
    /// at the domain value corresponding to the FRI query challenge `index` and its symmetric
    /// element.
    fn open_composition_poly(
        composition_poly_merkle_tree: &BatchedMerkleTree<A::FieldExtension>,
        lde_composition_poly_evaluations: &[Vec<FieldElement<A::FieldExtension>>],
        index: usize,
    ) -> PolynomialOpenings<A::FieldExtension>
    where
        FieldElement<A::Field>: AsBytes + Sync + Send,
        FieldElement<A::FieldExtension>: AsBytes + Sync + Send,
    {
        let proof = composition_poly_merkle_tree
            .get_proof_by_pos(index)
            .unwrap();

        let lde_composition_poly_parts_evaluation: Vec<_> = lde_composition_poly_evaluations
            .iter()
            .flat_map(|part| {
                vec![
                    part[reverse_index(index * 2, part.len() as u64)].clone(),
                    part[reverse_index(index * 2 + 1, part.len() as u64)].clone(),
                ]
            })
            .collect();

        PolynomialOpenings {
            proof: proof.clone(),
            proof_sym: proof,
            evaluations: lde_composition_poly_parts_evaluation
                .clone()
                .into_iter()
                .step_by(2)
                .collect(),
            evaluations_sym: lde_composition_poly_parts_evaluation
                .into_iter()
                .skip(1)
                .step_by(2)
                .collect(),
        }
    }

    /// Computes values and validity proofs of the evaluations of the trace polynomials
    /// at the domain value corresponding to the FRI query challenge `index` and its symmetric
    /// element.
    fn open_trace_polys<E>(
        domain: &Domain<A::Field>,
        tree: &BatchedMerkleTree<E>,
        lde_trace: &Table<E>,
        challenge: usize,
    ) -> PolynomialOpenings<E>
    where
        FieldElement<A::Field>: AsBytes + Sync + Send,
        FieldElement<E>: AsBytes + Sync + Send,
        A::Field: IsSubFieldOf<E>,
        E: IsField,
    {
        let domain_size = domain.lde_roots_of_unity_coset.len();

        let index = challenge * 2;
        let index_sym = challenge * 2 + 1;
        PolynomialOpenings {
            proof: tree.get_proof_by_pos(index).unwrap(),
            proof_sym: tree.get_proof_by_pos(index_sym).unwrap(),
            evaluations: lde_trace
                .get_row(reverse_index(index, domain_size as u64))
                .to_vec(),
            evaluations_sym: lde_trace
                .get_row(reverse_index(index_sym, domain_size as u64))
                .to_vec(),
        }
    }

    /// Open the deep composition polynomial on a list of indexes and their symmetric elements.
    fn open_deep_composition_poly(
        domain: &Domain<A::Field>,
        round_1_result: &Round1<A>,
        round_2_result: &Round2<A::FieldExtension>,
        indexes_to_open: &[usize],
    ) -> DeepPolynomialOpenings<A::Field, A::FieldExtension>
    where
        FieldElement<A::Field>: AsBytes + Send + Sync,
        FieldElement<A::FieldExtension>: AsBytes + Send + Sync,
    {
        let mut openings = Vec::new();

        for index in indexes_to_open.iter() {
            let main_trace_opening = Self::open_trace_polys::<A::Field>(
                domain,
                &round_1_result.main.lde_trace_merkle_tree,
                &round_1_result.lde_trace.main_table,
                *index,
            );

            let composition_openings = Self::open_composition_poly(
                &round_2_result.composition_poly_merkle_tree,
                &round_2_result.lde_composition_poly_evaluations,
                *index,
            );

            let aux_trace_polys = round_1_result.aux.as_ref().map(|aux| {
                Self::open_trace_polys::<A::FieldExtension>(
                    domain,
                    &aux.lde_trace_merkle_tree,
                    &round_1_result.lde_trace.aux_table,
                    *index,
                )
            });

            openings.push(DeepPolynomialOpening {
                composition_poly: composition_openings,
                main_trace_polys: main_trace_opening,
                aux_trace_polys,
            });
        }

        openings
    }

    // FIXME remove unwrap() calls and return errors
    /// Generates a STARK proof for the trace `main_trace` with public inputs `pub_inputs`.
    /// Warning: the transcript must be safely initializated before passing it to this method.
    fn prove(
        trace: &mut TraceTable<A::Field, A::FieldExtension>,
        pub_inputs: &A::PublicInputs,
        proof_options: &ProofOptions,
        mut transcript: impl IsStarkTranscript<A::FieldExtension>,
    ) -> Result<StarkProof<A::Field, A::FieldExtension>, ProvingError>
    where
        A: Send + Sync,
        FieldElement<A::Field>: AsBytes + Send + Sync,
        A::FieldExtension: IsFFTField,
        FieldElement<A::FieldExtension>: AsBytes + Send + Sync,
    {
        info!("Started proof generation...");
        #[cfg(feature = "instruments")]
        println!("- Started round 0: Air Initialization");
        #[cfg(feature = "instruments")]
        let timer0 = Instant::now();

        let air = A::new(trace.num_rows(), pub_inputs, proof_options);
        let domain = Domain::new(&air);

        #[cfg(feature = "instruments")]
        let elapsed0 = timer0.elapsed();
        #[cfg(feature = "instruments")]
        println!("  Time spent: {:?}", elapsed0);

        // ===================================
        // ==========|   Round 1   |==========
        // ===================================

        #[cfg(feature = "instruments")]
        println!("- Started round 1: RAP");
        #[cfg(feature = "instruments")]
        let timer1 = Instant::now();

        let round_1_result =
            Self::round_1_randomized_air_with_preprocessing(&air, trace, &domain, &mut transcript)?;

        #[cfg(debug_assertions)]
        validate_trace(
            &air,
            &round_1_result.main.trace_polys,
            round_1_result
                .aux
                .as_ref()
                .map(|a| &a.trace_polys)
                .unwrap_or(&vec![]),
            &domain,
            &round_1_result.rap_challenges,
        );

        #[cfg(feature = "instruments")]
        let elapsed1 = timer1.elapsed();
        #[cfg(feature = "instruments")]
        println!("  Time spent: {:?}", elapsed1);

        // ===================================
        // ==========|   Round 2   |==========
        // ===================================

        #[cfg(feature = "instruments")]
        println!("- Started round 2: Compute composition polynomial");
        #[cfg(feature = "instruments")]
        let timer2 = Instant::now();

        // <<<< Receive challenge: 𝛽
        let beta = transcript.sample_field_element();
        let num_boundary_constraints = air
            .boundary_constraints(&round_1_result.rap_challenges)
            .constraints
            .len();

        let num_transition_constraints = air.context().num_transition_constraints;

        let mut coefficients: Vec<_> =
            core::iter::successors(Some(FieldElement::one()), |x| Some(x * &beta))
                .take(num_boundary_constraints + num_transition_constraints)
                .collect();

        let transition_coefficients: Vec<_> =
            coefficients.drain(..num_transition_constraints).collect();
        let boundary_coefficients = coefficients;

        let round_2_result = Self::round_2_compute_composition_polynomial(
            &air,
            &domain,
            &round_1_result,
            &transition_coefficients,
            &boundary_coefficients,
        );

        // >>>> Send commitments: [H₁], [H₂]
        transcript.append_bytes(&round_2_result.composition_poly_root);

        #[cfg(feature = "instruments")]
        let elapsed2 = timer2.elapsed();
        #[cfg(feature = "instruments")]
        println!("  Time spent: {:?}", elapsed2);

        // ===================================
        // ==========|   Round 3   |==========
        // ===================================

        #[cfg(feature = "instruments")]
        println!("- Started round 3: Evaluate polynomial in out of domain elements");
        #[cfg(feature = "instruments")]
        let timer3 = Instant::now();

        // <<<< Receive challenge: z
        let z = transcript.sample_z_ood(
            &domain.lde_roots_of_unity_coset,
            &domain.trace_roots_of_unity,
        );

        let round_3_result = Self::round_3_evaluate_polynomials_in_out_of_domain_element(
            &air,
            &domain,
            &round_1_result,
            &round_2_result,
            &z,
        );

        // >>>> Send values: tⱼ(zgᵏ)
        let trace_ood_evaluations_columns = round_3_result.trace_ood_evaluations.columns();
        for col in trace_ood_evaluations_columns.iter() {
            for elem in col.iter() {
                transcript.append_field_element(elem);
            }
        }

        // >>>> Send values: Hᵢ(z^N)
        for element in round_3_result.composition_poly_parts_ood_evaluation.iter() {
            transcript.append_field_element(element);
        }

        #[cfg(feature = "instruments")]
        let elapsed3 = timer3.elapsed();
        #[cfg(feature = "instruments")]
        println!("  Time spent: {:?}", elapsed3);

        // ===================================
        // ==========|   Round 4   |==========
        // ===================================

        #[cfg(feature = "instruments")]
        println!("- Started round 4: FRI");
        #[cfg(feature = "instruments")]
        let timer4 = Instant::now();

        // Part of this round is running FRI, which is an interactive
        // protocol on its own. Therefore we pass it the transcript
        // to simulate the interactions with the verifier.
        let round_4_result = Self::round_4_compute_and_run_fri_on_the_deep_composition_polynomial(
            &air,
            &domain,
            &round_1_result,
            &round_2_result,
            &round_3_result,
            &z,
            &mut transcript,
        );

        #[cfg(feature = "instruments")]
        let elapsed4 = timer4.elapsed();
        #[cfg(feature = "instruments")]
        println!("  Time spent: {:?}", elapsed4);

        #[cfg(feature = "instruments")]
        {
            let total_time = elapsed1 + elapsed2 + elapsed3 + elapsed4;
            println!(
                " Fraction of proving time per round: {:.4} {:.4} {:.4} {:.4} {:.4}",
                elapsed0.as_nanos() as f64 / total_time.as_nanos() as f64,
                elapsed1.as_nanos() as f64 / total_time.as_nanos() as f64,
                elapsed2.as_nanos() as f64 / total_time.as_nanos() as f64,
                elapsed3.as_nanos() as f64 / total_time.as_nanos() as f64,
                elapsed4.as_nanos() as f64 / total_time.as_nanos() as f64
            );
        }

        info!("End proof generation");

        Ok(StarkProof::<A::Field, A::FieldExtension> {
            // [t]
            lde_trace_main_merkle_root: round_1_result.main.lde_trace_merkle_root,
            // [t]
            lde_trace_aux_merkle_root: round_1_result.aux.map(|x| x.lde_trace_merkle_root),
            // tⱼ(zgᵏ)
            trace_ood_evaluations: round_3_result.trace_ood_evaluations,
            // [H₁] and [H₂]
            composition_poly_root: round_2_result.composition_poly_root,
            // Hᵢ(z^N)
            composition_poly_parts_ood_evaluation: round_3_result
                .composition_poly_parts_ood_evaluation,
            // [pₖ]
            fri_layers_merkle_roots: round_4_result.fri_layers_merkle_roots,
            // pₙ
            fri_last_value: round_4_result.fri_last_value,
            // Open(p₀(D₀), 𝜐ₛ), Open(pₖ(Dₖ), −𝜐ₛ^(2ᵏ))
            query_list: round_4_result.query_list,
            // Open(H₁(D_LDE, 𝜐₀), Open(H₂(D_LDE, 𝜐₀), Open(tⱼ(D_LDE), 𝜐₀)
            // Open(H₁(D_LDE, -𝜐ᵢ), Open(H₂(D_LDE, -𝜐ᵢ), Open(tⱼ(D_LDE), -𝜐ᵢ)
            deep_poly_openings: round_4_result.deep_poly_openings,
            // nonce obtained from grinding
            nonce: round_4_result.nonce,

            trace_length: air.trace_length(),
        })
    }
}

#[cfg(test)]
mod tests {
    use std::num::ParseIntError;

    fn decode_hex(s: &str) -> Result<Vec<u8>, ParseIntError> {
        (0..s.len())
            .step_by(2)
            .map(|i| u8::from_str_radix(&s[i..i + 2], 16))
            .collect()
    }

    use crate::{
        examples::{
            fibonacci_2_cols_shifted::{self, Fibonacci2ColsShifted},
            simple_fibonacci::{self, FibonacciPublicInputs},
        },
        proof::options::ProofOptions,
        transcript::StoneProverTranscript,
        verifier::{Challenges, IsStarkVerifier, Verifier},
        Felt252,
    };

    use super::*;
    use lambdaworks_math::{
        field::{
            element::FieldElement, fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
            traits::IsFFTField,
        },
        polynomial::Polynomial,
    };

    #[test]
    fn test_domain_constructor() {
        let pub_inputs = FibonacciPublicInputs {
            a0: Felt252::one(),
            a1: Felt252::one(),
        };
        let trace = simple_fibonacci::fibonacci_trace([Felt252::from(1), Felt252::from(1)], 8);
        let trace_length = trace.num_rows();
        let coset_offset = 3;
        let blowup_factor: usize = 2;
        let grinding_factor = 20;

        let proof_options = ProofOptions {
            blowup_factor: blowup_factor as u8,
            fri_number_of_queries: 1,
            coset_offset,
            grinding_factor,
        };

        let domain = Domain::new(&simple_fibonacci::FibonacciAIR::new(
            trace_length,
            &pub_inputs,
            &proof_options,
        ));
        assert_eq!(domain.blowup_factor, 2);
        assert_eq!(domain.interpolation_domain_size, trace_length);
        assert_eq!(domain.root_order, trace_length.trailing_zeros());
        assert_eq!(domain.coset_offset, FieldElement::from(coset_offset));

        let primitive_root = Stark252PrimeField::get_primitive_root_of_unity(
            (trace_length * blowup_factor).trailing_zeros() as u64,
        )
        .unwrap();

        assert_eq!(
            domain.trace_primitive_root,
            primitive_root.pow(blowup_factor)
        );
        for i in 0..(trace_length * blowup_factor) {
            assert_eq!(
                domain.lde_roots_of_unity_coset[i],
                primitive_root.pow(i) * FieldElement::from(coset_offset)
            );
        }
    }

    #[test]
    fn test_evaluate_polynomial_on_lde_domain_on_trace_polys() {
        let trace = simple_fibonacci::fibonacci_trace([Felt252::from(1), Felt252::from(1)], 8);

        let trace_length = trace.num_rows();

        let trace_polys = trace.compute_trace_polys_main::<Stark252PrimeField>();
        let coset_offset = Felt252::from(3);
        let blowup_factor: usize = 2;
        let domain_size = 8;

        let primitive_root = Stark252PrimeField::get_primitive_root_of_unity(
            (trace_length * blowup_factor).trailing_zeros() as u64,
        )
        .unwrap();

        for poly in trace_polys.iter() {
            let lde_evaluation =
                evaluate_polynomial_on_lde_domain(poly, blowup_factor, domain_size, &coset_offset)
                    .unwrap();
            assert_eq!(lde_evaluation.len(), trace_length * blowup_factor);
            for (i, evaluation) in lde_evaluation.iter().enumerate() {
                assert_eq!(
                    *evaluation,
                    poly.evaluate(&(coset_offset * primitive_root.pow(i)))
                );
            }
        }
    }

    #[test]
    fn test_evaluate_polynomial_on_lde_domain_edge_case() {
        let poly = Polynomial::new_monomial(Felt252::one(), 8);
        let blowup_factor: usize = 4;
        let domain_size: usize = 8;
        let offset = Felt252::from(3);
        let evaluations =
            evaluate_polynomial_on_lde_domain(&poly, blowup_factor, domain_size, &offset).unwrap();
        assert_eq!(evaluations.len(), domain_size * blowup_factor);

        let primitive_root: Felt252 = Stark252PrimeField::get_primitive_root_of_unity(
            (domain_size * blowup_factor).trailing_zeros() as u64,
        )
        .unwrap();
        for (i, eval) in evaluations.iter().enumerate() {
            assert_eq!(*eval, poly.evaluate(&(offset * primitive_root.pow(i))));
        }
    }

    fn proof_parts_stone_compatibility_case_1() -> (
        StarkProof<Stark252PrimeField, Stark252PrimeField>,
        fibonacci_2_cols_shifted::PublicInputs<Stark252PrimeField>,
        ProofOptions,
        [u8; 4],
    ) {
        let mut trace = fibonacci_2_cols_shifted::compute_trace(FieldElement::one(), 4);

        let claimed_index = 3;
        let col = 0;
        let claimed_value = *trace.get_main(claimed_index, col);
        let mut proof_options = ProofOptions::default_test_options();
        proof_options.blowup_factor = 4;
        proof_options.coset_offset = 3;
        proof_options.grinding_factor = 0;
        proof_options.fri_number_of_queries = 1;

        let pub_inputs = fibonacci_2_cols_shifted::PublicInputs {
            claimed_value,
            claimed_index,
        };

        let transcript_init_seed = [0xca, 0xfe, 0xca, 0xfe];

        let proof = Prover::<Fibonacci2ColsShifted<_>>::prove(
            &mut trace,
            &pub_inputs,
            &proof_options,
            StoneProverTranscript::new(&transcript_init_seed),
        )
        .unwrap();
        (proof, pub_inputs, proof_options, transcript_init_seed)
    }

    fn stone_compatibility_case_1_proof() -> StarkProof<Stark252PrimeField, Stark252PrimeField> {
        let (proof, _, _, _) = proof_parts_stone_compatibility_case_1();
        proof
    }

    fn stone_compatibility_case_1_challenges(
    ) -> Challenges<Fibonacci2ColsShifted<Stark252PrimeField>> {
        let (proof, public_inputs, options, seed) = proof_parts_stone_compatibility_case_1();

        let air = Fibonacci2ColsShifted::new(proof.trace_length, &public_inputs, &options);
        let domain = Domain::new(&air);
        Verifier::step_1_replay_rounds_and_recover_challenges(
            &air,
            &proof,
            &domain,
            &mut StoneProverTranscript::new(&seed),
        )
    }

    #[test]
    fn stone_compatibility_case_1_proof_is_valid() {
        let (proof, public_inputs, options, seed) = proof_parts_stone_compatibility_case_1();
        assert!(Verifier::<Fibonacci2ColsShifted<_>>::verify(
            &proof,
            &public_inputs,
            &options,
            StoneProverTranscript::new(&seed)
        ));
    }

    #[test]
    fn stone_compatibility_case_1_trace_commitment() {
        let proof = stone_compatibility_case_1_proof();

        assert_eq!(
            proof.lde_trace_main_merkle_root.to_vec(),
            decode_hex("0eb9dcc0fb1854572a01236753ce05139d392aa3aeafe72abff150fe21175594").unwrap()
        );
    }

    #[test]
    fn stone_compatibility_case_1_composition_poly_challenges() {
        let challenges = stone_compatibility_case_1_challenges();

        assert_eq!(challenges.transition_coeffs[0], FieldElement::one());
        let beta = challenges.transition_coeffs[1];
        assert_eq!(
            beta,
            FieldElement::from_hex_unchecked(
                "86105fff7b04ed4068ecccb8dbf1ed223bd45cd26c3532d6c80a818dbd4fa7"
            ),
        );

        assert_eq!(challenges.boundary_coeffs[0], beta.pow(2u64));
        assert_eq!(challenges.boundary_coeffs[1], beta.pow(3u64));
    }

    #[test]
    fn stone_compatibility_case_1_composition_poly_commitment() {
        let proof = stone_compatibility_case_1_proof();
        // Composition polynomial commitment
        assert_eq!(
            proof.composition_poly_root.to_vec(),
            decode_hex("7cdd8d5fe3bd62254a417e2e260e0fed4fccdb6c9005e828446f645879394f38").unwrap()
        );
    }

    #[test]
    fn stone_compatibility_case_1_out_of_domain_challenge() {
        let challenges = stone_compatibility_case_1_challenges();
        assert_eq!(
            challenges.z,
            FieldElement::from_hex_unchecked(
                "317629e783794b52cd27ac3a5e418c057fec9dd42f2b537cdb3f24c95b3e550"
            )
        );
    }

    #[test]
    fn stone_compatibility_case_1_out_of_domain_trace_evaluation() {
        let proof = stone_compatibility_case_1_proof();

        assert_eq!(
            proof.trace_ood_evaluations.get_row(0)[0],
            FieldElement::from_hex_unchecked(
                "70d8181785336cc7e0a0a1078a79ee6541ca0803ed3ff716de5a13c41684037",
            )
        );
        assert_eq!(
            proof.trace_ood_evaluations.get_row(1)[0],
            FieldElement::from_hex_unchecked(
                "29808fc8b7480a69295e4b61600480ae574ca55f8d118100940501b789c1630",
            )
        );
        assert_eq!(
            proof.trace_ood_evaluations.get_row(0)[1],
            FieldElement::from_hex_unchecked(
                "7d8110f21d1543324cc5e472ab82037eaad785707f8cae3d64c5b9034f0abd2",
            )
        );
        assert_eq!(
            proof.trace_ood_evaluations.get_row(1)[1],
            FieldElement::from_hex_unchecked(
                "1b58470130218c122f71399bf1e04cf75a6e8556c4751629d5ce8c02cc4e62d",
            )
        );
    }

    #[test]
    fn stone_compatibility_case_1_out_of_domain_composition_poly_evaluation() {
        let proof = stone_compatibility_case_1_proof();

        assert_eq!(
            proof.composition_poly_parts_ood_evaluation[0],
            FieldElement::from_hex_unchecked(
                "1c0b7c2275e36d62dfb48c791be122169dcc00c616c63f8efb2c2a504687e85",
            )
        );
    }

    #[test]
    fn stone_compatibility_case_1_deep_composition_poly_challenges() {
        let challenges = stone_compatibility_case_1_challenges();

        // Trace terms coefficients
        assert_eq!(challenges.trace_term_coeffs[0][0], FieldElement::one());
        let gamma = challenges.trace_term_coeffs[0][1];
        assert_eq!(
            &gamma,
            &FieldElement::from_hex_unchecked(
                "a0c79c1c77ded19520873d9c2440451974d23302e451d13e8124cf82fc15dd"
            )
        );
        assert_eq!(&challenges.trace_term_coeffs[1][0], &gamma.pow(2_u64));
        assert_eq!(&challenges.trace_term_coeffs[1][1], &gamma.pow(3_u64));

        // Composition polynomial parts terms coefficient
        assert_eq!(&challenges.gammas[0], &gamma.pow(4_u64));
    }

    #[test]
    fn stone_compatibility_case_1_fri_commit_phase_challenge_0() {
        let challenges = stone_compatibility_case_1_challenges();

        // Challenge to fold FRI polynomial
        assert_eq!(
            challenges.zetas[0],
            FieldElement::from_hex_unchecked(
                "5c6b5a66c9fda19f583f0b10edbaade98d0e458288e62c2fa40e3da2b293cef"
            )
        );
    }

    #[test]
    fn stone_compatibility_case_1_fri_commit_phase_layer_1_commitment() {
        let proof = stone_compatibility_case_1_proof();

        // Commitment of first layer of FRI
        assert_eq!(
            proof.fri_layers_merkle_roots[0].to_vec(),
            decode_hex("327d47da86f5961ee012b2b0e412de16023ffba97c82bfe85102f00daabd49fb").unwrap()
        );
    }

    #[test]
    fn stone_compatibility_case_1_fri_commit_phase_challenge_1() {
        let challenges = stone_compatibility_case_1_challenges();
        assert_eq!(
            challenges.zetas[1],
            FieldElement::from_hex_unchecked(
                "13c337c9dc727bea9eef1f82cab86739f17acdcef562f9e5151708f12891295"
            )
        );
    }

    #[test]
    fn stone_compatibility_case_1_fri_commit_phase_last_value() {
        let proof = stone_compatibility_case_1_proof();

        assert_eq!(
            proof.fri_last_value,
            FieldElement::from_hex_unchecked(
                "43fedf9f9e3d1469309862065c7d7ca0e7e9ce451906e9c01553056f695aec9"
            )
        );
    }

    #[test]
    fn stone_compatibility_case_1_fri_query_iota_challenge() {
        let challenges = stone_compatibility_case_1_challenges();
        assert_eq!(challenges.iotas[0], 1);
    }

    #[test]
    fn stone_compatibility_case_1_fri_query_phase_trace_openings() {
        let proof = stone_compatibility_case_1_proof();

        // Trace Col 0
        assert_eq!(
            proof.deep_poly_openings[0].main_trace_polys.evaluations[0],
            FieldElement::from_hex_unchecked(
                "4de0d56f9cf97dff326c26592fbd4ae9ee756080b12c51cfe4864e9b8734f43"
            )
        );

        // Trace Col 1
        assert_eq!(
            proof.deep_poly_openings[0].main_trace_polys.evaluations[1],
            FieldElement::from_hex_unchecked(
                "1bc1aadf39f2faee64d84cb25f7a95d3dceac1016258a39fc90c9d370e69e8e"
            )
        );

        // Trace Col 0 symmetric
        assert_eq!(
            proof.deep_poly_openings[0].main_trace_polys.evaluations_sym[0],
            FieldElement::from_hex_unchecked(
                "321f2a9063068310cd93d9a6d042b516118a9f7f4ed3ae301b79b16478cb0c6"
            )
        );

        // Trace Col 1 symmetric
        assert_eq!(
            proof.deep_poly_openings[0].main_trace_polys.evaluations_sym[1],
            FieldElement::from_hex_unchecked(
                "643e5520c60d06219b27b34da0856a2c23153efe9da75c6036f362c8f196186"
            )
        );
    }

    #[test]
    fn stone_compatibility_case_1_fri_query_phase_trace_terms_authentication_path() {
        let proof = stone_compatibility_case_1_proof();

        // Trace poly auth path level 1
        assert_eq!(
            proof.deep_poly_openings[0]
                .main_trace_polys
                .proof
                .merkle_path[1]
                .to_vec(),
            decode_hex("91b0c0b24b9d00067b0efab50832b76cf97192091624d42b86740666c5d369e6").unwrap()
        );

        // Trace poly auth path level 2
        assert_eq!(
            proof.deep_poly_openings[0]
                .main_trace_polys
                .proof
                .merkle_path[2]
                .to_vec(),
            decode_hex("993b044db22444c0c0ebf1095b9a51faeb001c9b4dea36abe905f7162620dbbd").unwrap()
        );

        // Trace poly auth path level 3
        assert_eq!(
            proof.deep_poly_openings[0]
                .main_trace_polys
                .proof
                .merkle_path[3]
                .to_vec(),
            decode_hex("5017abeca33fa82576b5c5c2c61792693b48c9d4414a407eef66b6029dae07ea").unwrap()
        );
    }

    #[test]
    fn stone_compatibility_case_1_fri_query_phase_composition_poly_openings() {
        let proof = stone_compatibility_case_1_proof();

        // Composition poly
        assert_eq!(
            proof.deep_poly_openings[0].composition_poly.evaluations[0],
            FieldElement::from_hex_unchecked(
                "2b54852557db698e97253e9d110d60e9bf09f1d358b4c1a96f9f3cf9d2e8755"
            )
        );
        // Composition poly sym
        assert_eq!(
            proof.deep_poly_openings[0].composition_poly.evaluations_sym[0],
            FieldElement::from_hex_unchecked(
                "190f1b0acb7858bd3f5285b68befcf32b436a5f1e3a280e1f42565c1f35c2c3"
            )
        );
    }

    #[test]
    fn stone_compatibility_case_1_fri_query_phase_composition_poly_authentication_path() {
        let proof = stone_compatibility_case_1_proof();

        // Composition poly auth path level 0
        assert_eq!(
            proof.deep_poly_openings[0]
                .composition_poly
                .proof
                .merkle_path[0]
                .to_vec(),
            decode_hex("403b75a122eaf90a298e5d3db2cc7ca096db478078122379a6e3616e72da7546").unwrap()
        );

        // Composition poly auth path level 1
        assert_eq!(
            proof.deep_poly_openings[0]
                .composition_poly
                .proof
                .merkle_path[1]
                .to_vec(),
            decode_hex("07950888c0355c204a1e83ecbee77a0a6a89f93d41cc2be6b39ddd1e727cc965").unwrap()
        );

        // Composition poly auth path level 2
        assert_eq!(
            proof.deep_poly_openings[0]
                .composition_poly
                .proof
                .merkle_path[2]
                .to_vec(),
            decode_hex("58befe2c5de74cc5a002aa82ea219c5b242e761b45fd266eb95521e9f53f44eb").unwrap()
        );
    }

    #[test]
    fn stone_compatibility_case_1_fri_query_phase_query_lengths() {
        let proof = stone_compatibility_case_1_proof();

        assert_eq!(proof.query_list.len(), 1);

        assert_eq!(proof.query_list[0].layers_evaluations_sym.len(), 1);

        assert_eq!(
            proof.query_list[0].layers_auth_paths[0].merkle_path.len(),
            2
        );
    }

    #[test]
    fn stone_compatibility_case_1_fri_query_phase_layer_1_evaluation_symmetric() {
        let proof = stone_compatibility_case_1_proof();

        assert_eq!(
            proof.query_list[0].layers_evaluations_sym[0],
            FieldElement::from_hex_unchecked(
                "0684991e76e5c08db17f33ea7840596be876d92c143f863e77cad10548289fd0"
            )
        );
    }

    #[test]
    fn stone_compatibility_case_1_fri_query_phase_layer_1_authentication_path() {
        let proof = stone_compatibility_case_1_proof();

        // FRI layer 1 auth path level 0
        assert_eq!(
            proof.query_list[0].layers_auth_paths[0].merkle_path[0].to_vec(),
            decode_hex("0683622478e9e93cc2d18754872f043619f030b494d7ec8e003b1cbafe83b67b").unwrap()
        );

        // FRI layer 1 auth path level 1
        assert_eq!(
            proof.query_list[0].layers_auth_paths[0].merkle_path[1].to_vec(),
            decode_hex("7985d945abe659a7502698051ec739508ed6bab594984c7f25e095a0a57a2e55").unwrap()
        );
    }

    fn proof_parts_stone_compatibility_case_2() -> (
        StarkProof<Stark252PrimeField, Stark252PrimeField>,
        fibonacci_2_cols_shifted::PublicInputs<Stark252PrimeField>,
        ProofOptions,
        [u8; 4],
    ) {
        let mut trace = fibonacci_2_cols_shifted::compute_trace(FieldElement::from(12345), 512);

        let claimed_index = 420;
        let col = 0;
        let claimed_value = *trace.get_main(claimed_index, col);
        let mut proof_options = ProofOptions::default_test_options();
        proof_options.blowup_factor = 1 << 6;
        proof_options.coset_offset = 3;
        proof_options.grinding_factor = 0;
        proof_options.fri_number_of_queries = 1;

        let pub_inputs = fibonacci_2_cols_shifted::PublicInputs {
            claimed_value,
            claimed_index,
        };

        let transcript_init_seed = [0xfa, 0xfa, 0xfa, 0xee];

        let proof = Prover::<Fibonacci2ColsShifted<_>>::prove(
            &mut trace,
            &pub_inputs,
            &proof_options,
            StoneProverTranscript::new(&transcript_init_seed),
        )
        .unwrap();
        (proof, pub_inputs, proof_options, transcript_init_seed)
    }

    fn stone_compatibility_case_2_proof() -> StarkProof<Stark252PrimeField, Stark252PrimeField> {
        let (proof, _, _, _) = proof_parts_stone_compatibility_case_2();
        proof
    }

    fn stone_compatibility_case_2_challenges(
    ) -> Challenges<Fibonacci2ColsShifted<Stark252PrimeField>> {
        let (proof, public_inputs, options, seed) = proof_parts_stone_compatibility_case_2();

        let air = Fibonacci2ColsShifted::new(proof.trace_length, &public_inputs, &options);
        let domain = Domain::new(&air);
        Verifier::step_1_replay_rounds_and_recover_challenges(
            &air,
            &proof,
            &domain,
            &mut StoneProverTranscript::new(&seed),
        )
    }

    #[test]
    fn stone_compatibility_case_2_trace_commitment() {
        let proof = stone_compatibility_case_2_proof();

        assert_eq!(
            proof.lde_trace_main_merkle_root.to_vec(),
            decode_hex("6d31dd00038974bde5fe0c5e3a765f8ddc822a5df3254fca85a1950ae0208cbe").unwrap()
        );
    }

    #[test]
    fn stone_compatibility_case_2_fri_query_iota_challenge() {
        let challenges = stone_compatibility_case_2_challenges();
        assert_eq!(challenges.iotas[0], 4239);
    }

    #[test]
    fn stone_compatibility_case_2_fri_query_phase_layer_7_evaluation_symmetric() {
        let proof = stone_compatibility_case_2_proof();

        assert_eq!(
            proof.query_list[0].layers_evaluations_sym[7],
            FieldElement::from_hex_unchecked(
                "7aa40c5a4e30b44fee5bcc47c54072a435aa35c1a31b805cad8126118cc6860"
            )
        );
    }

    #[test]
    fn stone_compatibility_case_2_fri_query_phase_layer_8_authentication_path() {
        let proof = stone_compatibility_case_2_proof();

        // FRI layer 7 auth path level 5
        assert_eq!(
            proof.query_list[0].layers_auth_paths[7].merkle_path[5].to_vec(),
            decode_hex("f12f159b548ca2c571a270870d43e7ec2ead78b3e93b635738c31eb9bcda3dda").unwrap()
        );
    }
}
