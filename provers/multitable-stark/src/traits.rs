use std::collections::HashMap;

use lambdaworks_crypto::fiat_shamir::is_transcript::IsTranscript;
use lambdaworks_math::{
    field::{
        element::FieldElement,
        traits::{IsFFTField, IsField, IsSubFieldOf},
    },
    polynomial::Polynomial,
};

use crate::{constraints::transition::TransitionConstraint, domain::Domain};

use super::{
    constraints::boundary::BoundaryConstraints, context::AirContext, frame::Frame,
    proof::options::ProofOptions, trace::TraceTable,
};

type ZerofierGroupKey = (usize, usize, Option<usize>, Option<usize>, usize);

/// This enum is necessary because, while both the prover and verifier perform the same operations
///  to compute transition constraints, their frames differ.
///  The prover uses a frame containing elements from both the base field and its extension
/// (common when working with small fields and challengers in the extension).
/// In contrast, the verifier, lacking access to the trace and relying solely on evaluations at the challengers,
/// works with a frame that contains only elements from the extension.
pub enum TransitionEvaluationContext<'a, F, E>
where
    F: IsSubFieldOf<E>,
    E: IsField,
{
    Prover {
        frame: &'a Frame<'a, F, E>,
        periodic_values: &'a [FieldElement<F>],
        rap_challenges: &'a [FieldElement<E>],
    },
    Verifier {
        frame: &'a Frame<'a, E, E>,
        periodic_values: &'a [FieldElement<E>],
        rap_challenges: &'a [FieldElement<E>],
    },
}

impl<'a, F, E> TransitionEvaluationContext<'a, F, E>
where
    F: IsSubFieldOf<E>,
    E: IsField,
{
    pub fn new_prover(
        frame: &'a Frame<'a, F, E>,
        periodic_values: &'a [FieldElement<F>],
        rap_challenges: &'a [FieldElement<E>],
    ) -> Self {
        Self::Prover {
            frame,
            periodic_values,
            rap_challenges,
        }
    }

    pub fn new_verifier(
        frame: &'a Frame<'a, E, E>,
        periodic_values: &'a [FieldElement<E>],
        rap_challenges: &'a [FieldElement<E>],
    ) -> Self {
        Self::Verifier {
            frame,
            periodic_values,
            rap_challenges,
        }
    }
}

/// AIR is a representation of the Constraints
pub trait AIR {
    type Field: IsFFTField + IsSubFieldOf<Self::FieldExtension> + Send + Sync;
    type FieldExtension: IsField + Send + Sync;
    type PublicInputs;

    // How many rows are one row in Cairo.
    // TODO: Remove this.
    const STEP_SIZE: usize;

    fn new(
        trace_length: Vec<usize>,
        pub_inputs: Vec<&Self::PublicInputs>,
        // TODO: Same proof options for all the Tables?
        proof_options: &ProofOptions,
    ) -> Self;

    fn build_auxiliary_trace(
        &self,
        _main_trace: Vec<&mut TraceTable<Self::Field, Self::FieldExtension>>,
        _rap_challenges: Vec<&[FieldElement<Self::FieldExtension>]>,
    ) where
        Self::FieldExtension: IsFFTField,
    {
    }

    fn build_rap_challenges(
        &self,
        _transcript: &mut impl IsTranscript<Self::FieldExtension>,
    ) -> Vec<FieldElement<Self::FieldExtension>> {
        Vec::new()
    }

    /// Returns the amount main trace columns and auxiliary trace columns
    fn trace_layout(&self) -> Vec<(usize, usize)>;
    // check if the option (Vec<usize>,Vec<usize>) is better

    fn has_trace_interaction(&self) -> Vec<bool> {
        self.trace_layout()
            .into_iter()
            .map(|(_, aux_trace_columns)| aux_trace_columns != 0)
            .collect()
    }

    // // One table version of has_trace_interaction. You can find it in stark > traits.rs:
    // fn has_trace_interaction(&self) -> bool {
    //     let (_main_trace_columns, aux_trace_columns) = self.trace_layout();
    //     aux_trace_columns != 0
    // }

    fn num_auxiliary_rap_columns(&self) -> Vec<usize> {
        self.trace_layout()
            .into_iter()
            .map(|(_, aux_trace_columns)| aux_trace_columns)
            .collect()
    }

    fn composition_poly_degree_bound(&self) -> Vec<usize>;

    /// The method called by the prover to evaluate the transitions corresponding to an evaluation frame.
    /// In the case of the prover, the main evaluation table of the frame takes values in
    /// `Self::Field`, since they are the evaluations of the main trace at the LDE domain.
    /// In the case of the verifier, the frame take elements of Self::FieldExtension.
    fn compute_transition(
        &self,
        evaluation_context: Vec<&TransitionEvaluationContext<Self::Field, Self::FieldExtension>>,
    ) -> Vec<FieldElement<Self::FieldExtension>> {
        let mut evaluations =
            vec![FieldElement::<Self::FieldExtension>::zero(); self.num_transition_constraints()];
        self.transition_constraints()
            .iter()
            .zip(evaluation_context)
            .for_each(|(c, ctx)| c.evaluate(ctx, &mut evaluations));

        evaluations
    }

    fn boundary_constraints(
        &self,
        rap_challenges: &[FieldElement<Self::FieldExtension>],
    ) -> BoundaryConstraints<Self::FieldExtension>;

    fn context(&self) -> &AirContext;

    fn trace_length(&self) -> usize;

    fn options(&self) -> &ProofOptions {
        &self.context().proof_options
    }

    fn blowup_factor(&self) -> u8 {
        self.options().blowup_factor
    }

    fn coset_offset(&self) -> FieldElement<Self::Field> {
        FieldElement::from(self.options().coset_offset)
    }

    fn trace_primitive_root(&self) -> FieldElement<Self::Field> {
        let trace_length = self.trace_length();
        let root_of_unity_order = u64::from(trace_length.trailing_zeros());

        Self::Field::get_primitive_root_of_unity(root_of_unity_order).unwrap()
    }

    fn num_transition_constraints(&self) -> usize {
        self.context().num_transition_constraints
    }

    fn pub_inputs(&self) -> &Self::PublicInputs;

    fn get_periodic_column_values(&self) -> Vec<Vec<FieldElement<Self::Field>>> {
        vec![]
    }

    fn get_periodic_column_polynomials(&self) -> Vec<Polynomial<FieldElement<Self::Field>>> {
        let mut result = Vec::new();
        for periodic_column in self.get_periodic_column_values() {
            let values: Vec<_> = periodic_column
                .iter()
                .cycle()
                .take(self.trace_length())
                .cloned()
                .collect();
            let poly =
                Polynomial::<FieldElement<Self::Field>>::interpolate_fft::<Self::Field>(&values)
                    .unwrap();
            result.push(poly);
        }
        result
    }

    fn transition_constraints(
        &self,
    ) -> &Vec<Vec<Box<dyn TransitionConstraint<Self::Field, Self::FieldExtension>>>>;

    fn transition_zerofier_evaluations(
        &self,
        domain: &Domain<Self::Field>,
    ) -> Vec<Vec<FieldElement<Self::Field>>> {
        let mut evals = vec![Vec::new(); self.num_transition_constraints()];

        let mut zerofier_groups: HashMap<ZerofierGroupKey, Vec<FieldElement<Self::Field>>> =
            HashMap::new();

        self.transition_constraints().iter().for_each(|c| {
            let period = c.period();
            let offset = c.offset();
            let exemptions_period = c.exemptions_period();
            let periodic_exemptions_offset = c.periodic_exemptions_offset();
            let end_exemptions = c.end_exemptions();

            // This hashmap is used to avoid recomputing with an fft the same zerofier evaluation
            // If there are multiple domain and subdomains it can be further optimized
            // as to share computation between them

            let zerofier_group_key = (
                period,
                offset,
                exemptions_period,
                periodic_exemptions_offset,
                end_exemptions,
            );
            zerofier_groups
                .entry(zerofier_group_key)
                .or_insert_with(|| c.zerofier_evaluations_on_extended_domain(domain));

            let zerofier_evaluations = zerofier_groups.get(&zerofier_group_key).unwrap();
            evals[c.constraint_idx()] = zerofier_evaluations.clone();
        });

        evals
    }
}
