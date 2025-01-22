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

pub trait Chip1AIR: AIR {
    /// Chip1-specific transition constraints
    fn chip1_transition_constraints(
        &self,
    ) -> Vec<Box<dyn TransitionConstraint<Self::Field, Self::FieldExtension>>>;
}

pub trait Chip2AIR: AIR {
    /// Chip2-specific transition constraints
    fn chip2_transition_constraints(
        &self,
    ) -> Vec<Box<dyn TransitionConstraint<Self::Field, Self::FieldExtension>>>;
}

pub struct CombinedAIR<C1, C2, F, E>
where
    C1: Chip1AIR<Field = F, FieldExtension = E>,
    C2: Chip2AIR<Field = F, FieldExtension = E>,
    F: IsFFTField + IsSubFieldOf<E> + Send + Sync,
    E: IsField + Send + Sync,
{
    chip1: C1,
    chip2: C2,
    interaction_constraints: Vec<Box<dyn TransitionConstraint<F, E>>>,
    //context: AirContext,
    pub_inputs: CombinedPublicInputs<C1::PublicInputs, C2::PublicInputs>,
    trace_length: usize,
}

impl<C1, C2, F, E> CombinedAIR<C1, C2, F, E>
where
    C1: Chip1AIR<Field = F, FieldExtension = E>,
    C2: Chip2AIR<Field = F, FieldExtension = E>,
    F: IsFFTField + IsSubFieldOf<E> + Send + Sync,
    E: IsField + Send + Sync,
{
    pub fn new(
        chip1: C1,
        chip2: C2,
        interaction_constraints: Vec<Box<dyn TransitionConstraint<F, E>>>,
        trace_length: usize,
        pub_inputs: &CombinedPublicInputs<C1::PublicInputs, C2::PublicInputs>,
        proof_options: &ProofOptions,
    ) -> Self {
        let num_constraints = chip1.num_transition_constraints()
            + chip2.num_transition_constraints()
            + interaction_constraints.len();

        let context = AirContext {
            proof_options: proof_options.clone(),
            num_transition_constraints: num_constraints,
        };

        Self {
            chip1,
            chip2,
            interaction_constraints,
            context,
            pub_inputs: pub_inputs.clone(),
            trace_length,
        }
    }
}
impl<C1, C2, F, E> AIR for CombinedAIR<C1, C2, F, E>
where
    C1: Chip1AIR<Field = F, FieldExtension = E>,
    C2: Chip2AIR<Field = F, FieldExtension = E>,
    F: IsFFTField + IsSubFieldOf<E> + Send + Sync,
    E: IsField + Send + Sync,
{
    type Field = F;
    type FieldExtension = E;
    type PublicInputs = CombinedPublicInputs<C1::PublicInputs, C2::PublicInputs>;

    const STEP_SIZE: usize = 1;

    fn new(
        trace_length: usize,
        pub_inputs: &Self::PublicInputs,
        proof_options: &ProofOptions,
    ) -> Self {
        unimplemented!("Use CombinedAIR::new instead")
    }

    fn trace_layout(&self) -> (usize, usize) {
        let (main1, aux1) = self.chip1.trace_layout();
        let (main2, aux2) = self.chip2.trace_layout();
        (main1 + main2, aux1 + aux2)
    }

    fn transition_constraints(
        &self,
    ) -> &Vec<Box<dyn TransitionConstraint<Self::Field, Self::FieldExtension>>> {
        &self.interaction_constraints
    }

    fn boundary_constraints(
        &self,
        rap_challenges: &[FieldElement<Self::FieldExtension>],
    ) -> BoundaryConstraints<Self::FieldExtension> {
        let mut boundary_constraints = self.chip1.boundary_constraints(rap_challenges);
        boundary_constraints.extend(self.chip2.boundary_constraints(rap_challenges));
        boundary_constraints
    }

    fn context(&self) -> &AirContext {
        &self.context
    }

    fn trace_length(&self) -> usize {
        self.trace_length
    }

    fn pub_inputs(&self) -> &Self::PublicInputs {
        &self.pub_inputs
    }
}
pub struct InteractionConstraint<F, E> {
    chip1_column: usize,
    chip2_column: usize,
    _phantom: std::marker::PhantomData<(F, E)>,
}

impl<F, E> TransitionConstraint<F, E> for InteractionConstraint<F, E>
where
    F: IsFFTField + IsSubFieldOf<E> + Send + Sync,
    E: IsField + Send + Sync,
{
    fn evaluate(
        &self,
        evaluation_context: &TransitionEvaluationContext<F, E>,
        evaluations: &mut Vec<FieldElement<E>>,
    ) {
        match evaluation_context {
            TransitionEvaluationContext::Prover {
                frame,
                periodic_values,
                rap_challenges,
            } => {
                let chip1_value = frame.get_main_trace_column(self.chip1_column);
                let chip2_value = frame.get_main_trace_column(self.chip2_column);
                evaluations[0] = chip1_value[0] - chip2_value[0];
            }
            TransitionEvaluationContext::Verifier {
                frame,
                periodic_values,
                rap_challenges,
            } => {
                // Similar logic for the verifier
            }
        }
    }

    // Implement other required methods (period, offset, etc.)
}
pub trait InteractionAIR {
    type Field: IsFFTField + IsSubFieldOf<Self::FieldExtension> + Send + Sync;
    type FieldExtension: IsField + Send + Sync;

    /// Returns interaction-specific transition constraints
    fn interaction_transition_constraints(
        &self,
    ) -> Vec<Box<dyn TransitionConstraint<Self::Field, Self::FieldExtension>>>;

    /// Combines constraints from Chip1, Chip2, and interaction
    fn combined_transition_constraints(
        &self,
    ) -> Vec<Box<dyn TransitionConstraint<Self::Field, Self::FieldExtension>>> {
        let mut constraints = self.chip1_transition_constraints();
        constraints.extend(self.chip2_transition_constraints());
        constraints.extend(self.interaction_transition_constraints());
        constraints
    }

    /// Returns Chip1's transition constraints
    fn chip1_transition_constraints(
        &self,
    ) -> Vec<Box<dyn TransitionConstraint<Self::Field, Self::FieldExtension>>>;

    /// Returns Chip2's transition constraints
    fn chip2_transition_constraints(
        &self,
    ) -> Vec<Box<dyn TransitionConstraint<Self::Field, Self::FieldExtension>>>;
}

/// AIR is a representation of the Constraints
pub trait AIR {
    type Field: IsFFTField + IsSubFieldOf<Self::FieldExtension> + Send + Sync;
    type FieldExtension: IsField + Send + Sync;
    type PublicInputs;

    const STEP_SIZE: usize;

    fn new(
        trace_length: usize,
        pub_inputs: &Self::PublicInputs,
        proof_options: &ProofOptions,
    ) -> Self;

    fn build_auxiliary_trace(
        &self,
        _main_trace: &mut TraceTable<Self::Field, Self::FieldExtension>,
        _rap_challenges: &[FieldElement<Self::FieldExtension>],
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
    fn trace_layout(&self) -> (usize, usize);

    fn has_trace_interaction(&self) -> bool {
        let (_main_trace_columns, aux_trace_columns) = self.trace_layout();
        aux_trace_columns != 0
    }

    fn num_auxiliary_rap_columns(&self) -> usize {
        self.trace_layout().1
    }

    fn composition_poly_degree_bound(&self) -> usize;

    /// The method called by the prover to evaluate the transitions corresponding to an evaluation frame.
    /// In the case of the prover, the main evaluation table of the frame takes values in
    /// `Self::Field`, since they are the evaluations of the main trace at the LDE domain.
    /// In the case of the verifier, the frame take elements of Self::FieldExtension.
    fn compute_transition(
        &self,
        evaluation_context: &TransitionEvaluationContext<Self::Field, Self::FieldExtension>,
    ) -> Vec<FieldElement<Self::FieldExtension>> {
        let mut evaluations =
            vec![FieldElement::<Self::FieldExtension>::zero(); self.num_transition_constraints()];
        self.transition_constraints()
            .iter()
            .for_each(|c| c.evaluate(evaluation_context, &mut evaluations));

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
    ) -> &Vec<Box<dyn TransitionConstraint<Self::Field, Self::FieldExtension>>>;

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
