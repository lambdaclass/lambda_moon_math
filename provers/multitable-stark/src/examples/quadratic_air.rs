use std::marker::PhantomData;

use crate::{
    constraints::{
        boundary::{BoundaryConstraint, BoundaryConstraints},
        transition::TransitionConstraint,
    },
    context::AirContext,
    proof::options::ProofOptions,
    trace::TraceTable,
    traits::{TransitionEvaluationContext, AIR},
};
use lambdaworks_math::field::{element::FieldElement, traits::IsFFTField};

#[derive(Clone)]
struct QuadraticConstraint<F: IsFFTField> {
    phantom: PhantomData<F>,
}

impl<F: IsFFTField> QuadraticConstraint<F> {
    pub fn new() -> Self {
        Self {
            phantom: PhantomData,
        }
    }
}

impl<F> TransitionConstraint<F, F> for QuadraticConstraint<F>
where
    F: IsFFTField + Send + Sync,
{
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        0
    }

    fn end_exemptions(&self) -> usize {
        1
    }

    fn evaluate(
        &self,
        evaluation_context: &TransitionEvaluationContext<F, F>,
        transition_evaluations: &mut [FieldElement<F>],
    ) {
        let (frame, _periodic_values, _rap_challenges) = match evaluation_context {
            TransitionEvaluationContext::Prover {
                frame,
                periodic_values,
                rap_challenges,
            }
            | TransitionEvaluationContext::Verifier {
                frame,
                periodic_values,
                rap_challenges,
            } => (frame, periodic_values, rap_challenges),
        };

        let first_step = frame.get_evaluation_step(0);
        let second_step = frame.get_evaluation_step(1);

        let x = first_step.get_main_evaluation_element(0, 0);
        let x_squared = second_step.get_main_evaluation_element(0, 0);

        let res = x_squared - x * x;

        transition_evaluations[self.constraint_idx()] = res;
    }
}

pub struct QuadraticAIR<F>
where
    F: IsFFTField,
{
    context: AirContext,
    trace_length: usize,
    pub_inputs: QuadraticPublicInputs<F>,
    constraints: Vec<Box<dyn TransitionConstraint<F, F>>>,
}

#[derive(Clone, Debug)]
pub struct QuadraticPublicInputs<F>
where
    F: IsFFTField,
{
    pub a0: FieldElement<F>,
}

impl<F> AIR for QuadraticAIR<F>
where
    F: IsFFTField + Send + Sync + 'static,
{
    type Field = F;
    type FieldExtension = F;
    type PublicInputs = QuadraticPublicInputs<Self::Field>;

    const STEP_SIZE: usize = 1;

    fn new(
        trace_length: usize,
        pub_inputs: &Self::PublicInputs,
        proof_options: &ProofOptions,
    ) -> Self {
        let constraints: Vec<Box<dyn TransitionConstraint<Self::Field, Self::FieldExtension>>> =
            vec![Box::new(QuadraticConstraint::new())];

        let context = AirContext {
            proof_options: proof_options.clone(),
            trace_columns: 1,
            transition_offsets: vec![0, 1],
            num_transition_constraints: constraints.len(),
        };

        Self {
            trace_length,
            context,
            pub_inputs: pub_inputs.clone(),
            constraints,
        }
    }

    fn boundary_constraints(
        &self,
        _rap_challenges: &[FieldElement<Self::Field>],
    ) -> BoundaryConstraints<Self::Field> {
        let a0 = BoundaryConstraint::new_simple_main(0, self.pub_inputs.a0.clone());

        BoundaryConstraints::from_constraints(vec![a0])
    }

    fn transition_constraints(
        &self,
    ) -> &Vec<Box<dyn TransitionConstraint<Self::Field, Self::FieldExtension>>> {
        &self.constraints
    }

    fn context(&self) -> &AirContext {
        &self.context
    }

    fn composition_poly_degree_bound(&self) -> usize {
        2 * self.trace_length()
    }

    fn trace_layout(&self) -> (usize, usize) {
        (1, 0)
    }

    fn trace_length(&self) -> usize {
        self.trace_length
    }

    fn pub_inputs(&self) -> &Self::PublicInputs {
        &self.pub_inputs
    }
}

pub fn quadratic_trace<F: IsFFTField>(
    initial_value: FieldElement<F>,
    trace_length: usize,
) -> TraceTable<F, F> {
    let mut ret: Vec<FieldElement<F>> = vec![];

    ret.push(initial_value);

    for i in 1..(trace_length) {
        ret.push(ret[i - 1].clone() * ret[i - 1].clone());
    }

    TraceTable::from_columns_main(vec![ret], 1)
}
