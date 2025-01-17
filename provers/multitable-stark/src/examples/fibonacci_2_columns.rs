use std::marker::PhantomData;

use super::simple_fibonacci::FibonacciPublicInputs;
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
struct FibTransition1<F: IsFFTField> {
    phantom: PhantomData<F>,
}

impl<F: IsFFTField> FibTransition1<F> {
    pub fn new() -> Self {
        Self {
            phantom: PhantomData,
        }
    }
}

impl<F> TransitionConstraint<F, F> for FibTransition1<F>
where
    F: IsFFTField + Send + Sync,
{
    fn degree(&self) -> usize {
        1
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

        // s_{0, i+1} = s_{0, i} + s_{1, i}
        let s0_0 = first_step.get_main_evaluation_element(0, 0);
        let s0_1 = first_step.get_main_evaluation_element(0, 1);
        let s1_0 = second_step.get_main_evaluation_element(0, 0);

        let res = s1_0 - s0_0 - s0_1;

        transition_evaluations[self.constraint_idx()] = res;
    }
}

#[derive(Clone)]
struct FibTransition2<F: IsFFTField> {
    phantom: PhantomData<F>,
}

impl<F: IsFFTField> FibTransition2<F> {
    pub fn new() -> Self {
        Self {
            phantom: PhantomData,
        }
    }
}

impl<F> TransitionConstraint<F, F> for FibTransition2<F>
where
    F: IsFFTField + Send + Sync,
{
    fn degree(&self) -> usize {
        1
    }

    fn constraint_idx(&self) -> usize {
        1
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

        // s_{1, i+1} = s_{1, i} + s_{0, i+1}
        let s0_1 = first_step.get_main_evaluation_element(0, 1);
        let s1_0 = second_step.get_main_evaluation_element(0, 0);
        let s1_1 = second_step.get_main_evaluation_element(0, 1);

        let res = s1_1 - s0_1 - s1_0;

        transition_evaluations[self.constraint_idx()] = res;
    }
}

pub struct Fibonacci2ColsAIR<F>
where
    F: IsFFTField,
{
    context: AirContext,
    trace_length: usize,
    pub_inputs: FibonacciPublicInputs<F>,
    constraints: Vec<Box<dyn TransitionConstraint<F, F>>>,
}

/// The AIR for to a 2 column trace, where the columns form a Fibonacci sequence when
/// stacked in row-major order.
impl<F> AIR for Fibonacci2ColsAIR<F>
where
    F: IsFFTField + Send + Sync + 'static,
{
    type Field = F;
    type FieldExtension = F;
    type PublicInputs = FibonacciPublicInputs<Self::Field>;

    const STEP_SIZE: usize = 1;

    fn new(
        trace_length: usize,
        pub_inputs: &Self::PublicInputs,
        proof_options: &ProofOptions,
    ) -> Self {
        let constraints: Vec<Box<dyn TransitionConstraint<Self::Field, Self::FieldExtension>>> = vec![
            Box::new(FibTransition1::new()),
            Box::new(FibTransition2::new()),
        ];

        let context = AirContext {
            proof_options: proof_options.clone(),
            transition_offsets: vec![0, 1],
            num_transition_constraints: constraints.len(),
            trace_columns: 2,
        };

        Self {
            trace_length,
            context,
            constraints,
            pub_inputs: pub_inputs.clone(),
        }
    }

    fn boundary_constraints(
        &self,
        _rap_challenges: &[FieldElement<Self::Field>],
    ) -> BoundaryConstraints<Self::Field> {
        let a0 = BoundaryConstraint::new_main(0, 0, self.pub_inputs.a0.clone());
        let a1 = BoundaryConstraint::new_main(1, 0, self.pub_inputs.a1.clone());

        BoundaryConstraints::from_constraints(vec![a0, a1])
    }

    fn transition_constraints(&self) -> &Vec<Box<dyn TransitionConstraint<F, F>>> {
        &self.constraints
    }

    fn context(&self) -> &AirContext {
        &self.context
    }

    fn composition_poly_degree_bound(&self) -> usize {
        self.trace_length()
    }

    fn trace_length(&self) -> usize {
        self.trace_length
    }

    fn trace_layout(&self) -> (usize, usize) {
        (2, 0)
    }

    fn pub_inputs(&self) -> &Self::PublicInputs {
        &self.pub_inputs
    }
}

pub fn compute_trace<F: IsFFTField>(
    initial_values: [FieldElement<F>; 2],
    trace_length: usize,
) -> TraceTable<F, F> {
    let mut ret1: Vec<FieldElement<F>> = vec![];
    let mut ret2: Vec<FieldElement<F>> = vec![];

    ret1.push(initial_values[0].clone());
    ret2.push(initial_values[1].clone());

    for i in 1..(trace_length) {
        let new_val = ret1[i - 1].clone() + ret2[i - 1].clone();
        ret1.push(new_val.clone());
        ret2.push(new_val + ret2[i - 1].clone());
    }

    TraceTable::from_columns_main(vec![ret1, ret2], 1)
}
