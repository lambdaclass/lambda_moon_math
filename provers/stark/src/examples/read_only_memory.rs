use std::{marker::PhantomData, ops::Div};

use crate::{
    constraints::{
        boundary::{BoundaryConstraint, BoundaryConstraints},
        transition::TransitionConstraint,
    },
    context::AirContext,
    frame::Frame,
    proof::options::ProofOptions,
    trace::TraceTable,
    traits::AIR,
};
use lambdaworks_crypto::fiat_shamir::is_transcript::IsTranscript;
use lambdaworks_math::{
    field::{element::FieldElement, traits::IsFFTField},
    helpers::resize_to_next_power_of_two,
    traits::ByteConversion,
};

#[derive(Clone)]
struct FibConstraint<F: IsFFTField> {
    phantom: PhantomData<F>,
}

impl<F: IsFFTField> FibConstraint<F> {
    pub fn new() -> Self {
        Self {
            phantom: PhantomData,
        }
    }
}

impl<F> TransitionConstraint<F, F> for FibConstraint<F>
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
        // NOTE: This is hard-coded for the example of steps = 16 in the integration tests.
        // If that number changes in the test, this should be changed too or the test will fail.
        3 + 32 - 16 - 1
    }

    fn evaluate(
        &self,
        frame: &Frame<F, F>,
        transition_evaluations: &mut [FieldElement<F>],
        _periodic_values: &[FieldElement<F>],
        _rap_challenges: &[FieldElement<F>],
    ) {
        let first_step = frame.get_evaluation_step(0);
        let second_step = frame.get_evaluation_step(1);
        let third_step = frame.get_evaluation_step(2);

        let a0 = first_step.get_main_evaluation_element(0, 0);
        let a1 = second_step.get_main_evaluation_element(0, 0);
        let a2 = third_step.get_main_evaluation_element(0, 0);

        let res = a2 - a1 - a0;

        transition_evaluations[self.constraint_idx()] = res;
    }
}

#[derive(Clone)]
struct PermutationConstraint<F: IsFFTField> {
    phantom: PhantomData<F>,
}

impl<F: IsFFTField> PermutationConstraint<F> {
    pub fn new() -> Self {
        Self {
            phantom: PhantomData,
        }
    }
}

impl<F> TransitionConstraint<F, F> for PermutationConstraint<F>
where
    F: IsFFTField + Send + Sync,
{
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        1
    }

    fn end_exemptions(&self) -> usize {
        1
    }

    fn evaluate(
        &self,
        frame: &Frame<F, F>,
        transition_evaluations: &mut [FieldElement<F>],
        _periodic_values: &[FieldElement<F>],
        rap_challenges: &[FieldElement<F>],
    ) {
        let first_step = frame.get_evaluation_step(0);
        let second_step = frame.get_evaluation_step(1);

        // Auxiliary constraints
        let z_i = first_step.get_aux_evaluation_element(0, 0);
        let z_i_plus_one = second_step.get_aux_evaluation_element(0, 0);
        let gamma = &rap_challenges[0];

        let a_i = first_step.get_main_evaluation_element(0, 0);
        let b_i = first_step.get_main_evaluation_element(0, 1);

        let res = z_i_plus_one * (b_i + gamma) - z_i * (a_i + gamma);

        transition_evaluations[self.constraint_idx()] = res;
    }
}

pub struct FibonacciRAP<F>
where
    F: IsFFTField,
{
    context: AirContext,
    trace_length: usize,
    pub_inputs: FibonacciRAPPublicInputs<F>,
    transition_constraints: Vec<Box<dyn TransitionConstraint<F, F>>>,
}

#[derive(Clone, Debug)]
pub struct FibonacciRAPPublicInputs<F>
where
    F: IsFFTField,
{
    pub steps: usize,
    pub a0: FieldElement<F>,
    pub a1: FieldElement<F>,
}

impl<F> AIR for FibonacciRAP<F>
where
    F: IsFFTField + Send + Sync + 'static,
    FieldElement<F>: ByteConversion,
{
    type Field = F;
    type FieldExtension = F;
    type PublicInputs = FibonacciRAPPublicInputs<Self::Field>;

    const STEP_SIZE: usize = 1;

    fn new(
        trace_length: usize,
        pub_inputs: &Self::PublicInputs,
        proof_options: &ProofOptions,
    ) -> Self {
        let transition_constraints: Vec<
            Box<dyn TransitionConstraint<Self::Field, Self::FieldExtension>>,
        > = vec![
            Box::new(FibConstraint::new()),
            Box::new(PermutationConstraint::new()),
        ];

        let exemptions = 3 + trace_length - pub_inputs.steps - 1;

        let context = AirContext {
            proof_options: proof_options.clone(),
            trace_columns: 3,
            transition_offsets: vec![0, 1, 2],
            transition_exemptions: vec![exemptions, 1],
            num_transition_constraints: transition_constraints.len(),
        };

        Self {
            context,
            trace_length,
            pub_inputs: pub_inputs.clone(),
            transition_constraints,
        }
    }

    fn build_auxiliary_trace(
        &self,
        main_trace: &TraceTable<Self::Field>,
        challenges: &[FieldElement<F>],
    ) -> TraceTable<Self::Field> {
        let main_segment_cols = main_trace.columns();
        let not_perm = &main_segment_cols[0];
        let perm = &main_segment_cols[1];
        let gamma = &challenges[0];

        let trace_len = main_trace.n_rows();

        let mut aux_col = Vec::new();
        for i in 0..trace_len {
            if i == 0 {
                aux_col.push(FieldElement::<Self::Field>::one());
            } else {
                let z_i = &aux_col[i - 1];
                let n_p_term = not_perm[i - 1].clone() + gamma;
                let p_term = &perm[i - 1] + gamma;

                aux_col.push(z_i * n_p_term.div(p_term));
            }
        }
        TraceTable::from_columns(vec![aux_col], 0, 1)
    }

    fn build_rap_challenges(
        &self,
        transcript: &mut impl IsTranscript<Self::Field>,
    ) -> Vec<FieldElement<Self::FieldExtension>> {
        vec![transcript.sample_field_element()]
    }

    fn trace_layout(&self) -> (usize, usize) {
        (2, 1)
    }

    fn boundary_constraints(
        &self,
        _rap_challenges: &[FieldElement<Self::FieldExtension>],
    ) -> BoundaryConstraints<Self::FieldExtension> {
        // Main boundary constraints
        let a0 =
            BoundaryConstraint::new_simple_main(0, FieldElement::<Self::FieldExtension>::one());
        let a1 =
            BoundaryConstraint::new_simple_main(1, FieldElement::<Self::FieldExtension>::one());

        // Auxiliary boundary constraints
        let a0_aux = BoundaryConstraint::new_aux(0, 0, FieldElement::<Self::FieldExtension>::one());

        BoundaryConstraints::from_constraints(vec![a0, a1, a0_aux])
        // BoundaryConstraints::from_constraints(vec![a0, a1])
    }

    fn transition_constraints(
        &self,
    ) -> &Vec<Box<dyn TransitionConstraint<Self::Field, Self::FieldExtension>>> {
        &self.transition_constraints
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

    fn pub_inputs(&self) -> &Self::PublicInputs {
        &self.pub_inputs
    }

    fn compute_transition_verifier(
        &self,
        frame: &Frame<Self::FieldExtension, Self::FieldExtension>,
        periodic_values: &[FieldElement<Self::FieldExtension>],
        rap_challenges: &[FieldElement<Self::FieldExtension>],
    ) -> Vec<FieldElement<Self::Field>> {
        self.compute_transition_prover(frame, periodic_values, rap_challenges)
    }
}

pub fn fibonacci_rap_trace<F: IsFFTField>(
    initial_values: [FieldElement<F>; 2],
    trace_length: usize,
) -> TraceTable<F> {
    let mut fib_seq: Vec<FieldElement<F>> = vec![];

    fib_seq.push(initial_values[0].clone());
    fib_seq.push(initial_values[1].clone());

    for i in 2..(trace_length) {
        fib_seq.push(fib_seq[i - 1].clone() + fib_seq[i - 2].clone());
    }

    let last_value = fib_seq[trace_length - 1].clone();
    let mut fib_permuted = fib_seq.clone();
    fib_permuted[0] = last_value;
    fib_permuted[trace_length - 1] = initial_values[0].clone();

    fib_seq.push(FieldElement::<F>::zero());
    fib_permuted.push(FieldElement::<F>::zero());
    let mut trace_cols = vec![fib_seq, fib_permuted];
    resize_to_next_power_of_two(&mut trace_cols);

    TraceTable::from_columns(trace_cols, 2, 1)
}

#[cfg(test)]
mod test {
    use super::*;
    use lambdaworks_math::field::fields::u64_prime_field::FE17;

    #[test]
    fn test_build_fibonacci_rap_trace() {
        // The fibonacci RAP trace should have two columns:
        //     * The usual fibonacci sequence column
        //     * The permuted fibonacci sequence column. The first and last elements are permuted.
        // Also, a 0 is appended at the end of both columns. The reason for this can be read in
        // https://hackmd.io/@aztec-network/plonk-arithmetiization-air#RAPs---PAIRs-with-interjected-verifier-randomness

        let trace = fibonacci_rap_trace([FE17::from(1), FE17::from(1)], 8);
        let mut expected_trace = vec![
            vec![
                FE17::one(),
                FE17::one(),
                FE17::from(2),
                FE17::from(3),
                FE17::from(5),
                FE17::from(8),
                FE17::from(13),
                FE17::from(21),
                FE17::zero(),
            ],
            vec![
                FE17::from(21),
                FE17::one(),
                FE17::from(2),
                FE17::from(3),
                FE17::from(5),
                FE17::from(8),
                FE17::from(13),
                FE17::one(),
                FE17::zero(),
            ],
        ];
        resize_to_next_power_of_two(&mut expected_trace);

        assert_eq!(trace.columns(), expected_trace);
    }

    #[test]
    fn aux_col() {
        let trace = fibonacci_rap_trace([FE17::from(1), FE17::from(1)], 64);
        let trace_cols = trace.columns();

        let not_perm = trace_cols[0].clone();
        let perm = trace_cols[1].clone();
        let gamma = FE17::from(10);

        assert_eq!(perm.len(), not_perm.len());
        let trace_len = not_perm.len();

        let mut aux_col = Vec::new();
        for i in 0..trace_len {
            if i == 0 {
                aux_col.push(FE17::one());
            } else {
                let z_i = aux_col[i - 1];
                let n_p_term = not_perm[i - 1] + gamma;
                let p_term = perm[i - 1] + gamma;

                aux_col.push(z_i * n_p_term.div(p_term));
            }
        }

        assert_eq!(aux_col.last().unwrap(), &FE17::one());
    }
}
