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
struct ContinuityConstraint<F: IsFFTField> {
    phantom: PhantomData<F>,
}

impl<F: IsFFTField> ContinuityConstraint<F> {
    pub fn new() -> Self {
        Self {
            phantom: PhantomData,
        }
    }
}

impl<F> TransitionConstraint<F, F> for ContinuityConstraint<F>
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
        // NOTE: We are assuming that hte trace has as length a power of 2.
        1
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

        let a_perm0 = first_step.get_main_evaluation_element(0, 2);
        let a_perm1 = second_step.get_main_evaluation_element(0, 2);
        let res = (a_perm1 - a_perm0) * (a_perm1 - a_perm0 - FieldElement::<F>::one());

        transition_evaluations[self.constraint_idx()] = res;
    }
}

#[derive(Clone)]
struct SingleValueConstraint<F: IsFFTField> {
    phantom: PhantomData<F>,
}

impl<F: IsFFTField> SingleValueConstraint<F> {
    pub fn new() -> Self {
        Self {
            phantom: PhantomData,
        }
    }
}

impl<F> TransitionConstraint<F, F> for SingleValueConstraint<F>
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
        // NOTE: We are assuming that hte trace has as length a power of 2.
        1
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

        let a_perm0 = first_step.get_main_evaluation_element(0, 2);
        let a_perm1 = second_step.get_main_evaluation_element(0, 2);
        let v_perm0 = first_step.get_main_evaluation_element(0, 3);
        let v_perm1 = second_step.get_main_evaluation_element(0, 3);

        let res = (v_perm1 - v_perm0) * (a_perm1 - a_perm0 - FieldElement::<F>::one());

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
        2
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
        let p0 = first_step.get_aux_evaluation_element(0, 0);
        let p1 = second_step.get_aux_evaluation_element(0, 0);
        let alpha = &rap_challenges[0];
        let z = &rap_challenges[1];
        let a1 = second_step.get_main_evaluation_element(0, 0);
        let v1 = second_step.get_main_evaluation_element(0, 1);
        let a_perm_1 = second_step.get_main_evaluation_element(0, 2);
        let v_perm_1 = second_step.get_main_evaluation_element(0, 3);

        let res = (z - (a_perm_1 + alpha * v_perm_1)) * p1 - (z - (a1 + alpha * v1)) * p0;

        transition_evaluations[self.constraint_idx()] = res;
    }
}

pub struct ReadOnlyRAP<F>
where
    F: IsFFTField,
{
    context: AirContext,
    trace_length: usize,
    transition_constraints: Vec<Box<dyn TransitionConstraint<F, F>>>,
}

impl<F> AIR for ReadOnlyRAP<F>
where
    F: IsFFTField + Send + Sync + 'static,
    FieldElement<F>: ByteConversion,
{
    type Field = F;
    type FieldExtension = F;
    type PublicInputs = ();

    const STEP_SIZE: usize = 1;

    fn new(
        trace_length: usize,
        pub_inputs: &Self::PublicInputs,
        proof_options: &ProofOptions,
    ) -> Self {
        let transition_constraints: Vec<
            Box<dyn TransitionConstraint<Self::Field, Self::FieldExtension>>,
        > = vec![
            Box::new(ContinuityConstraint::new()),
            Box::new(PermutationConstraint::new()),
        ];

        let context = AirContext {
            proof_options: proof_options.clone(),
            trace_columns: 5,
            transition_offsets: vec![0, 1],
            transition_exemptions: vec![1],
            num_transition_constraints: transition_constraints.len(),
        };

        Self {
            context,
            trace_length,
            transition_constraints,
        }
    }

    fn build_auxiliary_trace(
        &self,
        main_trace: &TraceTable<Self::Field>,
        challenges: &[FieldElement<F>],
    ) -> TraceTable<Self::Field> {
        let main_segment_cols = main_trace.columns();
        let a = &main_segment_cols[0];
        let v = &main_segment_cols[1];
        let a_perm = &main_segment_cols[2];
        let v_perm = &main_segment_cols[3];
        let z = &challenges[0];
        let alpha = &challenges[1];

        let trace_len = main_trace.n_rows();

        let mut aux_col = Vec::new();
        let num = z - (&a[0] + alpha * &v[0]);
        let den = z - (&a_perm[0] + alpha * &v_perm[0]);
        aux_col.push(num / den);

        for i in 0..trace_len - 1 {
            let num = (z - (&a[i + 1] + alpha * &v[i + 1])) * &aux_col[i];
            let den = z - (&a_perm[i + 1] + alpha * &v_perm[i + 1]);
            aux_col.push(num / den);
        }

        TraceTable::from_columns(vec![aux_col], 0, 1)
    }

    fn build_rap_challenges(
        &self,
        transcript: &mut impl IsTranscript<Self::Field>,
    ) -> Vec<FieldElement<Self::FieldExtension>> {
        vec![
            transcript.sample_field_element(),
            transcript.sample_field_element(),
        ]
    }

    fn trace_layout(&self) -> (usize, usize) {
        (4, 1)
    }

    fn boundary_constraints(
        &self,
        _rap_challenges: &[FieldElement<Self::FieldExtension>],
    ) -> BoundaryConstraints<Self::FieldExtension> {
        // Auxiliary boundary constraints
        let a0_aux = BoundaryConstraint::new_aux(0, 0, FieldElement::<Self::FieldExtension>::one());

        BoundaryConstraints::from_constraints(vec![a0_aux])
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
        &()
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
