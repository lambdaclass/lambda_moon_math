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
use lambdaworks_crypto::fiat_shamir::is_transcript::IsTranscript;
use lambdaworks_math::field::traits::IsPrimeField;
use lambdaworks_math::{
    field::{element::FieldElement, traits::IsFFTField},
    traits::ByteConversion,
};

/// This condition ensures the continuity in a read-only memory structure, preserving strict ordering.
/// Equation based on Cairo Whitepaper section 9.7.2
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
        // NOTE: We are assuming that the trace has as length a power of 2.
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

        let a_sorted_0 = first_step.get_main_evaluation_element(0, 2);
        let a_sorted_1 = second_step.get_main_evaluation_element(0, 2);
        // (a'_{i+1} - a'_i)(a'_{i+1} - a'_i - 1) = 0 where a' is the sorted address
        let res = (a_sorted_1 - a_sorted_0) * (a_sorted_1 - a_sorted_0 - FieldElement::<F>::one());

        // The eval always exists, except if the constraint idx were incorrectly defined.
        if let Some(eval) = transition_evaluations.get_mut(self.constraint_idx()) {
            *eval = res;
        }
    }
}
/// Transition constraint that ensures that same addresses have same values, making the memory read-only.
/// Equation based on Cairo Whitepaper section 9.7.2
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
        // NOTE: We are assuming that the trace has as length a power of 2.
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

        let a_sorted0 = first_step.get_main_evaluation_element(0, 2);
        let a_sorted1 = second_step.get_main_evaluation_element(0, 2);
        let v_sorted0 = first_step.get_main_evaluation_element(0, 3);
        let v_sorted1 = second_step.get_main_evaluation_element(0, 3);
        // (v'_{i+1} - v'_i) * (a'_{i+1} - a'_i - 1) = 0
        let res = (v_sorted1 - v_sorted0) * (a_sorted1 - a_sorted0 - FieldElement::<F>::one());

        // The eval always exists, except if the constraint idx were incorrectly defined.
        if let Some(eval) = transition_evaluations.get_mut(self.constraint_idx()) {
            *eval = res;
        }
    }
}
/// Permutation constraint ensures that the values are permuted in the memory.
/// Equation based on Cairo Whitepaper section 9.7.2
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
        evaluation_context: &TransitionEvaluationContext<F, F>,
        transition_evaluations: &mut [FieldElement<F>],
    ) {
        let (frame, _periodic_values, rap_challenges) = match evaluation_context {
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

        // Auxiliary constraints
        let p0 = first_step.get_aux_evaluation_element(0, 0);
        let p1 = second_step.get_aux_evaluation_element(0, 0);
        let z = &rap_challenges[0];
        let alpha = &rap_challenges[1];
        let a1 = second_step.get_main_evaluation_element(0, 0);
        let v1 = second_step.get_main_evaluation_element(0, 1);
        let a_sorted_1 = second_step.get_main_evaluation_element(0, 2);
        let v_sorted_1 = second_step.get_main_evaluation_element(0, 3);
        // (z - (a'_{i+1} + α * v'_{i+1})) * p_{i+1} = (z - (a_{i+1} + α * v_{i+1})) * p_i
        let res = (z - (a_sorted_1 + alpha * v_sorted_1)) * p1 - (z - (a1 + alpha * v1)) * p0;

        // The eval always exists, except if the constraint idx were incorrectly defined.
        if let Some(eval) = transition_evaluations.get_mut(self.constraint_idx()) {
            *eval = res;
        }
    }
}

pub struct ReadOnlyRAP<F>
where
    F: IsFFTField,
{
    context: AirContext,
    trace_length: usize,
    pub_inputs: ReadOnlyPublicInputs<F>,
    transition_constraints: Vec<Box<dyn TransitionConstraint<F, F>>>,
}

#[derive(Clone, Debug)]
pub struct ReadOnlyPublicInputs<F>
where
    F: IsFFTField,
{
    pub a0: FieldElement<F>,
    pub v0: FieldElement<F>,
    pub a_sorted0: FieldElement<F>,
    pub v_sorted0: FieldElement<F>,
}

impl<F> AIR for ReadOnlyRAP<F>
where
    F: IsFFTField + Send + Sync + 'static,
    FieldElement<F>: ByteConversion,
{
    type Field = F;
    type FieldExtension = F;
    type PublicInputs = ReadOnlyPublicInputs<F>;

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
            Box::new(SingleValueConstraint::new()),
            Box::new(PermutationConstraint::new()),
        ];

        let context = AirContext {
            proof_options: proof_options.clone(),
            trace_columns: 5,
            transition_offsets: vec![0, 1],
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
        trace: &mut TraceTable<Self::Field, Self::FieldExtension>,
        challenges: &[FieldElement<F>],
    ) {
        let main_segment_cols = trace.columns_main();
        let a = &main_segment_cols[0];
        let v = &main_segment_cols[1];
        let a_sorted = &main_segment_cols[2];
        let v_sorted = &main_segment_cols[3];
        let z = &challenges[0];
        let alpha = &challenges[1];

        let trace_len = trace.num_rows();

        let mut aux_col = Vec::new();
        let num = z - (&a[0] + alpha * &v[0]);
        let den = z - (&a_sorted[0] + alpha * &v_sorted[0]);
        aux_col.push(num / den);
        // Apply the same equation given in the permutation case to the rest of the trace
        for i in 0..trace_len - 1 {
            let num = (z - (&a[i + 1] + alpha * &v[i + 1])) * &aux_col[i];
            let den = z - (&a_sorted[i + 1] + alpha * &v_sorted[i + 1]);
            aux_col.push(num / den);
        }

        for (i, aux_elem) in aux_col.iter().enumerate().take(trace.num_rows()) {
            trace.set_aux(i, 0, aux_elem.clone())
        }
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
        rap_challenges: &[FieldElement<Self::FieldExtension>],
    ) -> BoundaryConstraints<Self::FieldExtension> {
        let a0 = &self.pub_inputs.a0;
        let v0 = &self.pub_inputs.v0;
        let a_sorted0 = &self.pub_inputs.a_sorted0;
        let v_sorted0 = &self.pub_inputs.v_sorted0;
        let z = &rap_challenges[0];
        let alpha = &rap_challenges[1];

        // Main boundary constraints
        let c1 = BoundaryConstraint::new_main(0, 0, a0.clone());
        let c2 = BoundaryConstraint::new_main(1, 0, v0.clone());
        let c3 = BoundaryConstraint::new_main(2, 0, a_sorted0.clone());
        let c4 = BoundaryConstraint::new_main(3, 0, v_sorted0.clone());

        // Auxiliary boundary constraints
        let num = z - (a0 + alpha * v0);
        let den = z - (a_sorted0 + alpha * v_sorted0);
        let p0_value = num / den;

        let c_aux1 = BoundaryConstraint::new_aux(0, 0, p0_value);
        let c_aux2 = BoundaryConstraint::new_aux(
            0,
            self.trace_length - 1,
            FieldElement::<Self::FieldExtension>::one(),
        );

        BoundaryConstraints::from_constraints(vec![c1, c2, c3, c4, c_aux1, c_aux2])
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
}

/// Given the adress and value columns, it returns the trace table with 5 columns, which are:
/// Addres, Value, Adress Sorted, Value Sorted and a Column of Zeroes (where we'll insert the auxiliary colunn).
pub fn sort_rap_trace<F: IsFFTField + IsPrimeField>(
    address: Vec<FieldElement<F>>,
    value: Vec<FieldElement<F>>,
) -> TraceTable<F, F> {
    let mut address_value_pairs: Vec<_> = address.iter().zip(value.iter()).collect();

    address_value_pairs.sort_by_key(|(addr, _)| addr.representative());

    let (sorted_address, sorted_value): (Vec<FieldElement<F>>, Vec<FieldElement<F>>) =
        address_value_pairs
            .into_iter()
            .map(|(addr, val)| (addr.clone(), val.clone()))
            .unzip();
    let main_columns = vec![address.clone(), value.clone(), sorted_address, sorted_value];
    // create a vector with zeros of the same length as the main columns
    let zero_vec = vec![FieldElement::<F>::zero(); main_columns[0].len()];
    TraceTable::from_columns(main_columns, vec![zero_vec], 1)
}

#[cfg(test)]
mod test {
    use super::*;
    use lambdaworks_math::field::fields::u64_prime_field::FE17;

    #[test]
    fn test_sort_rap_trace() {
        let address_col = vec![
            FE17::from(5),
            FE17::from(2),
            FE17::from(3),
            FE17::from(4),
            FE17::from(1),
            FE17::from(6),
            FE17::from(7),
            FE17::from(8),
        ];
        let value_col = vec![
            FE17::from(50),
            FE17::from(20),
            FE17::from(30),
            FE17::from(40),
            FE17::from(10),
            FE17::from(60),
            FE17::from(70),
            FE17::from(80),
        ];

        let sorted_trace = sort_rap_trace(address_col.clone(), value_col.clone());

        let expected_sorted_addresses = vec![
            FE17::from(1),
            FE17::from(2),
            FE17::from(3),
            FE17::from(4),
            FE17::from(5),
            FE17::from(6),
            FE17::from(7),
            FE17::from(8),
        ];
        let expected_sorted_values = vec![
            FE17::from(10),
            FE17::from(20),
            FE17::from(30),
            FE17::from(40),
            FE17::from(50),
            FE17::from(60),
            FE17::from(70),
            FE17::from(80),
        ];

        assert_eq!(sorted_trace.columns_main()[2], expected_sorted_addresses);
        assert_eq!(sorted_trace.columns_main()[3], expected_sorted_values);
    }
}
