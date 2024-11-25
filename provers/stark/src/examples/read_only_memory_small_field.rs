use std::marker::PhantomData;

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
use lambdaworks_math::field::{
    fields::fft_friendly::{
        babybear::Babybear31PrimeField, quartic_babybear::Degree4BabyBearExtensionField,
    },
    traits::{IsField, IsPrimeField, IsSubFieldOf},
};
use lambdaworks_math::{
    field::{element::FieldElement, traits::IsFFTField},
    traits::ByteConversion,
};
type F = Babybear31PrimeField;
type E = Degree4BabyBearExtensionField;

#[derive(Clone)]
struct ContinuityConstraint;

impl ContinuityConstraint {
    pub fn new() -> Self {
        Self {}
    }
}

impl TransitionConstraint<F, E> for ContinuityConstraint {
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

    fn evaluate_prover(
        &self,
        frame: &Frame<F, E>,
        transition_evaluations: &mut [FieldElement<E>],
        _periodic_values: &[FieldElement<F>],
        _rap_challenges: &[FieldElement<E>],
    ) {
        let first_step = frame.get_evaluation_step(0);
        let second_step = frame.get_evaluation_step(1);

        let a_perm0 = first_step.get_main_evaluation_element(0, 2);
        let a_perm1 = second_step.get_main_evaluation_element(0, 2);
        let res = (a_perm1 - a_perm0) * (a_perm1 - a_perm0 - FieldElement::<E>::one());

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn evaluate_verifier(
        &self,
        frame: &Frame<E, E>,
        transition_evaluations: &mut [FieldElement<E>],
        _periodic_values: &[FieldElement<E>],
        _rap_challenges: &[FieldElement<E>],
    ) {
        let first_step = frame.get_evaluation_step(0);
        let second_step = frame.get_evaluation_step(1);

        let a_perm0 = first_step.get_main_evaluation_element(0, 2);
        let a_perm1 = second_step.get_main_evaluation_element(0, 2);
        let res = (a_perm1 - a_perm0) * (a_perm1 - a_perm0 - FieldElement::<E>::one());

        transition_evaluations[self.constraint_idx()] = res;
    }
}

#[derive(Clone)]
struct SingleValueConstraint {}

impl SingleValueConstraint {
    pub fn new() -> Self {
        Self {}
    }
}

impl TransitionConstraint<F, E> for SingleValueConstraint {
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

    fn evaluate_prover(
        &self,
        frame: &Frame<F, E>,
        transition_evaluations: &mut [FieldElement<E>],
        _periodic_values: &[FieldElement<F>],
        _rap_challenges: &[FieldElement<E>],
    ) {
        let first_step = frame.get_evaluation_step(0);
        let second_step = frame.get_evaluation_step(1);

        let a_perm0 = first_step.get_main_evaluation_element(0, 2);
        let a_perm1 = second_step.get_main_evaluation_element(0, 2);
        let v_perm0 = first_step.get_main_evaluation_element(0, 3);
        let v_perm1 = second_step.get_main_evaluation_element(0, 3);

        let res = (v_perm1 - v_perm0) * (a_perm1 - a_perm0 - FieldElement::<E>::one());

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn evaluate_verifier(
        &self,
        frame: &Frame<E, E>,
        transition_evaluations: &mut [FieldElement<E>],
        _periodic_values: &[FieldElement<E>],
        _rap_challenges: &[FieldElement<E>],
    ) {
        let first_step = frame.get_evaluation_step(0);
        let second_step = frame.get_evaluation_step(1);

        let a_perm0 = first_step.get_main_evaluation_element(0, 2);
        let a_perm1 = second_step.get_main_evaluation_element(0, 2);
        let v_perm0 = first_step.get_main_evaluation_element(0, 3);
        let v_perm1 = second_step.get_main_evaluation_element(0, 3);

        let res = (v_perm1 - v_perm0) * (a_perm1 - a_perm0 - FieldElement::<E>::one());

        transition_evaluations[self.constraint_idx()] = res;
    }
}

#[derive(Clone)]
struct PermutationConstraint {}

impl PermutationConstraint {
    pub fn new() -> Self {
        Self {}
    }
}

impl TransitionConstraint<F, E> for PermutationConstraint {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        2
    }

    fn end_exemptions(&self) -> usize {
        1
    }

    fn evaluate_prover(
        &self,
        frame: &Frame<F, E>,
        transition_evaluations: &mut [FieldElement<E>],
        _periodic_values: &[FieldElement<F>],
        rap_challenges: &[FieldElement<E>],
    ) {
        let first_step = frame.get_evaluation_step(0);
        let second_step = frame.get_evaluation_step(1);

        // Auxiliary constraints
        let p0 = first_step.get_aux_evaluation_element(0, 0);
        let p1 = second_step.get_aux_evaluation_element(0, 0);
        let z = &rap_challenges[0];
        let alpha = &rap_challenges[1];
        let a1 = second_step.get_main_evaluation_element(0, 0);
        let v1 = second_step.get_main_evaluation_element(0, 1);
        let a_perm_1 = second_step.get_main_evaluation_element(0, 2);
        let v_perm_1 = second_step.get_main_evaluation_element(0, 3);

        let res = (z - (a_perm_1 + v_perm_1 * alpha)) * p1 - (z - (a1 + v1 * alpha)) * p0;

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn evaluate_verifier(
        &self,
        frame: &Frame<E, E>,
        transition_evaluations: &mut [FieldElement<E>],
        _periodic_values: &[FieldElement<E>],
        rap_challenges: &[FieldElement<E>],
    ) {
        let first_step = frame.get_evaluation_step(0);
        let second_step = frame.get_evaluation_step(1);

        // Auxiliary constraints
        let p0 = first_step.get_aux_evaluation_element(0, 0);
        let p1 = second_step.get_aux_evaluation_element(0, 0);
        let z = &rap_challenges[0];
        let alpha = &rap_challenges[1];
        let a1 = second_step.get_main_evaluation_element(0, 0);
        let v1 = second_step.get_main_evaluation_element(0, 1);
        let a_perm_1 = second_step.get_main_evaluation_element(0, 2);
        let v_perm_1 = second_step.get_main_evaluation_element(0, 3);

        let res = (z - (a_perm_1 + v_perm_1 * alpha)) * p1 - (z - (a1 + alpha * v1)) * p0;

        transition_evaluations[self.constraint_idx()] = res;
    }
}

pub struct ReadOnlyRAP {
    context: AirContext,
    trace_length: usize,
    pub_inputs: ReadOnlyPublicInputs,
    transition_constraints: Vec<Box<dyn TransitionConstraint<F, E>>>,
}

#[derive(Clone, Debug)]
pub struct ReadOnlyPublicInputs {
    pub a0: FieldElement<F>,
    pub v0: FieldElement<F>,
    pub a_perm0: FieldElement<F>,
    pub v_perm0: FieldElement<F>,
}

impl AIR for ReadOnlyRAP {
    type Field = F;
    type FieldExtension = E;
    type PublicInputs = ReadOnlyPublicInputs;

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
        challenges: &[FieldElement<E>],
    ) {
        let main_segment_cols = trace.columns_main();
        let a = &main_segment_cols[0];
        let v = &main_segment_cols[1];
        let a_perm = &main_segment_cols[2];
        let v_perm = &main_segment_cols[3];
        let z = &challenges[0];
        let alpha = &challenges[1];

        let trace_len = trace.num_rows();

        let mut aux_col = Vec::new();
        let num = z - (&a[0] + &v[0] * alpha);
        let den = z - (&a_perm[0] + &v_perm[0] * alpha);
        aux_col.push(num / den);

        for i in 0..trace_len - 1 {
            let num = (z - (&a[i + 1] + &v[i + 1]) * alpha) * &aux_col[i];
            let den = z - (&a_perm[i + 1] + &v_perm[i + 1]) * alpha;
            aux_col.push(num / den);
        }

        for (i, aux_elem) in aux_col.iter().enumerate().take(trace.num_rows()) {
            trace.set_aux(i, 0, aux_elem.clone())
        }
    }

    fn build_rap_challenges(
        &self,
        transcript: &mut impl IsTranscript<Self::FieldExtension>,
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
        let a_perm0 = &self.pub_inputs.a_perm0;
        let v_perm0 = &self.pub_inputs.v_perm0;
        let z = &rap_challenges[0];
        let alpha = &rap_challenges[1];

        // TODO: We shouldn't have to embedd tha main boundary constraints into the extension field.
        // Main boundary constraints
        let c1 = BoundaryConstraint::new_main(0, 0, a0.clone().to_extension());
        let c2 = BoundaryConstraint::new_main(1, 0, v0.clone().to_extension());
        let c3 = BoundaryConstraint::new_main(2, 0, a_perm0.clone().to_extension());
        let c4 = BoundaryConstraint::new_main(3, 0, v_perm0.clone().to_extension());

        // Auxiliary boundary constraints
        let num = -(a0 + v0 * alpha) + z;
        let den = -(a_perm0 + v_perm0 * alpha) + z;
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

pub fn sort_rap_trace(
    address: Vec<FieldElement<F>>,
    value: Vec<FieldElement<F>>,
) -> TraceTable<F, E> {
    let mut address_value_pairs: Vec<_> = address.iter().zip(value.iter()).collect();

    address_value_pairs.sort_by_key(|(addr, _)| addr.representative());

    let (sorted_address, sorted_value): (Vec<FieldElement<F>>, Vec<FieldElement<F>>) =
        address_value_pairs
            .into_iter()
            .map(|(addr, val)| (addr.clone(), val.clone()))
            .unzip();
    let main_columns = vec![address.clone(), value.clone(), sorted_address, sorted_value];
    // create a vector with zeros of the same length as the main columns
    let zero_vec = vec![FieldElement::<E>::zero(); main_columns[0].len()];
    TraceTable::from_columns(main_columns, vec![zero_vec], 1)
}
#[cfg(test)]
mod test {
    use super::*;
    use lambdaworks_math::field::fields::fft_friendly::babybear::Babybear31PrimeField;
    use lambdaworks_math::field::fields::fft_friendly::quadratic_babybear::QuadraticBabybearField;

    type FE = FieldElement<Babybear31PrimeField>;
    type QFE = FieldElement<QuadraticBabybearField>;

    #[test]
    fn test_sort_rap_trace() {
        let address_col = vec![
            FE::from(5),
            FE::from(2),
            FE::from(3),
            FE::from(4),
            FE::from(1),
            FE::from(6),
            FE::from(7),
            FE::from(8),
        ];
        let value_col = vec![
            FE::from(50),
            FE::from(20),
            FE::from(30),
            FE::from(40),
            FE::from(10),
            FE::from(60),
            FE::from(70),
            FE::from(80),
        ];

        let sorted_trace = sort_rap_trace(address_col.clone(), value_col.clone());

        let expected_sorted_addresses = vec![
            FE::from(1),
            FE::from(2),
            FE::from(3),
            FE::from(4),
            FE::from(5),
            FE::from(6),
            FE::from(7),
            FE::from(8),
        ];
        let expected_sorted_values = vec![
            FE::from(10),
            FE::from(20),
            FE::from(30),
            FE::from(40),
            FE::from(50),
            FE::from(60),
            FE::from(70),
            FE::from(80),
        ];

        assert_eq!(sorted_trace.columns_main()[2], expected_sorted_addresses);
        assert_eq!(sorted_trace.columns_main()[3], expected_sorted_values);
    }
}
