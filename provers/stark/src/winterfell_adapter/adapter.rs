use crate::trace::TraceTable as LambdaTraceTable;
use crate::{
    constraints::boundary::{BoundaryConstraint, BoundaryConstraints},
    traits::AIR,
};
use lambdaworks_math::field::{
    element::FieldElement, fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
};
use winterfell::{Air, EvaluationFrame, FieldExtension, ProofOptions, TraceInfo, TraceTable};

#[derive(Clone)]
pub struct AirAdapterPublicInputs<A>
where
    A: Air<BaseField = FieldElement<Stark252PrimeField>>,
    A::PublicInputs: Clone,
{
    winterfell_public_inputs: A::PublicInputs,
    transition_degrees: Vec<usize>,
    transition_exemptions: Vec<usize>,
    transition_offsets: Vec<usize>,
    trace_columns: usize,
    composition_poly_degree_bound: usize,
}

#[derive(Clone)]
pub struct AirAdapter<A>
where
    A: Air<BaseField = FieldElement<Stark252PrimeField>>,
    A::PublicInputs: Clone,
{
    winterfell_air: A,
    public_inputs: AirAdapterPublicInputs<A>,
    air_context: crate::context::AirContext,
}

impl<A> AirAdapter<A>
where
    A: Air<BaseField = FieldElement<Stark252PrimeField>> + Clone,
    A::PublicInputs: Clone,
{
    pub fn convert_winterfell_trace_table(
        trace: TraceTable<FieldElement<Stark252PrimeField>>,
    ) -> LambdaTraceTable<Stark252PrimeField> {
        let mut columns = Vec::new();
        for i in 0..trace.width() {
            columns.push(trace.get_column(i).to_owned());
        }

        LambdaTraceTable::from_columns(&columns)
    }
}

impl<A> AIR for AirAdapter<A>
where
    A: Air<BaseField = FieldElement<Stark252PrimeField>> + Clone,
    A::PublicInputs: Clone,
{
    type Field = Stark252PrimeField;
    type RAPChallenges = (); // RAP Challenges not supported?
    type PublicInputs = AirAdapterPublicInputs<A>;

    fn new(
        trace_length: usize,
        pub_inputs: &Self::PublicInputs,
        lambda_proof_options: &crate::proof::options::ProofOptions,
    ) -> Self {
        let winter_trace_info = TraceInfo::new(pub_inputs.trace_columns, trace_length);
        let trace_columns = winter_trace_info.layout().main_trace_width() + winter_trace_info.layout().aux_trace_width();
        let winter_proof_options = ProofOptions::new(
            lambda_proof_options.fri_number_of_queries,
            lambda_proof_options.blowup_factor as usize,
            lambda_proof_options.grinding_factor as u32,
            FieldExtension::None,
            2,
            0,
        );

        let winterfell_air = A::new(
            winter_trace_info,
            pub_inputs.winterfell_public_inputs.clone(),
            winter_proof_options,
        );
        let winterfell_context = winterfell_air.context();

        let lambda_context = crate::context::AirContext {
            proof_options: lambda_proof_options.clone(),
            transition_degrees: pub_inputs.transition_degrees.to_owned(),
            transition_exemptions: pub_inputs.transition_exemptions.to_owned(),
            transition_offsets: pub_inputs.transition_offsets.to_owned(),
            num_transition_constraints: winterfell_context.num_transition_constraints(),
            trace_columns,
            num_transition_exemptions: winterfell_context.num_transition_exemptions(),
        };

        Self {
            winterfell_air,
            public_inputs: pub_inputs.clone(),
            air_context: lambda_context,
        }
    }

    fn build_auxiliary_trace(
        &self,
        _main_trace: &crate::trace::TraceTable<Self::Field>,
        _rap_challenges: &Self::RAPChallenges,
    ) -> crate::trace::TraceTable<Self::Field> {
        // Not supported
        crate::trace::TraceTable::empty()
    }

    fn build_rap_challenges(
        &self,
        _transcript: &mut impl crate::transcript::IsStarkTranscript<Self::Field>,
    ) -> Self::RAPChallenges {
        // Not supported
    }

    fn number_auxiliary_rap_columns(&self) -> usize {
        // Not supported
        0
    }

    fn composition_poly_degree_bound(&self) -> usize {
        self.public_inputs.composition_poly_degree_bound
    }

    fn compute_transition(
        &self,
        frame: &crate::frame::Frame<Self::Field>,
        _rap_challenges: &Self::RAPChallenges,
    ) -> Vec<FieldElement<Self::Field>> {
        let frame =
            EvaluationFrame::from_rows(frame.get_row(0).to_vec(), frame.get_row(1).to_vec());
        let mut result = vec![FieldElement::zero(); self.num_transition_constraints()];
        self.winterfell_air
            .evaluate_transition::<FieldElement<Stark252PrimeField>>(&frame, &[], &mut result); // Periodic values not supported
        result
    }

    fn boundary_constraints(
        &self,
        _rap_challenges: &Self::RAPChallenges,
    ) -> crate::constraints::boundary::BoundaryConstraints<Self::Field> {
        let mut result = Vec::new();
        for assertion in self.winterfell_air.get_assertions() {
            assert!(assertion.is_single());
            result.push(BoundaryConstraint::new(
                assertion.column(),
                assertion.first_step(),
                assertion.values()[0],
            ));
        }
        BoundaryConstraints::from_constraints(result)
    }

    fn context(&self) -> &crate::context::AirContext {
        &self.air_context
    }

    fn trace_length(&self) -> usize {
        self.winterfell_air.context().trace_len()
    }

    fn pub_inputs(&self) -> &Self::PublicInputs {
        &self.public_inputs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        proof::options::ProofOptions,
        prover::{IsStarkProver, Prover},
        transcript::StoneProverTranscript,
        verifier::{IsStarkVerifier, Verifier},
        winterfell_adapter::examples::{fibonacci_2_terms::{FibAir2Terms, self}, fibonacci_rap::{FibonacciRAP, self}},
    };

    #[test]
    fn prove_and_verify_a_winterfell_fibonacci_2_terms_air() {
        let lambda_proof_options = ProofOptions::default_test_options();
        let trace = AirAdapter::<FibAir2Terms>::convert_winterfell_trace_table(fibonacci_2_terms::build_trace(16));
        let pub_inputs = AirAdapterPublicInputs {
            winterfell_public_inputs: trace.columns()[1][7],
            transition_degrees: vec![1, 1],
            transition_exemptions: vec![1, 1],
            transition_offsets: vec![0, 1],
            composition_poly_degree_bound: 8,
            trace_columns: 2
        };

        let proof = Prover::prove::<AirAdapter<FibAir2Terms>>(
            &trace,
            &pub_inputs,
            &lambda_proof_options,
            StoneProverTranscript::new(&[]),
        )
        .unwrap();
        assert!(Verifier::verify::<AirAdapter<FibAir2Terms>>(
            &proof,
            &pub_inputs,
            &lambda_proof_options,
            StoneProverTranscript::new(&[]),
        ));
    }


    #[test]
    fn prove_and_verify_a_winterfell_fibonacci_rap_air() {
        let lambda_proof_options = ProofOptions::default_test_options();
        let trace = AirAdapter::<FibonacciRAP>::convert_winterfell_trace_table(fibonacci_rap::build_trace(16));
        let pub_inputs = AirAdapterPublicInputs {
            winterfell_public_inputs: trace.columns()[1][15],
            transition_degrees: vec![1, 1],
            transition_exemptions: vec![1, 1],
            transition_offsets: vec![0, 1],
            composition_poly_degree_bound: 16,
            trace_columns: 2
        };

        let proof = Prover::prove::<AirAdapter<FibonacciRAP>>(
            &trace,
            &pub_inputs,
            &lambda_proof_options,
            StoneProverTranscript::new(&[]),
        )
        .unwrap();
        assert!(Verifier::verify::<AirAdapter<FibonacciRAP>>(
            &proof,
            &pub_inputs,
            &lambda_proof_options,
            StoneProverTranscript::new(&[]),
        ));
    }
}
