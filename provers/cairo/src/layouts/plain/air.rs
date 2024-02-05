use crate::{
    cairo_mem::CairoMemory,
    execution_trace::{set_mem_permutation_column, set_rc_permutation_column, CairoTraceTable},
    register_states::RegisterStates,
    transition_constraints::*,
};
use cairo_vm::{air_public_input::MemorySegmentAddresses, without_std::collections::HashMap};
use itertools::Itertools;
use lambdaworks_math::{
    errors::DeserializationError,
    field::{
        element::FieldElement, fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
    },
    traits::{AsBytes, ByteConversion, Deserializable},
};
use stark_platinum_prover::constraints::transition::TransitionConstraint;
use stark_platinum_prover::{
    constraints::boundary::{BoundaryConstraint, BoundaryConstraints},
    context::AirContext,
    frame::Frame,
    proof::{options::ProofOptions, stark::StarkProof},
    prover::{IsStarkProver, Prover, ProvingError},
    trace::TraceTable,
    traits::AIR,
    transcript::{IsStarkTranscript, StoneProverTranscript},
    verifier::{IsStarkVerifier, Verifier},
    Felt252,
};

#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum SegmentName {
    RangeCheck,
    Output,
    Program,
    Execution,
    Ecdsa,
    Pedersen,
}

impl From<&str> for SegmentName {
    fn from(value: &str) -> Self {
        match value {
            "range_check" => SegmentName::RangeCheck,
            "output" => SegmentName::Output,
            "program" => SegmentName::Program,
            "execution" => SegmentName::Execution,
            "ecdsa" => SegmentName::Ecdsa,
            "pedersen" => SegmentName::Pedersen,
            n => panic!("Invalid segment name {n}"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Segment {
    pub begin_addr: usize,
    pub stop_ptr: usize,
}

impl Segment {
    pub fn new(begin_addr: u64, stop_ptr: u64) -> Self {
        let begin_addr: usize = begin_addr.try_into().unwrap();
        let stop_ptr: usize = stop_ptr.try_into().unwrap();

        stop_ptr.checked_sub(begin_addr).unwrap();

        Self {
            begin_addr,
            stop_ptr,
        }
    }

    pub fn segment_size(&self) -> usize {
        self.stop_ptr - self.begin_addr - 1
    }
}

impl From<&MemorySegmentAddresses> for Segment {
    fn from(value: &MemorySegmentAddresses) -> Self {
        Self {
            begin_addr: value.begin_addr,
            stop_ptr: value.stop_ptr,
        }
    }
}

pub type MemorySegmentMap = HashMap<SegmentName, Segment>;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PublicInputs {
    pub pc_init: Felt252,
    pub ap_init: Felt252,
    pub fp_init: Felt252,
    pub pc_final: Felt252,
    pub ap_final: Felt252,
    // These are Option because they're not known until
    // the trace is obtained. They represent the minimum
    // and maximum offsets used during program execution.
    // TODO: A possible refactor is moving them to the proof.
    // minimum range check value (0 < range_check_min < range_check_max < 2^16)
    pub range_check_min: Option<u16>,
    // maximum range check value
    pub range_check_max: Option<u16>,
    // Range-check builtin address range
    pub memory_segments: MemorySegmentMap,
    pub public_memory: HashMap<Felt252, Felt252>,
    pub num_steps: usize, // number of execution steps
}

impl PublicInputs {
    /// Creates a Public Input from register states and memory
    /// - In the future we should use the output of the Cairo Runner. This is not currently supported in Cairo RS
    /// - RangeChecks are not filled, and the prover mutates them inside the prove function. This works but also should be loaded from the Cairo RS output
    pub fn from_regs_and_mem(
        register_states: &RegisterStates,
        memory: &CairoMemory,
        codelen: usize,
    ) -> Self {
        let public_memory = (1..=codelen as u64)
            .map(|i| (Felt252::from(i), *memory.get(&i).unwrap()))
            .collect::<HashMap<Felt252, Felt252>>();

        let last_step = &register_states.rows[register_states.steps() - 1];

        PublicInputs {
            pc_init: Felt252::from(register_states.rows[0].pc),
            ap_init: Felt252::from(register_states.rows[0].ap),
            fp_init: Felt252::from(register_states.rows[0].fp),
            pc_final: FieldElement::from(last_step.pc),
            ap_final: FieldElement::from(last_step.ap),
            range_check_min: None,
            range_check_max: None,
            memory_segments: MemorySegmentMap::new(),
            public_memory,
            num_steps: register_states.steps(),
        }
    }
}

impl AsBytes for PublicInputs {
    fn as_bytes(&self) -> Vec<u8> {
        let mut bytes = vec![];
        let pc_init_bytes = self.pc_init.to_bytes_be();
        let felt_length = pc_init_bytes.len();
        bytes.extend(felt_length.to_be_bytes());
        bytes.extend(pc_init_bytes);
        bytes.extend(self.ap_init.to_bytes_be());
        bytes.extend(self.fp_init.to_bytes_be());
        bytes.extend(self.pc_final.to_bytes_be());
        bytes.extend(self.ap_final.to_bytes_be());

        if let Some(range_check_min) = self.range_check_min {
            bytes.extend(1u8.to_be_bytes());
            bytes.extend(range_check_min.to_be_bytes());
        } else {
            bytes.extend(0u8.to_be_bytes());
        }

        if let Some(range_check_max) = self.range_check_max {
            bytes.extend(1u8.to_be_bytes());
            bytes.extend(range_check_max.to_be_bytes());
        } else {
            bytes.extend(0u8.to_be_bytes());
        }

        let mut memory_segment_bytes = vec![];
        for (segment, range) in self.memory_segments.iter() {
            let segment_type = match segment {
                SegmentName::RangeCheck => 0u8,
                SegmentName::Output => 1u8,
                SegmentName::Program => 2u8,
                SegmentName::Execution => 3u8,
                SegmentName::Ecdsa => 4u8,
                SegmentName::Pedersen => 5u8,
            };
            memory_segment_bytes.extend(segment_type.to_be_bytes());
            memory_segment_bytes.extend(range.begin_addr.to_be_bytes());
            memory_segment_bytes.extend(range.stop_ptr.to_be_bytes());
        }
        let memory_segment_length = self.memory_segments.len();
        bytes.extend(memory_segment_length.to_be_bytes());
        bytes.extend(memory_segment_bytes);

        let mut public_memory_bytes = vec![];
        for (address, value) in self.public_memory.iter() {
            public_memory_bytes.extend(address.to_bytes_be());
            public_memory_bytes.extend(value.to_bytes_be());
        }
        let public_memory_length = self.public_memory.len();
        bytes.extend(public_memory_length.to_be_bytes());
        bytes.extend(public_memory_bytes);

        bytes.extend(self.num_steps.to_be_bytes());

        bytes
    }
}

impl Deserializable for PublicInputs {
    fn deserialize(bytes: &[u8]) -> Result<Self, DeserializationError>
    where
        Self: Sized,
    {
        let mut bytes = bytes;
        let felt_len = usize::from_be_bytes(
            bytes
                .get(0..8)
                .ok_or(DeserializationError::InvalidAmountOfBytes)?
                .try_into()
                .map_err(|_| DeserializationError::InvalidAmountOfBytes)?,
        );
        bytes = &bytes[8..];
        let pc_init = Felt252::from_bytes_be(
            bytes
                .get(..felt_len)
                .ok_or(DeserializationError::InvalidAmountOfBytes)?,
        )?;
        bytes = &bytes[felt_len..];
        let ap_init = Felt252::from_bytes_be(
            bytes
                .get(..felt_len)
                .ok_or(DeserializationError::InvalidAmountOfBytes)?,
        )?;
        bytes = &bytes[felt_len..];
        let fp_init = Felt252::from_bytes_be(
            bytes
                .get(..felt_len)
                .ok_or(DeserializationError::InvalidAmountOfBytes)?,
        )?;
        bytes = &bytes[felt_len..];
        let pc_final = Felt252::from_bytes_be(
            bytes
                .get(..felt_len)
                .ok_or(DeserializationError::InvalidAmountOfBytes)?,
        )?;
        bytes = &bytes[felt_len..];
        let ap_final = Felt252::from_bytes_be(
            bytes
                .get(..felt_len)
                .ok_or(DeserializationError::InvalidAmountOfBytes)?,
        )?;
        bytes = &bytes[felt_len..];

        if bytes.is_empty() {
            return Err(DeserializationError::InvalidAmountOfBytes);
        }
        let range_check_min = match bytes[0] {
            0 => {
                bytes = &bytes[1..];
                None
            }
            1 => {
                bytes = &bytes[1..];
                let range_check_min = u16::from_be_bytes(
                    bytes[..2]
                        .try_into()
                        .map_err(|_| DeserializationError::InvalidAmountOfBytes)?,
                );
                bytes = &bytes[2..];
                Some(range_check_min)
            }
            _ => return Err(DeserializationError::FieldFromBytesError),
        };

        if bytes.is_empty() {
            return Err(DeserializationError::InvalidAmountOfBytes);
        }
        let range_check_max = match bytes[0] {
            0 => {
                bytes = &bytes[1..];
                None
            }
            1 => {
                bytes = &bytes[1..];
                let range_check_max = u16::from_be_bytes(
                    bytes[..2]
                        .try_into()
                        .map_err(|_| DeserializationError::InvalidAmountOfBytes)?,
                );
                bytes = &bytes[2..];
                Some(range_check_max)
            }
            _ => return Err(DeserializationError::FieldFromBytesError),
        };

        let mut memory_segments = MemorySegmentMap::new();
        let memory_segment_length = usize::from_be_bytes(
            bytes
                .get(0..8)
                .ok_or(DeserializationError::InvalidAmountOfBytes)?
                .try_into()
                .map_err(|_| DeserializationError::InvalidAmountOfBytes)?,
        );
        bytes = &bytes[8..];
        for _ in 0..memory_segment_length {
            if bytes.is_empty() {
                return Err(DeserializationError::InvalidAmountOfBytes);
            }
            let segment_type = match bytes[0] {
                0u8 => SegmentName::RangeCheck,
                1u8 => SegmentName::Output,
                2u8 => SegmentName::Program,
                3u8 => SegmentName::Execution,
                4u8 => SegmentName::Ecdsa,
                5u8 => SegmentName::Pedersen,
                _ => return Err(DeserializationError::FieldFromBytesError),
            };
            bytes = &bytes[1..];
            let start = u64::from_be_bytes(
                bytes
                    .get(0..8)
                    .ok_or(DeserializationError::InvalidAmountOfBytes)?
                    .try_into()
                    .map_err(|_| DeserializationError::InvalidAmountOfBytes)?,
            );
            bytes = &bytes[8..];
            let end = u64::from_be_bytes(
                bytes
                    .get(0..8)
                    .ok_or(DeserializationError::InvalidAmountOfBytes)?
                    .try_into()
                    .map_err(|_| DeserializationError::InvalidAmountOfBytes)?,
            );
            bytes = &bytes[8..];
            memory_segments.insert(segment_type, Segment::new(start, end));
        }

        let mut public_memory = HashMap::new();
        let public_memory_length = usize::from_be_bytes(
            bytes
                .get(0..8)
                .ok_or(DeserializationError::InvalidAmountOfBytes)?
                .try_into()
                .map_err(|_| DeserializationError::InvalidAmountOfBytes)?,
        );
        bytes = &bytes[8..];
        for _ in 0..public_memory_length {
            let address = Felt252::from_bytes_be(
                bytes
                    .get(..felt_len)
                    .ok_or(DeserializationError::InvalidAmountOfBytes)?,
            )?;
            bytes = &bytes[felt_len..];
            let value = Felt252::from_bytes_be(
                bytes
                    .get(..felt_len)
                    .ok_or(DeserializationError::InvalidAmountOfBytes)?,
            )?;
            bytes = &bytes[felt_len..];
            public_memory.insert(address, value);
        }

        let num_steps = usize::from_be_bytes(
            bytes
                .get(0..8)
                .ok_or(DeserializationError::InvalidAmountOfBytes)?
                .try_into()
                .map_err(|_| DeserializationError::InvalidAmountOfBytes)?,
        );

        Ok(Self {
            pc_init,
            ap_init,
            fp_init,
            pc_final,
            ap_final,
            range_check_min,
            range_check_max,
            memory_segments,
            public_memory,
            num_steps,
        })
    }
}

pub struct CairoAIR {
    pub context: AirContext,
    pub trace_length: usize,
    pub pub_inputs: PublicInputs,
    pub transition_constraints:
        Vec<Box<dyn TransitionConstraint<Stark252PrimeField, Stark252PrimeField>>>,
}

impl AIR for CairoAIR {
    type Field = Stark252PrimeField;
    type FieldExtension = Stark252PrimeField;
    type PublicInputs = PublicInputs;

    const STEP_SIZE: usize = 16;

    /// Creates a new CairoAIR from proof_options
    ///
    /// # Arguments
    ///
    /// * `trace_length` - Length of the Cairo execution trace. Must be a power fo two.
    /// * `pub_inputs` - Public inputs sent by the Cairo runner.
    /// * `proof_options` - STARK proving configuration options.
    fn new(
        trace_length: usize,
        pub_inputs: &Self::PublicInputs,
        proof_options: &ProofOptions,
    ) -> Self {
        debug_assert!(trace_length.is_power_of_two());
        let trace_columns = 8;

        let transition_constraints: Vec<
            Box<dyn TransitionConstraint<Stark252PrimeField, Stark252PrimeField>>,
        > = vec![
            Box::new(BitPrefixFlag::new()),
            Box::new(ZeroFlagConstraint::new()),
            Box::new(InstructionUnpacking::new()),
            Box::new(CpuOperandsMemDstAddr::new()),
            Box::new(CpuOperandsMem0Addr::new()),
            Box::new(CpuOperandsMem1Addr::new()),
            Box::new(CpuUpdateRegistersApUpdate::new()),
            Box::new(CpuUpdateRegistersFpUpdate::new()),
            Box::new(CpuUpdateRegistersPcCondPositive::new()),
            Box::new(CpuUpdateRegistersPcCondNegative::new()),
            Box::new(CpuUpdateRegistersUpdatePcTmp0::new()),
            Box::new(CpuUpdateRegistersUpdatePcTmp1::new()),
            Box::new(CpuOperandsOpsMul::new()),
            Box::new(CpuOperandsRes::new()),
            Box::new(CpuOpcodesCallPushFp::new()),
            Box::new(CpuOpcodesCallPushPc::new()),
            Box::new(CpuOpcodesAssertEq::new()),
            Box::new(MemoryDiffIsBit::new()),
            Box::new(MemoryIsFunc::new()),
            Box::new(MemoryMultiColumnPermStep0::new()),
            Box::new(Rc16DiffIsBit::new()),
            Box::new(Rc16PermStep0::new()),
            Box::new(FlagOp1BaseOp0BitConstraint::new()),
            Box::new(FlagResOp1BitConstraint::new()),
            Box::new(FlagPcUpdateRegularBit::new()),
            Box::new(FlagFpUpdateRegularBit::new()),
            Box::new(CpuOpcodesCallOff0::new()),
            Box::new(CpuOpcodesCallOff1::new()),
            Box::new(CpuOpcodesCallFlags::new()),
            Box::new(CpuOpcodesRetOff0::new()),
            Box::new(CpuOpcodesRetOff2::new()),
            Box::new(CpuOpcodesRetFlags::new()),
        ];

        #[cfg(debug_assertions)]
        {
            use std::collections::HashSet;
            let constraints_set: HashSet<_> = transition_constraints
                .iter()
                .map(|c| c.constraint_idx())
                .collect();
            debug_assert_eq!(
                constraints_set.len(),
                transition_constraints.len(),
                "There are repeated constraint indexes"
            );
            (0..transition_constraints.len())
                .for_each(|idx| debug_assert!(constraints_set.iter().contains(&idx)));

            assert_eq!(transition_constraints.len(), 32);
        }

        let transition_exemptions = transition_constraints
            .iter()
            .map(|c| c.end_exemptions())
            .collect();

        let context = AirContext {
            proof_options: proof_options.clone(),
            trace_columns,
            transition_exemptions,
            transition_offsets: vec![0, 1],
            num_transition_constraints: transition_constraints.len(),
        };

        // The number of the transition constraints
        // and transition exemptions should be the same always.
        debug_assert_eq!(
            context.transition_exemptions.len(),
            context.num_transition_constraints
        );

        Self {
            context,
            pub_inputs: pub_inputs.clone(),
            trace_length,
            transition_constraints,
        }
    }

    fn build_auxiliary_trace(
        &self,
        trace: &mut TraceTable<Self::Field, Self::FieldExtension>,
        rap_challenges: &[Felt252],
    ) {
        let alpha_mem = rap_challenges[0];
        let z_mem = rap_challenges[1];
        let z_rc = rap_challenges[2];

        set_rc_permutation_column(trace, &z_rc);
        set_mem_permutation_column(trace, &alpha_mem, &z_mem);
    }

    fn build_rap_challenges(
        &self,
        transcript: &mut impl IsStarkTranscript<Self::Field>,
    ) -> Vec<Felt252> {
        let alpha_memory = transcript.sample_field_element();
        let z_memory = transcript.sample_field_element();
        let z_rc = transcript.sample_field_element();

        vec![alpha_memory, z_memory, z_rc]
    }

    fn trace_layout(&self) -> (usize, usize) {
        (6, 2)
    }

    /// From the Cairo whitepaper, section 9.10.
    /// These are part of the register constraints.
    ///
    /// Boundary constraints:
    ///  * ap_0 = fp_0 = ap_i
    ///  * ap_t = ap_f
    ///  * pc_0 = pc_i
    ///  * pc_t = pc_f
    fn boundary_constraints(&self, rap_challenges: &[Felt252]) -> BoundaryConstraints<Self::Field> {
        let initial_pc = BoundaryConstraint::new_main(3, 0, self.pub_inputs.pc_init);
        let initial_ap = BoundaryConstraint::new_main(5, 0, self.pub_inputs.ap_init);

        let final_pc = BoundaryConstraint::new_main(
            3,
            self.trace_length - Self::STEP_SIZE,
            self.pub_inputs.pc_final,
        );
        let final_ap = BoundaryConstraint::new_main(
            5,
            self.trace_length - Self::STEP_SIZE,
            self.pub_inputs.ap_final,
        );

        let z_memory = rap_challenges[1];
        let alpha_memory = rap_challenges[0];
        let one: FieldElement<Self::Field> = FieldElement::one();

        let mem_cumul_prod_denominator_no_padding = self
            .pub_inputs
            .public_memory
            .iter()
            .fold(one, |product, (address, value)| {
                product * (z_memory - (address + alpha_memory * value))
            });

        const PUB_MEMORY_ADDR_OFFSET: usize = 8;
        let pad_addr = Felt252::one();
        let pad_value = self.pub_inputs.public_memory.get(&pad_addr).unwrap();
        let val = z_memory - (pad_addr + alpha_memory * pad_value);
        let mem_cumul_prod_denominator_pad = val
            .pow(self.trace_length / PUB_MEMORY_ADDR_OFFSET - self.pub_inputs.public_memory.len());
        let mem_cumul_prod_denominator = (mem_cumul_prod_denominator_no_padding
            * mem_cumul_prod_denominator_pad)
            .inv()
            .unwrap();
        let mem_cumul_prod_final =
            z_memory.pow(self.trace_length / PUB_MEMORY_ADDR_OFFSET) * mem_cumul_prod_denominator;

        let mem_cumul_prod_final_constraint =
            BoundaryConstraint::new_aux(1, self.trace_length - 2, mem_cumul_prod_final);

        let rc_cumul_prod_final_constraint =
            BoundaryConstraint::new_aux(0, self.trace_length - 1, one);

        let rc_min_constraint = BoundaryConstraint::new_main(
            2,
            0,
            FieldElement::from(self.pub_inputs.range_check_min.unwrap() as u64),
        );

        let rc_max_constraint = BoundaryConstraint::new_main(
            2,
            self.trace_length - 1,
            FieldElement::from(self.pub_inputs.range_check_max.unwrap() as u64),
        );

        let constraints = vec![
            initial_pc,
            initial_ap,
            final_pc,
            final_ap,
            mem_cumul_prod_final_constraint,
            rc_cumul_prod_final_constraint,
            rc_min_constraint,
            rc_max_constraint,
        ];

        BoundaryConstraints::from_constraints(constraints)
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
        2 * self.trace_length
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

/// Wrapper function for generating Cairo proofs without the need to specify
/// concrete types.
/// The field is set to Stark252PrimeField and the AIR to CairoAIR.
pub fn generate_cairo_proof(
    trace: &mut CairoTraceTable,
    pub_input: &PublicInputs,
    proof_options: &ProofOptions,
) -> Result<StarkProof<Stark252PrimeField, Stark252PrimeField>, ProvingError> {
    Prover::<CairoAIR>::prove(
        trace,
        pub_input,
        proof_options,
        StoneProverTranscript::new(&[]),
    )
}

/// Wrapper function for verifying Cairo proofs without the need to specify
/// concrete types.
/// The field is set to Stark252PrimeField and the AIR to CairoAIR.
pub fn verify_cairo_proof(
    proof: &StarkProof<Stark252PrimeField, Stark252PrimeField>,
    pub_input: &PublicInputs,
    proof_options: &ProofOptions,
) -> bool {
    Verifier::<CairoAIR>::verify(
        proof,
        pub_input,
        proof_options,
        StoneProverTranscript::new(&[]),
    )
}

#[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
#[cfg(test)]
mod prop_test {
    use lambdaworks_math::{
        field::fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
        traits::{AsBytes, Deserializable},
    };
    use proptest::{prelude::*, prop_compose, proptest};
    use stark_platinum_prover::proof::{options::ProofOptions, stark::StarkProof};

    use crate::{
        cairo_layout::CairoLayout,
        layouts::plain::air::{generate_cairo_proof, verify_cairo_proof},
        runner::run::generate_prover_args,
        tests::utils::cairo0_program_path,
        Felt252,
    };

    use super::{MemorySegmentMap, PublicInputs, Segment, SegmentName};

    prop_compose! {
        fn some_felt()(base in any::<u64>(), exponent in any::<u128>()) -> Felt252 {
            Felt252::from(base).pow(exponent)
        }
    }

    prop_compose! {
        fn some_public_inputs()(
            pc_init in some_felt(),
            ap_init in some_felt(),
            fp_init in some_felt(),
            pc_final in some_felt(),
            ap_final in some_felt(),
            public_memory in proptest::collection::hash_map(any::<u64>(), any::<u64>(), (8_usize, 16_usize)),
            range_check_max in proptest::option::of(any::<u16>()),
            range_check_min in proptest::option::of(any::<u16>()),
            num_steps in any::<usize>(),
        ) -> PublicInputs {
            let public_memory = public_memory.iter().map(|(k, v)| (Felt252::from(*k), Felt252::from(*v))).collect();
            let memory_segments = MemorySegmentMap::from([(SegmentName::Output, Segment::new(10u64, 16u64)), (SegmentName::RangeCheck, Segment::new(20u64, 71u64))]);
            PublicInputs {
                pc_init,
                ap_init,
                fp_init,
                pc_final,
                ap_final,
                public_memory,
                range_check_max,
                range_check_min,
                num_steps,
                memory_segments,
            }
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig {cases: 5, .. ProptestConfig::default()})]
        #[test]
        fn test_public_inputs_serialization(
            public_inputs in some_public_inputs(),
        ){
            let serialized = AsBytes::as_bytes(&public_inputs);
            let deserialized: PublicInputs = Deserializable::deserialize(&serialized).unwrap();
            prop_assert_eq!(public_inputs.pc_init, deserialized.pc_init);
            prop_assert_eq!(public_inputs.ap_init, deserialized.ap_init);
            prop_assert_eq!(public_inputs.fp_init, deserialized.fp_init);
            prop_assert_eq!(public_inputs.pc_final, deserialized.pc_final);
            prop_assert_eq!(public_inputs.ap_final, deserialized.ap_final);
            prop_assert_eq!(public_inputs.public_memory, deserialized.public_memory);
            prop_assert_eq!(public_inputs.range_check_max, deserialized.range_check_max);
            prop_assert_eq!(public_inputs.range_check_min, deserialized.range_check_min);
            prop_assert_eq!(public_inputs.num_steps, deserialized.num_steps);
            prop_assert_eq!(public_inputs.memory_segments, deserialized.memory_segments);
        }
    }

    #[test]
    fn deserialize_and_verify() {
        let program_content = std::fs::read(cairo0_program_path("fibonacci_10.json")).unwrap();
        let (mut main_trace, pub_inputs) =
            generate_prover_args(&program_content, CairoLayout::Plain).unwrap();

        let proof_options = ProofOptions::default_test_options();

        // The proof is generated and serialized.
        let proof = generate_cairo_proof(&mut main_trace, &pub_inputs, &proof_options).unwrap();
        let proof_bytes: Vec<u8> = serde_cbor::to_vec(&proof).unwrap();

        // The trace and original proof are dropped to show that they are decoupled from
        // the verifying process.
        drop(main_trace);
        drop(proof);

        // At this point, the verifier only knows about the serialized proof, the proof options
        // and the public inputs.
        let proof: StarkProof<Stark252PrimeField, Stark252PrimeField> =
            serde_cbor::from_slice(&proof_bytes).unwrap();

        // The proof is verified successfully.
        assert!(verify_cairo_proof(&proof, &pub_inputs, &proof_options));
    }
}
