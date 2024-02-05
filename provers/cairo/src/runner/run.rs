use super::vec_writer::VecWriter;
use crate::cairo_layout::CairoLayout;
use crate::cairo_mem::CairoMemory;
use crate::execution_trace::{build_cairo_execution_trace, CairoTraceTable};
use crate::layouts::plain::air::{PublicInputs, Segment, SegmentName};
use crate::register_states::RegisterStates;
use crate::Felt252;
use cairo_vm::cairo_run::{self, EncodeTraceError};
use cairo_vm::hint_processor::builtin_hint_processor::builtin_hint_processor_definition::BuiltinHintProcessor;
use cairo_vm::vm::errors::{
    cairo_run_errors::CairoRunError, trace_errors::TraceError, vm_errors::VirtualMachineError,
};
use cairo_vm::without_std::collections::HashMap;

#[derive(Debug)]
pub enum Error {
    IO(std::io::Error),
    Runner(CairoRunError),
    EncodeTrace(EncodeTraceError),
    VirtualMachine(VirtualMachineError),
    Trace(TraceError),
}

impl From<std::io::Error> for Error {
    fn from(err: std::io::Error) -> Error {
        Error::IO(err)
    }
}

impl From<CairoRunError> for Error {
    fn from(err: CairoRunError) -> Error {
        Error::Runner(err)
    }
}

impl From<EncodeTraceError> for Error {
    fn from(err: EncodeTraceError) -> Error {
        Error::EncodeTrace(err)
    }
}

impl From<VirtualMachineError> for Error {
    fn from(err: VirtualMachineError) -> Error {
        Error::VirtualMachine(err)
    }
}

impl From<TraceError> for Error {
    fn from(err: TraceError) -> Error {
        Error::Trace(err)
    }
}

/// Runs a cairo program in JSON format and returns trace, memory and program length.
/// Uses [cairo-rs](https://github.com/lambdaclass/cairo-rs/) project to run the program.
///
///  # Params
///
/// `entrypoint_function` - the name of the entrypoint function tu run. If `None` is provided, the default value is `main`.
/// `layout` - type of layout of Cairo.
/// `program_content` - content of the input file.
/// `trace_path` - path where to store the generated trace file.
/// `memory_path` - path where to store the generated memory file.
///
/// # Returns
///
/// Ok() in case of succes, with the following values:
/// - register_states
/// - cairo_mem
/// - data_len
/// - range_check: an Option<(usize, usize)> containing the start and end of range check.
/// `Error` indicating the type of error.
#[allow(clippy::type_complexity)]
pub fn run_program(
    entrypoint_function: Option<&str>,
    layout: CairoLayout,
    program_content: &[u8],
) -> Result<(RegisterStates, CairoMemory, PublicInputs), Error> {
    // default value for entrypoint is "main"
    let entrypoint = entrypoint_function.unwrap_or("main");

    let trace_enabled = true;
    let mut hint_executor = BuiltinHintProcessor::new_empty();
    let cairo_run_config = cairo_run::CairoRunConfig {
        entrypoint,
        trace_enabled,
        relocate_mem: true,
        layout: layout.as_str(),
        proof_mode: true,
        secure_run: None,
        disable_trace_padding: false,
    };

    let (runner, vm) =
        match cairo_run::cairo_run(program_content, &cairo_run_config, &mut hint_executor) {
            Ok(runner) => runner,
            Err(error) => {
                eprintln!("{error}");
                panic!();
            }
        };

    let relocated_trace = vm.get_relocated_trace().unwrap();

    let mut trace_vec = Vec::<u8>::new();
    let mut trace_writer = VecWriter::new(&mut trace_vec);
    trace_writer.write_encoded_trace(relocated_trace);

    let relocated_memory = &runner.relocated_memory;

    let mut memory_vec = Vec::<u8>::new();
    let mut memory_writer = VecWriter::new(&mut memory_vec);
    memory_writer.write_encoded_memory(relocated_memory);

    trace_writer.flush().unwrap();
    memory_writer.flush().unwrap();

    //TO DO: Better error handling
    let cairo_mem = CairoMemory::from_bytes_le(&memory_vec).unwrap();
    let register_states = RegisterStates::from_bytes_le(&trace_vec).unwrap();

    let vm_pub_inputs = runner.get_air_public_input(&vm).unwrap();

    let mut pub_memory: HashMap<Felt252, Felt252> = HashMap::new();
    vm_pub_inputs.public_memory.iter().for_each(|mem_cell| {
        let addr = Felt252::from(mem_cell.address as u64);
        let value = Felt252::from_hex_unchecked(&mem_cell.value.as_ref().unwrap().to_str_radix(16));
        pub_memory.insert(addr, value);
    });

    let mut memory_segments: HashMap<SegmentName, Segment> = HashMap::new();
    vm_pub_inputs.memory_segments.iter().for_each(|(k, v)| {
        memory_segments.insert(SegmentName::from(*k), Segment::from(v));
    });

    let num_steps = register_states.steps();
    let public_inputs = PublicInputs {
        pc_init: Felt252::from(register_states.rows[0].pc),
        ap_init: Felt252::from(register_states.rows[0].ap),
        fp_init: Felt252::from(register_states.rows[0].fp),
        pc_final: Felt252::from(register_states.rows[num_steps - 1].pc),
        ap_final: Felt252::from(register_states.rows[num_steps - 1].ap),
        range_check_min: Some(vm_pub_inputs.rc_min as u16),
        range_check_max: Some(vm_pub_inputs.rc_max as u16),
        memory_segments,
        public_memory: pub_memory,
        num_steps,
    };

    Ok((register_states, cairo_mem, public_inputs))
}

pub fn generate_prover_args(
    program_content: &[u8],
    layout: CairoLayout,
) -> Result<(CairoTraceTable, PublicInputs), Error> {
    let (register_states, memory, mut pub_inputs) = run_program(None, layout, program_content)?;

    let main_trace = build_cairo_execution_trace(&register_states, &memory, &mut pub_inputs);

    Ok((main_trace, pub_inputs))
}

pub fn generate_prover_args_from_trace(
    trace_bin_path: &str,
    memory_bin_path: &str,
) -> Result<(CairoTraceTable, PublicInputs), Error> {
    // ## Generating the prover args
    let register_states =
        RegisterStates::from_file(trace_bin_path).expect("Cairo trace bin file not found");
    let memory =
        CairoMemory::from_file(memory_bin_path).expect("Cairo memory binary file not found");

    // data length
    let data_len = 0_usize;
    let mut pub_inputs = PublicInputs::from_regs_and_mem(&register_states, &memory, data_len);

    let main_trace = build_cairo_execution_trace(&register_states, &memory, &mut pub_inputs);

    Ok((main_trace, pub_inputs))
}
