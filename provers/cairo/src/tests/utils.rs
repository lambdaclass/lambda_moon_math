use crate::{
    cairo_layout::CairoLayout,
    layouts::plain::air::{generate_cairo_proof, verify_cairo_proof},
    runner::run::generate_prover_args,
    runner::run::generate_prover_args_from_trace,
};
use stark_platinum_prover::proof::options::ProofOptions;
use std::time::Instant;

pub fn cairo0_program_path(program_name: &str) -> String {
    const CARGO_DIR: &str = env!("CARGO_MANIFEST_DIR");
    const CAIRO0_BASE_REL_PATH: &str = "/cairo_programs/cairo0/";
    let program_base_path = CARGO_DIR.to_string() + CAIRO0_BASE_REL_PATH;
    program_base_path + program_name
}

pub fn cairo1_program_path(program_name: &str) -> String {
    const CARGO_DIR: &str = env!("CARGO_MANIFEST_DIR");
    const CAIRO1_BASE_REL_PATH: &str = "/cairo_programs/cairo1/";
    let program_base_path = CARGO_DIR.to_string() + CAIRO1_BASE_REL_PATH;
    program_base_path + program_name
}

/// Loads the program in path, runs it with the Cairo VM, and makes a proof of it
pub fn test_prove_cairo_program(file_path: &str, layout: CairoLayout) {
    let proof_options = ProofOptions::default_test_options();
    let timer = Instant::now();
    println!("Making proof ...");

    let program_content = std::fs::read(file_path).unwrap();
    let (mut main_trace, pub_inputs) = generate_prover_args(&program_content, layout).unwrap();
    let proof = generate_cairo_proof(&mut main_trace, &pub_inputs, &proof_options).unwrap();
    println!("  Time spent in proving: {:?} \n", timer.elapsed());

    assert!(verify_cairo_proof(&proof, &pub_inputs, &proof_options));
}

pub fn test_prove_cairo_program_from_trace(trace_bin_path: &str, memory_bin_path: &str) {
    let proof_options = ProofOptions::default_test_options();
    let (mut main_trace, pub_inputs) =
        generate_prover_args_from_trace(trace_bin_path, memory_bin_path).unwrap();

    // println
    let timer = Instant::now();
    println!("Making proof ...");
    let proof = generate_cairo_proof(&mut main_trace, &pub_inputs, &proof_options).unwrap();
    println!("  Time spent in proving: {:?} \n", timer.elapsed());
    assert!(verify_cairo_proof(&proof, &pub_inputs, &proof_options));
}
