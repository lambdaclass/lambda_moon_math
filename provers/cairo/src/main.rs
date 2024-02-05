use lambdaworks_math::field::fields::fft_friendly::stark_252_prime_field::Stark252PrimeField;
use platinum_prover::cairo_layout::CairoLayout;
use platinum_prover::layouts::plain::air::{
    generate_cairo_proof, verify_cairo_proof, PublicInputs,
};
use platinum_prover::runner::run::generate_prover_args;
use platinum_prover::runner::run::generate_prover_args_from_trace;
use stark_platinum_prover::proof::options::{ProofOptions, SecurityLevel};
use stark_platinum_prover::proof::stark::StarkProof;
mod commands;
use clap::Parser;

use std::env;
use std::fs::File;
use std::io::{Error, ErrorKind};
use std::process::{Command, Stdio};
use std::time::Instant;

/// Get current directory and return it as a String
fn get_root_dir() -> Result<String, Error> {
    let path_buf = env::current_dir()?.canonicalize()?;
    if let Some(path) = path_buf.to_str() {
        return Ok(path.to_string());
    }

    Err(Error::new(ErrorKind::NotFound, "not found"))
}

/// Attemps to compile the Cairo program with `cairo-compile`
/// and then save it to the desired path.  
/// Returns `Ok` on success else returns `Error`
fn cairo_compile(program_path: &String, out_file_path: &String) -> Result<(), Error> {
    let out_file = File::create(out_file_path)?;

    match Command::new("cairo-compile")
        .arg("--proof_mode")
        .arg(program_path)
        .stderr(Stdio::null())
        .stdout(out_file)
        .spawn()
    {
        Ok(mut child) => {
            // wait for spawned proccess to finish
            match child.wait() {
                Ok(_) => Ok(()),
                Err(err) => Err(err),
            }
        }
        Err(err) => Err(err),
    }
}

/// Attemps to compile the Cairo program with `docker`
/// and then save it to the desired path.  
/// Returns `Ok` on success else returns `Error`
fn docker_compile(program_path: &String, out_file_path: &String) -> Result<(), Error> {
    let out_file = File::create(out_file_path)?;
    let root_dir = get_root_dir()?;
    match Command::new("docker")
        .arg("run")
        .arg("--rm")
        .arg("-v")
        .arg(format!("{}/:/pwd", root_dir))
        .arg("cairo")
        .arg("--proof_mode")
        .arg(format!("/pwd/{}", program_path))
        .stderr(Stdio::null())
        .stdout(out_file)
        .spawn()
    {
        Ok(mut child) => {
            // wait for spawned proccess to finish
            match child.wait() {
                Ok(status) => match status.code() {
                    Some(0) => Ok(()), // exit success
                    _ => Err(Error::new(
                        ErrorKind::Other,
                        "File provided is not a Cairo uncompiled",
                    )),
                },
                Err(err) => Err(err),
            }
        }
        Err(err) => Err(err),
    }
}

/// Attemps to compile the Cairo program
/// with either `cairo-compile` or `docker``
fn try_compile(program_path: &String, out_file_path: &String) -> Result<(), Error> {
    if !program_path.contains(".cairo") {
        return Err(Error::new(
            ErrorKind::Other,
            "Provided file is not a Cairo source file",
        ));
    }

    if cairo_compile(program_path, out_file_path).is_ok()
        || docker_compile(program_path, out_file_path).is_ok()
    {
        Ok(())
    } else {
        Err(Error::new(
            ErrorKind::Other,
            "Failed to compile cairo program, neither cairo-compile nor docker found",
        ))
    }
}

fn generate_proof(
    input_path: &String,
    proof_options: &ProofOptions,
) -> Option<(
    StarkProof<Stark252PrimeField, Stark252PrimeField>,
    PublicInputs,
)> {
    let timer = Instant::now();

    let Ok(program_content) = std::fs::read(input_path) else {
        eprintln!("Error opening {input_path} file");
        return None;
    };

    // FIXME: We should set this through the CLI in the future
    let layout = CairoLayout::Plain;

    let Ok((mut main_trace, pub_inputs)) = generate_prover_args(&program_content, layout) else {
        eprintln!("Error generating prover args");
        return None;
    };

    println!("  Time spent: {:?} \n", timer.elapsed());

    let timer = Instant::now();
    println!("Making proof ...");
    let proof = match generate_cairo_proof(&mut main_trace, &pub_inputs, proof_options) {
        Ok(p) => p,
        Err(err) => {
            eprintln!("Error generating proof: {:?}", err);
            return None;
        }
    };

    println!("  Time spent in proving: {:?} \n", timer.elapsed());

    Some((proof, pub_inputs))
}

fn generate_proof_from_trace(
    trace_bin_path: &str,
    memory_bin_path: &str,
    proof_options: &ProofOptions,
) -> Option<(
    StarkProof<Stark252PrimeField, Stark252PrimeField>,
    PublicInputs,
)> {
    // ## Generating the prover args
    let timer = Instant::now();
    let Ok((mut main_trace, pub_inputs)) =
        generate_prover_args_from_trace(trace_bin_path, memory_bin_path)
    else {
        eprintln!("Error generating prover args");
        return None;
    };
    println!("  Time spent: {:?} \n", timer.elapsed());

    // ## Prove
    let timer = Instant::now();
    println!("Making proof ...");
    let proof = match generate_cairo_proof(&mut main_trace, &pub_inputs, proof_options) {
        Ok(p) => p,
        Err(err) => {
            eprintln!("Error generating proof: {:?}", err);
            return None;
        }
    };
    println!("  Time spent in proving: {:?} \n", timer.elapsed());

    Some((proof, pub_inputs))
}

fn verify_proof(
    proof: StarkProof<Stark252PrimeField, Stark252PrimeField>,
    pub_inputs: PublicInputs,
    proof_options: &ProofOptions,
) -> bool {
    let timer = Instant::now();

    println!("Verifying ...");
    let proof_verified = verify_cairo_proof(&proof, &pub_inputs, proof_options);
    println!("  Time spent in verifying: {:?} \n", timer.elapsed());

    if proof_verified {
        println!("Verification succeeded");
    } else {
        println!("Verification failed");
    }

    proof_verified
}

fn write_proof(
    proof: StarkProof<Stark252PrimeField, Stark252PrimeField>,
    pub_inputs: PublicInputs,
    proof_path: String,
) {
    let mut bytes = vec![];
    let proof_bytes: Vec<u8> =
        bincode::serde::encode_to_vec(proof, bincode::config::standard()).unwrap();

    let pub_inputs_bytes: Vec<u8> =
        bincode::serde::encode_to_vec(&pub_inputs, bincode::config::standard()).unwrap();

    // This should be reworked
    // Public inputs shouldn't be stored in the proof if the verifier wants to check them

    // An u32 is enough for storing proofs up to 32 GiB
    // They shouldn't exceed the order of kbs
    // Reading an usize leads to problem in WASM (32 bit vs 64 bit architecture)

    bytes.extend((proof_bytes.len() as u32).to_le_bytes());
    bytes.extend(proof_bytes);
    bytes.extend(pub_inputs_bytes);

    let Ok(()) = std::fs::write(&proof_path, bytes) else {
        eprintln!("Error writing proof to file: {}", &proof_path);
        return;
    };

    println!("Proof written to {}", &proof_path);
}

fn main() {
    let proof_options = ProofOptions::new_secure(SecurityLevel::Conjecturable100Bits, 3);
    let args: commands::ProverArgs = commands::ProverArgs::parse();
    match args.entity {
        commands::ProverEntity::Compile(args) => {
            let out_file_path = args.program_path.replace(".cairo", ".json");
            if let Err(err) = try_compile(&args.program_path, &out_file_path) {
                eprintln!("{}", err);
            } else {
                println!("Compiled cairo program");
            }
        }
        commands::ProverEntity::RunAndProve(args) => {
            // verify input file is .cairo
            if args.program_path.contains(".cairo") {
                eprintln!("\nYou are trying to prove a non compiled Cairo program. Please compile it before sending it to the prover.\n");
                return;
            }

            let Some((proof, pub_inputs)) = generate_proof(&args.program_path, &proof_options)
            else {
                return;
            };

            write_proof(proof, pub_inputs, args.proof_path);
        }
        commands::ProverEntity::Prove(args) => {
            let Some((proof, pub_inputs)) = generate_proof_from_trace(
                &args.trace_bin_path,
                &args.memory_bin_path,
                &proof_options,
            ) else {
                return;
            };

            write_proof(proof, pub_inputs, args.proof_path);
        }
        commands::ProverEntity::Verify(args) => {
            let Ok(program_content) = std::fs::read(&args.proof_path) else {
                eprintln!("Error opening {} file", args.proof_path);
                return;
            };
            let mut bytes = program_content.as_slice();
            if bytes.len() < 8 {
                eprintln!("Error reading proof from file: {}", args.proof_path);
                return;
            }

            // Proof len was stored as an u32, 4u8 needs to be read
            let proof_len = u32::from_le_bytes(bytes[0..4].try_into().unwrap()) as usize;

            bytes = &bytes[4..];
            if bytes.len() < proof_len {
                eprintln!("Error reading proof from file: {}", args.proof_path);
                return;
            }

            let Ok((proof, _)) = bincode::serde::decode_from_slice(
                &bytes[0..proof_len],
                bincode::config::standard(),
            ) else {
                println!("Error reading proof from file: {}", args.proof_path);
                return;
            };
            bytes = &bytes[proof_len..];

            let Ok((pub_inputs, _)) =
                bincode::serde::decode_from_slice(bytes, bincode::config::standard())
            else {
                println!("Error reading proof from file: {}", args.proof_path);
                return;
            };

            verify_proof(proof, pub_inputs, &proof_options);
        }
        commands::ProverEntity::ProveAndVerify(args) => {
            if args.program_path.contains(".cairo") {
                eprintln!("\nYou are trying to prove a non compiled Cairo program. Please compile it before sending it to the prover.\n");
                return;
            }

            let Some((proof, pub_inputs)) = generate_proof(&args.program_path, &proof_options)
            else {
                return;
            };
            verify_proof(proof, pub_inputs, &proof_options);
        }
        commands::ProverEntity::CompileAndProve(args) => {
            let out_file_path = args.program_path.replace(".cairo", ".json");
            match try_compile(&args.program_path, &out_file_path) {
                Ok(_) => {
                    let Some((proof, pub_inputs)) = generate_proof(&out_file_path, &proof_options)
                    else {
                        return;
                    };

                    write_proof(proof, pub_inputs, args.proof_path);
                }
                Err(err) => {
                    eprintln!("{}", err)
                }
            }
        }
        commands::ProverEntity::CompileProveAndVerify(args) => {
            let out_file_path = args.program_path.replace(".cairo", ".json");
            match try_compile(&args.program_path, &out_file_path) {
                Ok(_) => {
                    let Some((proof, pub_inputs)) = generate_proof(&out_file_path, &proof_options)
                    else {
                        return;
                    };
                    verify_proof(proof, pub_inputs, &proof_options);
                }
                Err(err) => {
                    eprintln!("{}", err)
                }
            }
        }
    }
}
