use super::*;
use crate::field::element::FieldElement;
use crate::field::fields::fft_friendly::stark_252_prime_field::MontgomeryConfigStark252PrimeField;
use crate::field::fields::montgomery_backed_prime_fields::U256PrimeField;
use lambdaworks_gpu::metal::abstractions::errors::MetalError;
use lambdaworks_gpu::metal::abstractions::state::*;
use metal::MTLSize;

#[test]
fn test_add_babybear_metal() -> Result<(), MetalError> {
    // Inicializar estado
    let state = MetalState::new(None)?;

    // N = 2013265921
    let lhs = vec![2013265920u32, 2];
    let rhs = vec![3u32, 1];
    let expected_output = vec![2u32, 3];

    // Buffers
    let lhs_buffer = state.alloc_buffer_data(&lhs);
    let rhs_buffer = state.alloc_buffer_data(&rhs);
    let out_buffer = state.alloc_buffer::<u32>(lhs.len());

    // El pipeline es el mismo, "add_babybear"
    let pipeline = state.setup_pipeline("add_babybear")?;

    objc::rc::autoreleasepool(|| {
        let (command_buffer, command_encoder) = state.setup_command(
            &pipeline,
            Some(&[(0, &lhs_buffer), (1, &rhs_buffer), (2, &out_buffer)]),
        );

        let grid_size = MTLSize::new(lhs.len() as u64, 1, 1);
        let threadgroup_size = MTLSize::new(pipeline.max_total_threads_per_threadgroup(), 1, 1);

        command_encoder.dispatch_threads(grid_size, threadgroup_size);
        command_encoder.end_encoding();
        command_buffer.commit();
        command_buffer.wait_until_completed();
    });

    let result: Vec<u32> = MetalState::retrieve_contents(&out_buffer);
    assert_eq!(result, expected_output);

    Ok(())
}

#[test]
fn test_add_babybear_metal_two_plus_one() -> Result<(), MetalError> {
    // Inicializar estado
    let state = MetalState::new(None)?;

    // Ahora probamos "2 + 1 = 3 mod 88,000,001"
    let lhs = vec![2u32];
    let rhs = vec![1u32];
    let expected_output = vec![3u32]; // 2 + 1 = 3 < 88,000,001 => no wrap

    // Buffers
    let lhs_buffer = state.alloc_buffer_data(&lhs);
    let rhs_buffer = state.alloc_buffer_data(&rhs);
    let out_buffer = state.alloc_buffer::<u32>(lhs.len());

    // El pipeline es el mismo, "add_babybear"
    let pipeline = state.setup_pipeline("add_babybear")?;

    objc::rc::autoreleasepool(|| {
        let (command_buffer, command_encoder) = state.setup_command(
            &pipeline,
            Some(&[(0, &lhs_buffer), (1, &rhs_buffer), (2, &out_buffer)]),
        );

        let grid_size = MTLSize::new(lhs.len() as u64, 1, 1);
        let threadgroup_size = MTLSize::new(pipeline.max_total_threads_per_threadgroup(), 1, 1);

        command_encoder.dispatch_threads(grid_size, threadgroup_size);
        command_encoder.end_encoding();
        command_buffer.commit();
        command_buffer.wait_until_completed();
    });

    let result: Vec<u32> = MetalState::retrieve_contents(&out_buffer);
    assert_eq!(result, expected_output);

    Ok(())
}

#[test]
fn test_sub_babybear_metal() -> Result<(), MetalError> {
    // Inicializar el estado de Metal
    let state = MetalState::new(None)?;

    // Probaremos:
    // 4 - 3 = 1,
    // 4 - 5 = - 1 = N - 1
    let lhs = vec![4u32, 4];
    let rhs = vec![3u32, 5];
    let expected_output = vec![1u32, 2013265920]; // 4 - 3 = 1

    // Buffers
    let lhs_buffer = state.alloc_buffer_data(&lhs);
    let rhs_buffer = state.alloc_buffer_data(&rhs);
    let out_buffer = state.alloc_buffer::<u32>(lhs.len());

    // *Nuevo* pipeline "sub_babybear"
    let pipeline = state.setup_pipeline("sub_babybear")?;

    objc::rc::autoreleasepool(|| {
        let (command_buffer, command_encoder) = state.setup_command(
            &pipeline,
            Some(&[(0, &lhs_buffer), (1, &rhs_buffer), (2, &out_buffer)]),
        );

        let grid_size = MTLSize::new(lhs.len() as u64, 1, 1);
        let threadgroup_size = MTLSize::new(pipeline.max_total_threads_per_threadgroup(), 1, 1);

        command_encoder.dispatch_threads(grid_size, threadgroup_size);
        command_encoder.end_encoding();
        command_buffer.commit();
        command_buffer.wait_until_completed();
    });

    // Recuperar resultados
    let result: Vec<u32> = MetalState::retrieve_contents(&out_buffer);
    assert_eq!(result, expected_output);

    Ok(())
}

#[test]
fn test_mul_babybear_metal() -> Result<(), MetalError> {
    // 1) Inicializar el estado de Metal
    let state = MetalState::new(None)?;

    // 2) Usamos un único par de valores:
    //    lhs = [2], rhs = [3] => resultado = [6].
    let lhs = vec![4u32];
    let rhs = vec![3u32];
    let expected_output = vec![12u32];

    // 3) Crear buffers
    let lhs_buffer = state.alloc_buffer_data(&lhs);
    let rhs_buffer = state.alloc_buffer_data(&rhs);
    let out_buffer = state.alloc_buffer::<u32>(lhs.len());
    // let debug_buffer = state.alloc_buffer::<u32>(4); // Buffer para depuración

    // 4) Configurar pipeline para "mul_babybear"
    let pipeline = state.setup_pipeline("mul_babybear")?;

    objc::rc::autoreleasepool(|| {
        let (command_buffer, command_encoder) = state.setup_command(
            &pipeline,
            Some(&[
                (0, &lhs_buffer),
                (1, &rhs_buffer),
                (2, &out_buffer),
                // (3, &debug_buffer), // Asocia el buffer de depuración
            ]),
        );

        // 5) Lanzar el kernel
        let grid_size = MTLSize::new(lhs.len() as u64, 1, 1);
        let threadgroup_size = MTLSize::new(pipeline.max_total_threads_per_threadgroup(), 1, 1);

        command_encoder.dispatch_threads(grid_size, threadgroup_size);
        command_encoder.end_encoding();
        command_buffer.commit();
        command_buffer.wait_until_completed();
    });

    // 6) Recuperar resultados y comparar
    let result: Vec<u32> = MetalState::retrieve_contents(&out_buffer);
    // let debug_result: Vec<u32> = MetalState::retrieve_contents(&debug_buffer);

    let result_representative = (result[0] as u64 * (2 as u64).pow(32)) % 2013265921;

    // // Agrega validación del debug_buffer si corresponde
    // println!("DEBUG VALUE (t): {:?}", debug_result[0]);
    // println!("DEBUG VALUE (u): {:?}", debug_result[1]);
    // println!("DEBUG VALUE (x_sub_u): {:?}", debug_result[2]);
    // println!("DEBUG VALUE (res): {:?}", debug_result[3]);

    assert_eq!(result_representative, expected_output[0] as u64);

    Ok(())
}

#[test]
fn test_pow_babybear_metal() -> Result<(), MetalError> {
    let state = MetalState::new(None)?;

    let base = vec![2u32];
    let exponent = vec![3u32];
    let expected_output = vec![8u32];

    let base_buffer = state.alloc_buffer_data(&base);
    let exp_buffer = state.alloc_buffer_data(&exponent);
    let out_buffer = state.alloc_buffer::<u32>(base.len());

    let pipeline = state.setup_pipeline("pow_babybear")?;

    objc::rc::autoreleasepool(|| {
        let (command_buffer, command_encoder) = state.setup_command(
            &pipeline,
            Some(&[(0, &base_buffer), (1, &exp_buffer), (2, &out_buffer)]),
        );

        let grid_size = MTLSize::new(base.len() as u64, 1, 1);
        let threadgroup_size = MTLSize::new(pipeline.max_total_threads_per_threadgroup(), 1, 1);

        command_encoder.dispatch_threads(grid_size, threadgroup_size);
        command_encoder.end_encoding();
        command_buffer.commit();
        command_buffer.wait_until_completed();
    });

    let result: Vec<u32> = MetalState::retrieve_contents(&out_buffer);

    let result_representative = (result[0] as u64 * (2 as u64).pow(32)) % 2013265921;

    assert_eq!(result_representative, expected_output[0] as u64);

    Ok(())
}

#[test]
fn test_inv_babybear_metal() -> Result<(), MetalError> {
    let state = MetalState::new(None)?;

    let input = vec![2u32];
    let expected_output = vec![1006632961u32];

    let input_buffer = state.alloc_buffer_data(&input);
    let output_buffer = state.alloc_buffer::<u32>(input.len());

    let pipeline = state.setup_pipeline("inv_babybear")?;

    objc::rc::autoreleasepool(|| {
        let (command_buffer, command_encoder) =
            state.setup_command(&pipeline, Some(&[(0, &input_buffer), (1, &output_buffer)]));

        let grid_size = MTLSize::new(input.len() as u64, 1, 1);
        let threadgroup_size = MTLSize::new(pipeline.max_total_threads_per_threadgroup(), 1, 1);

        command_encoder.dispatch_threads(grid_size, threadgroup_size);
        command_encoder.end_encoding();
        command_buffer.commit();
        command_buffer.wait_until_completed();
    });

    let result: Vec<u32> = MetalState::retrieve_contents(&output_buffer);

    let result_representative = (result[0] as u64 * (2 as u64).pow(32)) % 2013265921;

    assert_eq!(result_representative, expected_output[0] as u64);

    Ok(())
}
