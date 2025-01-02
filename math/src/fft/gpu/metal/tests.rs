use super::*;
use lambdaworks_gpu::metal::abstractions::errors::MetalError;
use lambdaworks_gpu::metal::abstractions::state::*;
use metal::MTLSize;

#[test]
fn test_add_babybear_metal() -> Result<(), MetalError> {
    // Inicializar el estado de Metal
    let state = MetalState::new(None)?;

    // Entradas y resultados esperados
    let lhs = vec![88000000u32, 3];
    let rhs = vec![1u32, 88000000];
    let expected_output = vec![0u32, 2];
    // Suma modular para el campo BabyBear

    // Crear buffers para las entradas y la salida
    let lhs_buffer = state.alloc_buffer_data(&lhs);
    let rhs_buffer = state.alloc_buffer_data(&rhs);
    let out_buffer = state.alloc_buffer::<u32>(lhs.len());

    // Configurar el pipeline para el kernel `add_babybear`
    let pipeline = state.setup_pipeline("add_babybear")?;

    objc::rc::autoreleasepool(|| {
        let (command_buffer, command_encoder) = state.setup_command(
            &pipeline,
            Some(&[(0, &lhs_buffer), (1, &rhs_buffer), (2, &out_buffer)]),
        );

        // Configurar el tamaño de la cuadrícula y ejecutar el kernel
        let grid_size = MTLSize::new(lhs.len() as u64, 1, 1);
        let threadgroup_size = MTLSize::new(pipeline.max_total_threads_per_threadgroup(), 1, 1);

        command_encoder.dispatch_threads(grid_size, threadgroup_size);
        command_encoder.end_encoding();
        command_buffer.commit();
        command_buffer.wait_until_completed();
    });

    // Recuperar y validar los resultados
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

    // Probaremos "4 - 3 = 1 mod 88,000,001"
    let lhs = vec![4u32];
    let rhs = vec![3u32];
    let expected_output = vec![1u32]; // 4 - 3 = 1

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
    let lhs = vec![2u32];
    let rhs = vec![3u32];
    let expected_output = vec![6u32];

    // 3) Crear buffers
    let lhs_buffer = state.alloc_buffer_data(&lhs);
    let rhs_buffer = state.alloc_buffer_data(&rhs);
    let out_buffer = state.alloc_buffer::<u32>(lhs.len());

    // 4) Configurar pipeline para "mul_babybear"
    let pipeline = state.setup_pipeline("mul_babybear")?;

    objc::rc::autoreleasepool(|| {
        let (command_buffer, command_encoder) = state.setup_command(
            &pipeline,
            Some(&[(0, &lhs_buffer), (1, &rhs_buffer), (2, &out_buffer)]),
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
    assert_eq!(result, expected_output);

    Ok(())
}
