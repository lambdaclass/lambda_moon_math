use super::*;
use crate::field::element::FieldElement;
use crate::field::fields::fft_friendly::{
    babybear_u32::Babybear31PrimeField, stark_252_prime_field::MontgomeryConfigStark252PrimeField,
};
use crate::field::fields::montgomery_backed_prime_fields::U256PrimeField;
use crate::field::{
    fields::u32_montgomery_backend_prime_field::U32MontgomeryBackendPrimeField, traits::IsFFTField,
};
use lambdaworks_gpu::metal::abstractions::errors::MetalError;
use lambdaworks_gpu::metal::abstractions::state::*;
use metal::MTLSize;

/////////////// Baby Bear Field Elements ///////////////////

#[cfg(test)]
mod tests {
    use super::*;
    type FE = FieldElement<Babybear31PrimeField>;

    #[test]
    fn test_add_babybear_metal_fe() -> Result<(), MetalError> {
        // Inicializar estado
        let state = MetalState::new(None)?;

        // N = 2013265921
        let lhs = vec![FE::from(2013265920), FE::from(2)];
        let rhs = vec![FE::from(3), FE::from(1)];
        let expected_output = vec![FE::from(2), FE::from(3)];

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

        let result: Vec<FE> = MetalState::retrieve_contents(&out_buffer);
        assert_eq!(result, expected_output);

        Ok(())
    }

    #[test]
    fn test_sub_babybear_metal_fe() -> Result<(), MetalError> {
        // Inicializar el estado de Metal
        let state = MetalState::new(None)?;

        // Probaremos:
        // 4 - 3 = 1,
        // 4 - 5 = - 1 = N - 1
        let lhs = vec![FE::from(4), FE::from(4)];
        let rhs = vec![FE::from(3), FE::from(5)];
        let expected_output = vec![FE::from(1), FE::from(2013265920)];

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
        let result: Vec<FE> = MetalState::retrieve_contents(&out_buffer);
        assert_eq!(result, expected_output);

        Ok(())
    }

    #[test]
    fn test_mul_babybear_metal_fe() -> Result<(), MetalError> {
        // 1) Inicializar el estado de Metal
        let state = MetalState::new(None)?;

        // 2) Usamos varios pares de valores:
        //    lhs = [0], rhs = [3] => resultado = [0].
        let lhs = vec![
            FE::from(3),
            FE::from(2013265921),
            FE::from(2013265919),
            FE::from(1006632961),
        ];
        let rhs = vec![FE::from(4), FE::from(3), FE::from(3), FE::from(2)];
        let expected_output = vec![FE::from(12), FE::zero(), FE::from(2013265915), FE::one()];

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
        let result: Vec<FE> = MetalState::retrieve_contents(&out_buffer);
        // let debug_result: Vec<u32> = MetalState::retrieve_contents(&debug_buffer);

        assert_eq!(result, expected_output);

        Ok(())
    }

    #[test]
    fn test_cube_babybear_metal_fe() -> Result<(), MetalError> {
        let state = MetalState::new(None)?;

        let input = vec![FE::from(2)];
        let expected_output = vec![FE::from(8)];

        let input_buffer = state.alloc_buffer_data(&input);
        let out_buffer = state.alloc_buffer::<u32>(input.len());

        let pipeline = state.setup_pipeline("cube_babybear")?;

        objc::rc::autoreleasepool(|| {
            let (command_buffer, command_encoder) =
                state.setup_command(&pipeline, Some(&[(0, &input_buffer), (1, &out_buffer)]));

            let grid_size = MTLSize::new(input.len() as u64, 1, 1);
            let threadgroup_size = MTLSize::new(pipeline.max_total_threads_per_threadgroup(), 1, 1);

            command_encoder.dispatch_threads(grid_size, threadgroup_size);
            command_encoder.end_encoding();
            command_buffer.commit();
            command_buffer.wait_until_completed();
        });

        let result: Vec<FE> = MetalState::retrieve_contents(&out_buffer);

        // res * r2 mod N
        //let result_representative = (result[0] as u64 * 1172168163) % 2013265921;

        assert_eq!(result, expected_output);

        Ok(())
    }

    #[test]
    fn test_quad_babybear_metal_fe() -> Result<(), MetalError> {
        let state = MetalState::new(None)?;

        let input = vec![FE::from(2)];
        let expected_output = vec![FE::from(16)];

        let input_buffer = state.alloc_buffer_data(&input);
        let out_buffer = state.alloc_buffer::<u32>(input.len());

        let pipeline = state.setup_pipeline("quad_babybear")?;

        objc::rc::autoreleasepool(|| {
            let (command_buffer, command_encoder) =
                state.setup_command(&pipeline, Some(&[(0, &input_buffer), (1, &out_buffer)]));

            let grid_size = MTLSize::new(input.len() as u64, 1, 1);
            let threadgroup_size = MTLSize::new(pipeline.max_total_threads_per_threadgroup(), 1, 1);

            command_encoder.dispatch_threads(grid_size, threadgroup_size);
            command_encoder.end_encoding();
            command_buffer.commit();
            command_buffer.wait_until_completed();
        });

        let result: Vec<FE> = MetalState::retrieve_contents(&out_buffer);

        assert_eq!(result, expected_output);

        Ok(())
    }

    #[test]
    fn test_pow_babybear_metal_fe() -> Result<(), MetalError> {
        let state = MetalState::new(None)?;

        let base = vec![FE::from(2), FE::from(3), FE::from(2), FE::from(2)];
        let exponent = vec![3u32, 2u32, 2013265920u32, 2013265919u32];
        let expected_output = vec![FE::from(8), FE::from(9), FE::one(), FE::from(1006632961)];

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

        let result: Vec<FE> = MetalState::retrieve_contents(&out_buffer);

        assert_eq!(result, expected_output);

        Ok(())
    }

    #[test]
    fn test_mul_by_inv_babybear_metal_fe() -> Result<(), MetalError> {
        let state = MetalState::new(None)?;

        let input = vec![FE::from(2)];
        let expected_output = vec![FE::one()];

        let input_buffer = state.alloc_buffer_data(&input);
        let output_buffer = state.alloc_buffer::<u32>(input.len());

        let pipeline = state.setup_pipeline("mul_by_inv_babybear")?;

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

        let result: Vec<FE> = MetalState::retrieve_contents(&output_buffer);

        assert_eq!(result, expected_output);

        Ok(())
    }

    #[test]
    fn test_inv_babybear_metal_fe() -> Result<(), MetalError> {
        let state = MetalState::new(None)?;

        let input = vec![FE::from(3), FE::from(2)];
        let expected_output = vec![FE::from(1342177281), FE::from(1006632961)];

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

        let result: Vec<FE> = MetalState::retrieve_contents(&output_buffer);

        assert_eq!(result, expected_output);

        Ok(())
    }

    #[test]
    fn test_pow_and_mul_babybear_metal_fe() -> Result<(), MetalError> {
        let state = MetalState::new(None)?;

        let base = vec![FE::from(2)];
        let exp = vec![2013265919u32];
        let expected_output = vec![FE::one()];

        let base_buffer = state.alloc_buffer_data(&base);
        let exp_buffer = state.alloc_buffer_data(&exp);
        let output_buffer = state.alloc_buffer::<u32>(base.len());

        let pipeline = state.setup_pipeline("pow_and_mul_babybear")?;

        objc::rc::autoreleasepool(|| {
            let (command_buffer, command_encoder) = state.setup_command(
                &pipeline,
                Some(&[(0, &base_buffer), (1, &exp_buffer), (2, &output_buffer)]),
            );

            let grid_size = MTLSize::new(base.len() as u64, 1, 1);
            let threadgroup_size = MTLSize::new(pipeline.max_total_threads_per_threadgroup(), 1, 1);

            command_encoder.dispatch_threads(grid_size, threadgroup_size);
            command_encoder.end_encoding();
            command_buffer.commit();
            command_buffer.wait_until_completed();
        });

        let result: Vec<FE> = MetalState::retrieve_contents(&output_buffer);

        assert_eq!(result, expected_output);

        Ok(())
    }
    // valores reales en montgomery:
    // FE::one() = 268435454
    // FE::from(2).inv() = 134217727

    // Valor que da mal usando en la inversa pow (no naive):
    // 2 * 2^{-1} = 641659045
}
