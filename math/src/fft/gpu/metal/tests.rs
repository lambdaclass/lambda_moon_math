use super::*;
use lambdaworks_gpu::metal::abstractions::errors::MetalError;
use lambdaworks_gpu::metal::abstractions::state::*;

#[test]
fn test_add_babybear_kernel() -> Result<(), MetalError> {
    let state = MetalState::new(None)?;

    // Example inputs and expected outputs
    let lhs = vec![1u32, 2u32, 3u32];
    let rhs = vec![2u32, 3u32, 4u32];
    let expected_output = vec![3u32, 5u32, 7u32];

    // Allocate buffers for the input and output
    let lhs_buffer = state.alloc_buffer_data(&lhs);
    let rhs_buffer = state.alloc_buffer_data(&rhs);
    let out_buffer = state.alloc_buffer::<u32>(lhs.len());

    // Setup and run the kernel
    let pipeline = state.setup_pipeline("test_add_babybear31")?;
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

    // Retrieve and verify the output
    let result = MetalState::retrieve_contents(&out_buffer);
    assert_eq!(result, expected_output);

    Ok(())
}
