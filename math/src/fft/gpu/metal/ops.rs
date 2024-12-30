use crate::field::{
    element::FieldElement,
    traits::{IsFFTField, IsField, IsSubFieldOf, RootsConfig},
};
use lambdaworks_gpu::metal::abstractions::{errors::MetalError, state::*};

use metal::MTLSize;

use core::mem;

/// Executes parallel ordered FFT over a slice of two-adic field elements, in Metal.
/// Twiddle factors are required to be in bit-reverse order.
///
/// "Ordered" means that the input is required to be in natural order, and the output will be
/// in this order too. Natural order means that input[i] corresponds to the i-th coefficient,
/// as opposed to bit-reverse order in which input[bit_rev(i)] corresponds to the i-th
/// coefficient.
///
/// It supports values in a field E and domain in a subfield F.
pub fn fft<F, E>(
    input: &[FieldElement<E>],
    twiddles: &[FieldElement<F>],
    state: &MetalState,
) -> Result<Vec<FieldElement<E>>, MetalError>
where
    F: IsFFTField + IsSubFieldOf<E>,
    E: IsField,
{
    // TODO: make a twiddle factor abstraction for handling invalid twiddles
    if !input.len().is_power_of_two() {
        return Err(MetalError::InputError(input.len()));
    }

    let pipeline = state.setup_pipeline(&format!("radix2_dit_butterfly_{}", F::field_name()))?;

    let input_buffer = state.alloc_buffer_data(input);
    let twiddles_buffer = state.alloc_buffer_data(twiddles);
    // TODO: twiddle factors security (right now anything can be passed as twiddle factors)

    objc::rc::autoreleasepool(|| {
        let (command_buffer, command_encoder) = state.setup_command(
            &pipeline,
            Some(&[(0, &input_buffer), (1, &twiddles_buffer)]), // index 2 is stage
        );

        let order = input.len().trailing_zeros();
        for stage in 0..order {
            command_encoder.set_bytes(2, mem::size_of_val(&stage) as u64, void_ptr(&stage));

            let grid_size = MTLSize::new(input.len() as u64 / 2, 1, 1); // one thread per butterfly
            let threadgroup_size = MTLSize::new(pipeline.thread_execution_width(), 1, 1);

            // WARN: Device should support non-uniform threadgroups (Metal3 and Apple4 or latter).
            command_encoder.dispatch_threads(grid_size, threadgroup_size);
        }
        command_encoder.end_encoding();

        command_buffer.commit();
        command_buffer.wait_until_completed();
    });

    let result = MetalState::retrieve_contents(&input_buffer);
    let result = bitrev_permutation::<F, _>(&result, state)?;
    Ok(result.into_iter().map(FieldElement::from_raw).collect())
}

/// Generates 2^{`order-1`} twiddle factors in parallel, with a certain `config`, in Metal.
pub fn gen_twiddles<F: IsFFTField>(
    order: u64,
    config: RootsConfig,
    state: &MetalState,
) -> Result<Vec<FieldElement<F>>, MetalError> {
    if order > 63 {
        return Err(MetalError::FunctionError(
            "Order should be less than or equal to 63".to_string(),
        ));
    }

    let len = (1 << order) / 2;

    let kernel = match config {
        RootsConfig::Natural => format!("calc_twiddles_{}", F::field_name()),
        RootsConfig::NaturalInversed => format!("calc_twiddles_inv_{}", F::field_name()),
        RootsConfig::BitReverse => format!("calc_twiddles_bitrev_{}", F::field_name()),
        RootsConfig::BitReverseInversed => {
            format!("calc_twiddles_bitrev_inv_{}", F::field_name())
        }
    };

    let pipeline = state.setup_pipeline(&kernel)?;

    let result_buffer = state.alloc_buffer::<F::BaseType>(len);

    objc::rc::autoreleasepool(|| {
        let (command_buffer, command_encoder) =
            state.setup_command(&pipeline, Some(&[(0, &result_buffer)]));

        let root = F::get_primitive_root_of_unity(order).unwrap();
        command_encoder.set_bytes(1, mem::size_of::<F::BaseType>() as u64, void_ptr(&root));

        let grid_size = MTLSize::new(len as u64, 1, 1);
        let threadgroup_size = MTLSize::new(pipeline.max_total_threads_per_threadgroup(), 1, 1);

        command_encoder.dispatch_threads(grid_size, threadgroup_size);
        command_encoder.end_encoding();

        command_buffer.commit();
        command_buffer.wait_until_completed();
    });

    let result = MetalState::retrieve_contents(&result_buffer);
    Ok(result.into_iter().map(FieldElement::from_raw).collect())
}

/// Executes a parallel bit-reverse permutation with the elements of `input`, in Metal.
pub fn bitrev_permutation<F: IsFFTField, T: Clone>(
    input: &[T],
    state: &MetalState,
) -> Result<Vec<T>, MetalError> {
    let pipeline = state.setup_pipeline(&format!("bitrev_permutation_{}", F::field_name()))?;

    let input_buffer = state.alloc_buffer_data(input);
    let result_buffer = state.alloc_buffer::<T>(input.len());

    objc::rc::autoreleasepool(|| {
        let (command_buffer, command_encoder) =
            state.setup_command(&pipeline, Some(&[(0, &input_buffer), (1, &result_buffer)]));

        let grid_size = MTLSize::new(input.len() as u64, 1, 1);
        let threadgroup_size = MTLSize::new(pipeline.max_total_threads_per_threadgroup(), 1, 1);

        command_encoder.dispatch_threads(grid_size, threadgroup_size);
        command_encoder.end_encoding();

        command_buffer.commit();
        command_buffer.wait_until_completed();
    });

    Ok(MetalState::retrieve_contents::<T>(&result_buffer))
}

fn void_ptr<T>(v: &T) -> *const core::ffi::c_void {
    v as *const T as *const core::ffi::c_void
}

#[cfg(test)]
mod tests {
    use crate::fft::cpu::roots_of_unity::get_twiddles;
    use crate::field::{
        fields::fft_friendly::stark_252_prime_field::Stark252PrimeField, traits::RootsConfig,
    };
    use lambdaworks_gpu::metal::abstractions::state::*;
    use proptest::{collection, prelude::*};

    use super::*;

    type F = Stark252PrimeField;
    type FE = FieldElement<F>;

    prop_compose! {
        fn powers_of_two(max_exp: u8)(exp in 1..max_exp) -> usize { 1 << exp }
        // max_exp cannot be multiple of the bits that represent a usize, generally 64 or 32.
        // also it can't exceed the test field's two-adicity.
    }

    prop_compose! {
        fn field_element()(num in any::<u64>().prop_filter("Avoid null polynomial", |x| x != &0)) -> FE {
            FE::from(num)
        }
    }

    fn field_vec(max_exp: u8) -> impl Strategy<Value = Vec<FE>> {
        powers_of_two(max_exp).prop_flat_map(|size| collection::vec(field_element(), size))
    }

    proptest! {
        // Property-based test that ensures Metal parallel FFT gives same result as a sequential one.
        #[test]
        fn test_metal_fft_matches_sequential(input in field_vec(6)) {
            let metal_state = MetalState::new(None).unwrap();
            let order = input.len().trailing_zeros();
            let twiddles = get_twiddles::<F>(order.into(), RootsConfig::BitReverse).unwrap();

            let metal_result = super::fft(&input, &twiddles, &metal_state).unwrap();
            let sequential_result = crate::fft::cpu::ops::fft(&input, &twiddles).unwrap();

            prop_assert_eq!(&metal_result, &sequential_result);
        }
    }

    // May want to modify the order constant, takes ~5s to run on a M1.
    #[test]
    fn test_metal_fft_matches_sequential_large_input() {
        const ORDER: usize = 20;
        let input = vec![FE::one(); 1 << ORDER];

        let metal_state = MetalState::new(None).unwrap();
        let order = input.len().trailing_zeros();
        let twiddles = get_twiddles::<F>(order.into(), RootsConfig::BitReverse).unwrap();

        let metal_result = super::fft(&input, &twiddles, &metal_state).unwrap();
        let sequential_result = crate::fft::cpu::ops::fft(&input, &twiddles).unwrap();

        assert_eq!(&metal_result, &sequential_result);
    }

    #[test]
    fn gen_twiddles_with_order_greater_than_63_should_fail() {
        let metal_state = MetalState::new(None).unwrap();
        let twiddles = gen_twiddles::<F>(64, RootsConfig::Natural, &metal_state);

        assert!(matches!(twiddles, Err(MetalError::FunctionError(_))));
    }

    #[test]
    fn test_babybear_addition() {
        let metal_state = MetalState::new(None).unwrap();

        // Define inputs for BabyBear field
        let lhs = vec![88000000u32, 3]; // Input elements in the field
        let rhs = vec![1u32, 88000000]; // Another set of elements in the field
        let expected_output = vec![0u32, 3]; // Expected results after modular addition

        // Allocate buffers
        let lhs_buffer = metal_state.alloc_buffer_data(&lhs);
        let rhs_buffer = metal_state.alloc_buffer_data(&rhs);
        let output_buffer = metal_state.alloc_buffer::<u32>(lhs.len());

        // Load the kernel
        let pipeline = metal_state.setup_pipeline("test_add_babybear").unwrap();

        objc::rc::autoreleasepool(|| {
            let (command_buffer, command_encoder) = metal_state.setup_command(
                &pipeline,
                Some(&[(0, &lhs_buffer), (1, &rhs_buffer), (2, &output_buffer)]),
            );

            let grid_size = MTLSize::new(lhs.len() as u64, 1, 1);
            let threadgroup_size = MTLSize::new(pipeline.max_total_threads_per_threadgroup(), 1, 1);

            command_encoder.dispatch_threads(grid_size, threadgroup_size);
            command_encoder.end_encoding();
            command_buffer.commit();
            command_buffer.wait_until_completed();
        });

        // Retrieve and validate results
        let result = MetalState::retrieve_contents::<u32>(&output_buffer);
        assert_eq!(result, expected_output);
    }
}
