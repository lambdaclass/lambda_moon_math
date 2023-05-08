use lambdaworks_math::field::{
    element::FieldElement,
    traits::{IsFFTField, RootsConfig},
};

use cudarc::{
    driver::{CudaDevice, LaunchAsync, LaunchConfig},
    nvrtc::safe::Ptx,
};

use crate::cuda::abstractions::{element::CUDAFieldElement, errors::CudaError, state::CudaState};

const SHADER_PTX_FFT: &str = include_str!("../shaders/fft.ptx");
const SHADER_PTX_TWIDDLES: &str = include_str!("../shaders/twiddles.ptx");

/// Executes parallel ordered FFT over a slice of two-adic field elements, in CUDA.
/// Twiddle factors are required to be in bit-reverse order.
///
/// "Ordered" means that the input is required to be in natural order, and the output will be
/// in this order too. Natural order means that input[i] corresponds to the i-th coefficient,
/// as opposed to bit-reverse order in which input[bit_rev(i)] corresponds to the i-th
/// coefficient.
pub fn fft<F>(
    input: &[FieldElement<F>],
    twiddles: &[FieldElement<F>],
    state: &CudaState,
) -> Result<Vec<FieldElement<F>>, CudaError>
where
    F: IsFFTField,
    F::BaseType: Unpin,
{
    let mut function = state.get_radix2_dit_butterfly(input, twiddles)?;

    let order = input.len().trailing_zeros();
    for stage in 0..order {
        let group_count = 1 << stage;
        let group_size = input.len() / group_count;

        function.launch(group_count, group_size)?;
    }

    let mut output = function.retrieve_result()?;

    in_place_bit_reverse_permute(&mut output);
    Ok(output)
}

pub fn gen_twiddles<F: IsFFTField>(
    order: u64,
    config: RootsConfig,
) -> Result<Vec<FieldElement<F>>, CudaError> {
    let count = (1 << order) / 2;
    let root: FieldElement<F> = F::get_primitive_root_of_unity(order)?;

    let (root, function_name) = match config {
        RootsConfig::Natural => (root, "calc_twiddles"),
        RootsConfig::NaturalInversed => (root.inv(), "calc_twiddles"),
        RootsConfig::BitReverse => (root, "calc_twiddles_bitrev"),
        RootsConfig::BitReverseInversed => (root.inv(), "calc_twiddles_bitrev"),
    };

    let device = CudaDevice::new(0).map_err(|err| CudaError::DeviceNotFound(err.to_string()))?;

    // d_ prefix is used to indicate device memory.
    let mut d_twiddles = device
        .htod_sync_copy(
            &(0..count)
                .map(|i| CUDAFieldElement::from(&FieldElement::from(i)))
                .collect::<Vec<_>>(),
        )
        .map_err(|err| CudaError::AllocateMemory(err.to_string()))?;
    let root = vec![CUDAFieldElement::from(&root)];
    let d_root = device
        .htod_sync_copy(&root)
        .map_err(|err| CudaError::AllocateMemory(err.to_string()))?;

    device
        .load_ptx(
            Ptx::from_src(SHADER_PTX_TWIDDLES),
            "twiddles",
            &["calc_twiddles_bitrev", "calc_twiddles"],
        )
        .map_err(|err| CudaError::PtxError(err.to_string()))?;
    let kernel = device
        .get_func("twiddles", function_name)
        .ok_or_else(|| CudaError::FunctionError(function_name.to_string()))?;

    let grid_dim = (1 as u32, 1, 1); // in blocks
    let block_dim = (count as u32, 1, 1);

    let config = LaunchConfig {
        grid_dim,
        block_dim,
        shared_mem_bytes: 0,
    };

    unsafe { kernel.clone().launch(config, (&mut d_twiddles, &d_root)) }
        .map_err(|err| CudaError::Launch(err.to_string()))?;

    let output = device
        .sync_reclaim(d_twiddles)
        .map_err(|err| CudaError::RetrieveMemory(err.to_string()))?;
    let output: Vec<_> = output
        .into_iter()
        .map(|cuda_elem| cuda_elem.into())
        .collect();

    Ok(output)
}

// TODO: remove after implementing in cuda
pub(crate) fn in_place_bit_reverse_permute<E>(input: &mut [E]) {
    for i in 0..input.len() {
        let bit_reversed_index = reverse_index(&i, input.len() as u64);
        if bit_reversed_index > i {
            input.swap(i, bit_reversed_index);
        }
    }
}

// TODO: remove after implementing in cuda
pub(crate) fn reverse_index(i: &usize, size: u64) -> usize {
    if size == 1 {
        *i
    } else {
        i.reverse_bits() >> (usize::BITS - size.trailing_zeros())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use lambdaworks_fft::roots_of_unity::get_twiddles;
    use lambdaworks_math::field::{
        element::FieldElement, fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
        traits::RootsConfig,
    };
    use proptest::{collection, prelude::*};

    // FFT related tests
    type F = Stark252PrimeField;
    type FE = FieldElement<F>;

    prop_compose! {
        fn powers_of_two(max_exp: u8)(exp in 1..max_exp) -> usize { 1 << exp }
        // max_exp cannot be multiple of the bits that represent a usize, generally 64 or 32.
        // also it can't exceed the test field's two-adicity.
    }

    prop_compose! {
        fn field_element()(num in any::<u64>().prop_filter("Avoid null coefficients", |x| x != &0)) -> FE {
            FE::from(num)
        }
    }

    fn field_vec(max_exp: u8) -> impl Strategy<Value = Vec<FE>> {
        powers_of_two(max_exp).prop_flat_map(|size| collection::vec(field_element(), size))
    }

    proptest! {
        #[test]
        fn test_cuda_fft_matches_sequential_fft(input in field_vec(4)) {
            let state = CudaState::new().unwrap();
            let order = input.len().trailing_zeros();
            let twiddles = get_twiddles(order.into(), RootsConfig::BitReverse).unwrap();

            let cuda_fft = fft(&input, &twiddles, &state).unwrap();
            let fft = lambdaworks_fft::ops::fft(&input, &twiddles).unwrap();

            prop_assert_eq!(cuda_fft, fft);
        }
    }
}
