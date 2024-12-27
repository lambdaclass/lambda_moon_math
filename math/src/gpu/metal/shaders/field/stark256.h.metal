#pragma once

#include "fp_u256.h.metal"          // For Fp256
#include "../fft/util.h.metal"      // For reverse_index
#include "../fft/twiddles.h.metal"  // Generic twiddle templates
#include "../fft/fft.h.metal"       // The generic DIT butterfly template function
#include "../fft/permutation.h.metal"

// Define Stark256 field parameters
namespace {
    typedef Fp256<
        /* =N **/ /*u256(*/ 576460752303423505, 0, 0, 1 /*)*/,
        /* =R_SQUARED **/ /*u256(*/ 576413109808302096, 18446744073700081664, 5151653887, 18446741271209837569 /*)*/,
        /* =N_PRIME **/ /*u256(*/ 576460752303423504, 18446744073709551615, 18446744073709551615, 18446744073709551615 /*)*/
    > FpStark256;
}

// --------------------------
// Twiddle Kernels (Stark256)
// --------------------------

[[kernel]] void calc_twiddles_stark256(
    device FpStark256* result [[ buffer(0) ]],
    constant FpStark256& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
) {
    calc_twiddles<FpStark256>(result, _omega, index);
}

[[kernel]] void calc_twiddles_inv_stark256(
    device FpStark256* result [[ buffer(0) ]],
    constant FpStark256& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
) {
    calc_twiddles_inv<FpStark256>(result, _omega, index);
}

[[kernel]] void calc_twiddles_bitrev_stark256(
    device FpStark256* result [[ buffer(0) ]],
    constant FpStark256& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]],
    uint size [[ threads_per_grid ]]
) {
    calc_twiddles_bitrev<FpStark256>(result, _omega, index, size);
}

[[kernel]] void calc_twiddles_bitrev_inv_stark256(
    device FpStark256* result [[ buffer(0) ]],
    constant FpStark256& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]],
    uint size [[ threads_per_grid ]]
) {
    calc_twiddles_bitrev_inv<FpStark256>(result, _omega, index, size);
}

// ------------------------------
// Montgomery DIT Butterfly (Stark256)
// ------------------------------

[[kernel]] void radix2_dit_butterfly_stark256(
    device FpStark256* input          [[ buffer(0) ]],
    constant FpStark256* twiddles     [[ buffer(1) ]],
    constant uint32_t& stage          [[ buffer(2) ]],
    uint32_t thread_count             [[ threads_per_grid ]],
    uint32_t thread_pos               [[ thread_position_in_grid ]]
) {
    // e.g., if your generic function is named `radix2_dit_butterfly_impl`
    radix2_dit_butterfly_impl<FpStark256>(input, twiddles, stage, thread_count, thread_pos);
}

[[kernel]] void bitrev_permutation_stark256(
    device FpStark256* input  [[ buffer(0) ]],
    device FpStark256* output [[ buffer(1) ]],
    uint index                [[ thread_position_in_grid ]],
    uint size                 [[ threads_per_grid ]]
)
{
    bitrev_permutation_impl<FpStark256>(input, output, index, size);
}
