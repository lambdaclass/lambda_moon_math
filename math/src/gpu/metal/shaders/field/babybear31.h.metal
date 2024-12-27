#pragma once

#include "fp_u32.h.metal"           // For Fp32
#include "../fft/util.h.metal"      // For reverse_index or other utility
#include "../fft/twiddles.h.metal"  // Generic twiddle templates
#include "../fft/fft.h.metal"       // The generic DIT butterfly template function
#include "../fft/permutation.h.metal"

// Define BabyBear field parameters
typedef Fp32<
    /* =N **/        88000001,      // BabyBear prime modulus
    /* =R_SQUARED **/ 1172168163,   // (R^2 mod N) precomputed
    /* =N_PRIME **/   2013265919    // Modular inverse of N under 2^32
> FpBabyBear;

// ---------------------------
// Twiddle Kernels (BabyBear)
// ---------------------------

[[kernel]] void calc_twiddles_babybear31(
    device FpBabyBear* result [[ buffer(0) ]],
    constant FpBabyBear& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
) {
    // Calls the template from `twiddles.h.metal`
    calc_twiddles<FpBabyBear>(result, _omega, index);
}

[[kernel]] void calc_twiddles_inv_babybear31(
    device FpBabyBear* result [[ buffer(0) ]],
    constant FpBabyBear& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
) {
    calc_twiddles_inv<FpBabyBear>(result, _omega, index);
}

[[kernel]] void calc_twiddles_bitrev_babybear31(
    device FpBabyBear* result [[ buffer(0) ]],
    constant FpBabyBear& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]],
    uint size [[ threads_per_grid ]]
) {
    calc_twiddles_bitrev<FpBabyBear>(result, _omega, index, size);
}

[[kernel]] void calc_twiddles_bitrev_inv_babybear31(
    device FpBabyBear* result [[ buffer(0) ]],
    constant FpBabyBear& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]],
    uint size [[ threads_per_grid ]]
) {
    calc_twiddles_bitrev_inv<FpBabyBear>(result, _omega, index, size);
}

// ------------------------------
// Montgomery DIT Butterfly (BabyBear)
// ------------------------------

[[kernel]] void radix2_dit_butterfly_babybear31(
    device FpBabyBear* input          [[ buffer(0) ]],
    constant FpBabyBear* twiddles     [[ buffer(1) ]],
    constant uint32_t& stage          [[ buffer(2) ]],
    uint32_t thread_count             [[ threads_per_grid ]],
    uint32_t thread_pos               [[ thread_position_in_grid ]]
) {
    // Calls the template from `fft.h.metal`
    // e.g., if your generic function is named `radix2_dit_butterfly_impl`
    // then do:
    radix2_dit_butterfly_impl<FpBabyBear>(input, twiddles, stage, thread_count, thread_pos);
}

// Specialized bit-reverse permutation for BabyBear
[[kernel]] void bitrev_permutation_babybear31(
    device FpBabyBear* input  [[ buffer(0) ]],
    device FpBabyBear* output [[ buffer(1) ]],
    uint index                [[ thread_position_in_grid ]],
    uint size                 [[ threads_per_grid ]]
)
{
    bitrev_permutation_impl<FpBabyBear>(input, output, index, size);
}
