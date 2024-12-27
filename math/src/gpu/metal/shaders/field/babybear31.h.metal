#pragma once

#include "fp_u32.h.metal"
#include "../fft/util.h.metal"
#include "../fft/twiddles.h.metal"

// Define BabyBear field parameters
typedef Fp32<
    /* =N **/        88000001,               // BabyBear prime modulus
    /* =R_SQUARED **/ 1172168163,            // (R^2 mod N) precomputed
    /* =N_PRIME **/   2013265919             // Modular inverse of N under 2^32
> FpBabyBear;

// Compute Twiddles for BabyBear
[[kernel]] void calc_twiddles_babybear31(
    device FpBabyBear* result [[ buffer(0) ]],
    constant FpBabyBear& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
) {
    calc_twiddles<FpBabyBear>(result, _omega, index);
}

// Compute Inverse Twiddles for BabyBear
[[kernel]] void calc_twiddles_inv_babybear31(
    device FpBabyBear* result [[ buffer(0) ]],
    constant FpBabyBear& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
) {
    calc_twiddles_inv<FpBabyBear>(result, _omega, index);
}

// Compute Bit-Reversed Twiddles for BabyBear
[[kernel]] void calc_twiddles_bitrev_babybear31(
    device FpBabyBear* result [[ buffer(0) ]],
    constant FpBabyBear& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]],
    uint size [[ threads_per_grid ]]
) {
    calc_twiddles_bitrev<FpBabyBear>(result, _omega, index, size);
}

// Compute Inverse Bit-Reversed Twiddles for BabyBear
[[kernel]] void calc_twiddles_bitrev_inv_babybear31(
    device FpBabyBear* result [[ buffer(0) ]],
    constant FpBabyBear& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]],
    uint size [[ threads_per_grid ]]
) {
    calc_twiddles_bitrev_inv<FpBabyBear>(result, _omega, index, size);
}
