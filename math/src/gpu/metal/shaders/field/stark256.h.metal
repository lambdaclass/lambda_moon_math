#pragma once

#include "fp_u256.h.metal"
#include "../fft/util.h.metal"
#include "../fft/twiddles.h.metal"

// Define Stark field parameters
typedef Fp256<
    /* =N **/ /*u256(*/ 576460752303423505, 0, 0, 1 /*)*/,
    /* =R_SQUARED **/ /*u256(*/ 576413109808302096, 18446744073700081664, 5151653887, 18446741271209837569 /*)*/,
    /* =N_PRIME **/ /*u256(*/ 576460752303423504, 18446744073709551615, 18446744073709551615, 18446744073709551615 /*)*/
> FpStark256;

// Compute Twiddles for Stark256
[[kernel]] void calc_twiddles_stark256(
    device FpStark256* result [[ buffer(0) ]],
    constant FpStark256& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
) {
    calc_twiddles<FpStark256>(result, _omega, index);
}

// Compute Inverse Twiddles for Stark256
[[kernel]] void calc_twiddles_inv_stark256(
    device FpStark256* result [[ buffer(0) ]],
    constant FpStark256& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
) {
    calc_twiddles_inv<FpStark256>(result, _omega, index);
}

// Compute Bit-Reversed Twiddles for Stark256
[[kernel]] void calc_twiddles_bitrev_stark256(
    device FpStark256* result [[ buffer(0) ]],
    constant FpStark256& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]],
    uint size [[ threads_per_grid ]]
) {
    calc_twiddles_bitrev<FpStark256>(result, _omega, index, size);
}

// Compute Inverse Bit-Reversed Twiddles for Stark256
[[kernel]] void calc_twiddles_bitrev_inv_stark256(
    device FpStark256* result [[ buffer(0) ]],
    constant FpStark256& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]],
    uint size [[ threads_per_grid ]]
) {
    calc_twiddles_bitrev_inv<FpStark256>(result, _omega, index, size);
}
