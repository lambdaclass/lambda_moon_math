#pragma once
#include "util.h.metal"

// Template for Twiddles
template<typename Fp>
void calc_twiddles(
    device Fp* result [[ buffer(0) ]],
    constant Fp& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
) {
    Fp omega = _omega;
    result[index] = omega.pow(index);
}

template<typename Fp>
void calc_twiddles_inv(
    device Fp* result [[ buffer(0) ]],
    constant Fp& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
) {
    Fp omega = _omega;
    result[index] = omega.pow(index).inverse();
}

template<typename Fp>
void calc_twiddles_bitrev(
    device Fp* result [[ buffer(0) ]],
    constant Fp& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]],
    uint size [[ threads_per_grid ]]
) {
    Fp omega = _omega;
    result[index] = omega.pow(reverse_index(index, size));
}

template<typename Fp>
void calc_twiddles_bitrev_inv(
    device Fp* result [[ buffer(0) ]],
    constant Fp& _omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]],
    uint size [[ threads_per_grid ]]
) {
    Fp omega = _omega;
    result[index] = omega.pow(reverse_index(index, size)).inverse();
}
