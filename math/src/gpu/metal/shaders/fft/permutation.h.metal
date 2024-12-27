#pragma once
#include "util.h.metal"

template<typename Fp>
void bitrev_permutation_impl(
    device Fp* input,
    device Fp* result,
    uint index,
    uint size
)
{
    result[index] = input[reverse_index(index, size)];
}
