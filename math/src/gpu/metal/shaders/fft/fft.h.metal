#pragma once
#include <metal_stdlib>
using namespace metal;

// Generic radix-2 DIT butterfly template (no kernel for specific types)
template<typename Fp>
void radix2_dit_butterfly_impl(
    device Fp* input,
    constant Fp* twiddles,
    constant uint32_t& stage,
    uint32_t thread_count,
    uint32_t thread_pos
)
{
    uint32_t half_group_size = thread_count >> stage;
    uint32_t group = thread_pos >> metal::ctz(half_group_size);

    uint32_t pos_in_group = thread_pos & (half_group_size - 1);
    uint32_t i = thread_pos * 2 - pos_in_group;

    Fp w = twiddles[group];
    Fp a = input[i];
    Fp b = input[i + half_group_size];

    Fp res_1 = a + w * b;
    Fp res_2 = a - w * b;

    input[i]                    = res_1;
    input[i + half_group_size]  = res_2;
}
