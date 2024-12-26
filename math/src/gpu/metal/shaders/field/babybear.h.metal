#pragma once

#include "felt_u32.h.metal" // Usa tu implementación Fp32 de 32 bits

#include "../fft/fft.h.metal"
#include "../fft/twiddles.h.metal"
#include "../fft/permutation.h.metal"

// Prime Field of BabyBear with modulus 0x8000001D, used for Starks
namespace {
    typedef Fp32<
        /* =N **/        0x8000001D,               // BabyBear prime modulus
        /* =R_SQUARED **/ 0x23700000,              // (R^2 mod N) precomputado
        /* =N_PRIME **/   0x7FFFFFFF               // Modular inverse of N under 2^32
    > Fp;
}

// Plantilla para operar con Fp32 en FFT
template [[ host_name("radix2_dit_butterfly_babybear") ]] 
[[kernel]] void radix2_dit_butterfly<Fp>(
    device Fp*,
    constant Fp*,
    constant uint32_t&,
    uint32_t,
    uint32_t
);

template [[ host_name("calc_twiddles_babybear") ]] 
[[kernel]] void calc_twiddles<Fp>(
    device Fp*,
    constant Fp&, 
    uint
);

template [[ host_name("calc_twiddles_inv_babybear") ]] 
[[kernel]] void calc_twiddles_inv<Fp>(
    device Fp*,
    constant Fp&, 
    uint
);

template [[ host_name("calc_twiddles_bitrev_babybear") ]] 
[[kernel]] void calc_twiddles_bitrev<Fp>(
    device Fp*,
    constant Fp&, 
    uint,
    uint
);

template [[ host_name("calc_twiddles_bitrev_inv_babybear") ]] 
[[kernel]] void calc_twiddles_bitrev_inv<Fp>(
    device Fp*,
    constant Fp&, 
    uint,
    uint
);

template [[ host_name("bitrev_permutation_babybear") ]] 
[[kernel]] void bitrev_permutation<Fp>(
    device Fp*, 
    device Fp*, 
    uint, 
    uint
);
