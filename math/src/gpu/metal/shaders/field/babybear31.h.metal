#pragma once

#include "fp_u32.h.metal" 
#include "../fft/fft.h.metal"
#include "../fft/twiddles.h.metal"
#include "../fft/permutation.h.metal"

// Define BabyBear field parameters
typedef Fp32<
    /* =N **/        88000001,               // BabyBear prime modulus
    /* =R_SQUARED **/ 1172168163,            // (R^2 mod N) precomputed
    /* =N_PRIME **/   2013265919             // Modular inverse of N under 2^32
> FpBabyBear;

// Radix-2 DIT Butterfly kernel
[[host_name("radix2_dit_butterfly_babybear31")]] 
[[kernel]] void radix2_dit_butterfly_babybear31(
    device FpBabyBear* values,
    constant FpBabyBear* twiddles,
    constant uint32_t& stage,
    uint32_t size,
    uint32_t stride
);

// Compute Twiddles kernel
[[host_name("calc_twiddles_babybear31")]] 
[[kernel]] void calc_twiddles_babybear31(
    device FpBabyBear* twiddles,
    constant FpBabyBear& root_of_unity,
    uint size
);

// Compute Inverse Twiddles kernel
[[host_name("calc_twiddles_inv_babybear31")]] 
[[kernel]] void calc_twiddles_inv_babybear31(
    device FpBabyBear* twiddles,
    constant FpBabyBear& root_of_unity,
    uint size
);

// Compute Bit-Reversed Twiddles kernel
[[host_name("calc_twiddles_bitrev_babybear31")]]    
[[kernel]] void calc_twiddles_bitrev_babybear31(
    device FpBabyBear* twiddles,
    constant FpBabyBear& root_of_unity,
    uint size,
    uint stride
);

// Compute Inverse Bit-Reversed Twiddles kernel
[[host_name("calc_twiddles_bitrev_inv_babybear31")]] 
[[kernel]] void calc_twiddles_bitrev_inv_babybear31(
    device FpBabyBear* twiddles,
    constant FpBabyBear& root_of_unity,
    uint size,
    uint stride
);

// Perform Bit-Reversal Permutation kernel
[[host_name("bitrev_permutation_babybear31")]] 
[[kernel]] void bitrev_permutation_babybear31(
    device FpBabyBear* input,
    device FpBabyBear* output,
    uint size,
    uint stride
);
