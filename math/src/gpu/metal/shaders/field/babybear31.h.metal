#pragma once

#include "fp_u32.h.metal"           // where Fp32 is defined
#include "../fft/util.h.metal"      // For reverse_index or other utility
#include "../fft/twiddles.h.metal"  // For calc_twiddles<...> etc.
#include "../fft/fft.h.metal"       // For radix2_dit_butterfly_impl<...>
#include "../fft/permutation.h.metal"

// --------------------------------------------------
// Definir BabyBear con la plantilla Fp32
// --------------------------------------------------
typedef Fp32<
    /* =N **/        2013265921,      // BabyBear prime modulus
    /* =R_SQUARED **/ 1172168163,   // (R^2 mod N) precomputed
    /* =N_PRIME **/   2281701377    // Modular inverse of N under 2^32
> FpBabyBear;

// --------------------------------------------------
// Twiddle Kernels (BabyBear)
// --------------------------------------------------

// En lugar de `constant FpBabyBear& _omega`, usamos `constant FpBabyBear* omega`.
[[kernel]] void calc_twiddles_babybear31(
    device FpBabyBear* result [[ buffer(0) ]],
    constant FpBabyBear* omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
) {
    // Llamamos la plantilla:
    // Ojo: si la firma de calc_twiddles<T> pide T por valor, pasamos `omega[0]`.
    calc_twiddles<FpBabyBear>(result, omega[0], index);
}

[[kernel]] void calc_twiddles_inv_babybear31(
    device FpBabyBear* result [[ buffer(0) ]],
    constant FpBabyBear* omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
) {
    calc_twiddles_inv<FpBabyBear>(result, omega[0], index);
}

[[kernel]] void calc_twiddles_bitrev_babybear31(
    device FpBabyBear* result [[ buffer(0) ]],
    constant FpBabyBear* omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]],
    uint size [[ threads_per_grid ]]
) {
    calc_twiddles_bitrev<FpBabyBear>(result, omega[0], index, size);
}

[[kernel]] void calc_twiddles_bitrev_inv_babybear31(
    device FpBabyBear* result [[ buffer(0) ]],
    constant FpBabyBear* omega [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]],
    uint size [[ threads_per_grid ]]
) {
    calc_twiddles_bitrev_inv<FpBabyBear>(result, omega[0], index, size);
}

// --------------------------------------------------
// Montgomery DIT Butterfly (BabyBear)
// --------------------------------------------------
[[kernel]] void radix2_dit_butterfly_babybear31(
    device FpBabyBear* input       [[ buffer(0) ]],
    constant FpBabyBear* twiddles  [[ buffer(1) ]],
    constant uint* stage           [[ buffer(2) ]],
    uint thread_count              [[ threads_per_grid ]],
    uint thread_pos                [[ thread_position_in_grid ]]
) {
    // Llamamos la genérica:
    radix2_dit_butterfly_impl<FpBabyBear>(
        input,
        twiddles,
        stage[0],
        thread_count,
        thread_pos
    );
}

// --------------------------------------------------
// bit-reverse permutation (BabyBear)
// --------------------------------------------------
[[kernel]] void bitrev_permutation_babybear31(
    device FpBabyBear* input  [[ buffer(0) ]],
    device FpBabyBear* output [[ buffer(1) ]],
    uint index                [[ thread_position_in_grid ]],
    uint size                 [[ threads_per_grid ]]
) {
    bitrev_permutation_impl<FpBabyBear>(input, output, index, size);
}

[[kernel]] void add_babybear(
    device FpBabyBear* lhs [[ buffer(0) ]],
    device FpBabyBear* rhs [[ buffer(1) ]],
    device FpBabyBear* out [[ buffer(2) ]],
    uint index [[ thread_position_in_grid ]]
) {
    // 1) Cargar en variables "thread-local"
    FpBabyBear a = lhs[index];
    FpBabyBear b = rhs[index];

    // 2) Sumar con tu operator+
    FpBabyBear c = a + b;  // <- aquí sí se usa "operator+(Fp32 rhs)"

    // 3) Escribir de vuelta
    out[index] = c;
}

[[kernel]] void sub_babybear(
    device FpBabyBear* lhs [[ buffer(0) ]],
    device FpBabyBear* rhs [[ buffer(1) ]],
    device FpBabyBear* out [[ buffer(2) ]],
    uint index [[ thread_position_in_grid ]]
)
{
    // 1) Carga a variables thread-local para poder usar operator- sin problemas
    FpBabyBear a = lhs[index];
    FpBabyBear b = rhs[index];

    // 2) Resta
    FpBabyBear c = a - b;

    // 3) Escribimos resultado de vuelta al device array
    out[index] = c;
}


[[kernel]] void mul_babybear(
    device FpBabyBear* lhs [[ buffer(0) ]],
    device FpBabyBear* rhs [[ buffer(1) ]],
    device FpBabyBear* out [[ buffer(2) ]],
    uint index [[ thread_position_in_grid ]]
)
{
    // Para que funcione el operator*, tenemos que cargar a variables "thread-local"
    FpBabyBear a = lhs[index];
    FpBabyBear b = rhs[index];

    // Multiplicamos modularmente (operator* ya definido en FpBabyBear)
    FpBabyBear c = a * b;

    // Guardamos el resultado de vuelta en `out`
    out[index] = c;
}
