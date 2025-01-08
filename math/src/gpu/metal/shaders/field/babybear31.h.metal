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

[[kernel]] void cube_babybear(
    device FpBabyBear* input [[ buffer(0) ]],
    device FpBabyBear* output [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
)
{
    // Para que funcione el operator*, tenemos que cargar a variables "thread-local"
    FpBabyBear a = input[index];

    // Multiplicamos modularmente (operator* ya definido en FpBabyBear)
    FpBabyBear c = a * a * a;

    // Guardamos el resultado de vuelta en `out`
    output[index] = c;
}

[[kernel]] void quad_babybear(
    device FpBabyBear* input [[ buffer(0) ]],
    device FpBabyBear* output [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
)
{
    // Para que funcione el operator*, tenemos que cargar a variables "thread-local"
    FpBabyBear a = input[index];

    // Multiplicamos modularmente (operator* ya definido en FpBabyBear)
    FpBabyBear c = a * a * a * a;

    // Guardamos el resultado de vuelta en `out`
    output[index] = c;
}

[[kernel]] void pow_babybear(
    device FpBabyBear* base [[ buffer(0) ]],
    device uint* exponent [[ buffer(1) ]],
    device FpBabyBear* out [[ buffer(2) ]],
    uint index [[ thread_position_in_grid ]]
)
{
    // Para que funcione el operator*, tenemos que cargar a variables "thread-local"
    FpBabyBear a = base[index];
    uint n = exponent[index];

    // Multiplicamos modularmente (operator* ya definido en FpBabyBear)
    FpBabyBear c = a.pow(n);

    // Guardamos el resultado de vuelta en `out`
    out[index] = c;
}

[[kernel]] void inv_babybear(
    device FpBabyBear* input [[ buffer(0) ]],
    device FpBabyBear* output [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
)
{
    // Para que funcione el operator*, tenemos que cargar a variables "thread-local"
    FpBabyBear a = input[index];

    // Multiplicamos modularmente (operator* ya definido en FpBabyBear)
    FpBabyBear result = a.inverse();

    // Guardamos el resultado de vuelta en `out`
    output[index] = result;
}

[[kernel]] void mul_by_inv_babybear(
    device FpBabyBear* input [[ buffer(0) ]],
    device FpBabyBear* output [[ buffer(1) ]],
    uint index [[ thread_position_in_grid ]]
)
{
    // Para que funcione el operator*, tenemos que cargar a variables "thread-local"
    FpBabyBear a = input[index];

    // Multiplicamos modularmente (operator* ya definido en FpBabyBear)
    FpBabyBear result = a * a.inverse();

    // Guardamos el resultado de vuelta en `out`
    output[index] = result;
}

[[kernel]] void pow_and_mul_babybear(
    device FpBabyBear* base [[ buffer(0) ]],
    device uint* exp [[ buffer(1) ]],
    device FpBabyBear* output [[ buffer(2) ]],
    uint index [[ thread_position_in_grid ]]
)
{
    // Para que funcione el operator*, tenemos que cargar a variables "thread-local"
    FpBabyBear a = base[index];
    uint n = exp[index];
    FpBabyBear b = base[index];
    
    // Multiplicamos modularmente (operator* ya definido en FpBabyBear)
    FpBabyBear result = b * a.pow(n);

    // Guardamos el resultado de vuelta en `out`
    output[index] = result;
}


// [[kernel]] void mul_babybear(
//     device uint* lhs [[ buffer(0) ]],
//     device uint* rhs [[ buffer(1) ]],
//     device uint* out [[ buffer(2) ]],
//     device uint* debug_buffer [[ buffer(3) ]], // Nuevo buffer para depuración
//     uint index [[ thread_position_in_grid ]]
// )
// {
//     // Para que funcione el operator*, tenemos que cargar a variables "thread-local"
//     uint a = lhs[index];
//     uint b = rhs[index];

//     // Constantes para la operación
//     const uint N = 2013265921;        // Cambiar según el módulo de BabyBear
//     const uint N_PRIME = 2281701377; // Cambiar según el preajuste

//     ulong x = (ulong)a * (ulong)b;
//     ulong t = (x * (ulong)N_PRIME) & 0xFFFFFFFFull; // mul and wrap
//     debug_buffer[0] = t;
//     ulong u = t * (ulong)N; // t * N mod 2^{32} (mul and wrap)
//     debug_buffer[1] = u;
//     ulong x_sub_u = x - u; // sub and wrap
//     debug_buffer[2] = x_sub_u;
//     uint res = (uint)(x_sub_u >> 32);
//     if (x < u) {
//         res += N; // add and wrap
//     }

//     // if (res >= N) {
//     // res -= N;
//     // }
//     debug_buffer[3] = res;
//     // Guardamos el resultado de vuelta en `out`
//     out[index] = res;
// }
