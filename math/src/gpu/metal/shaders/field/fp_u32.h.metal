#ifndef felt_u32_h
#define felt_u32_h

#include <metal_stdlib>
using namespace metal;

// -------------------------------------------------------
// Plantilla para un entero de 32 bits en Montgomery mod N_0
// -------------------------------------------------------
template <
    /* =N **/ uint N_0,
    /* =R_SQUARED **/ uint R_SQUARED_0,
    /* =N_PRIME **/ uint N_PRIME_0
>
class Fp32 {
public:
    // ------------- Constructors -------------
    inline Fp32() { inner = 0u; }
    inline Fp32(uint v) { inner = v; }

    // Conversion a uint (por si necesitas)
    inline explicit operator uint() const {
        return inner;
    }

    // ------------- Operadores (+, -, *) -------------
    // Se pasan los Fp32 por **valor**, NO por referencia
    inline Fp32 operator+(Fp32 rhs) const {
        // Suma modular con corrección
        uint sum = inner + rhs.inner;
        // “carryless” mod
        if (sum >= N_0) {
            sum -= N_0;
        }
        return Fp32(sum);
    }

    inline Fp32 operator-(Fp32 rhs) const {
        // Resta modular
        uint a = inner;
        uint b = rhs.inner;
        uint diff = (a >= b) ? (a - b) : (N_0 - (b - a));
        return Fp32(diff);
    }

    inline Fp32 operator*(Fp32 rhs) const {
        // Multiplicación Montgomery
        ulong x = (ulong)inner * (ulong)rhs.inner;
        ulong t = (x * (ulong)N_PRIME) & 0xFFFFFFFFull;
        ulong u = t * (ulong)N_0;
        ulong x_sub_u = x - u;
        bool over = (x < u);

        uint res = (uint)(x_sub_u >> 32);
        if (over) {
            res += N_0;
        }
        if (res >= N_0) {
            res -= N_0;
        }
        return Fp32(res);
    }

    // ------------- Exponenciación e Inversa -------------
    inline Fp32 pow(uint exp) const {
        // “1” en Montgomery => (1 * R_SQUARED) mod N,
        // pero como simplificación (si ya usas inner=1 en Monty),
        // ajusta según tu pipeline real.
        Fp32 ONE(mulMonty(1u, R_SQUARED_0));
        
        Fp32 result = ONE;
        Fp32 base = *this; // Copia local del “this”

        while (exp > 0) {
            if (exp & 1) {
                result = result * base;
            }
            exp >>= 1;
            base = base * base;
        }
        return result;
    }

    inline Fp32 inverse() const {
        // Para campo primo => inverso = this->pow(N_0 - 2)
        return pow(N_0 - 2u);
    }

    // ------------- Neg (opcional) -------------
    inline Fp32 neg() const {
        // devuleve -x mod N
        return (inner == 0u) ? Fp32(0u) : Fp32(N_0 - inner);
    }

private:
    // Valor interno (Montgomery en 32 bits)
    uint inner;

    // ------------------------------------------------------------
    // Si Metal se queja de "static constant" aquí, podrías quitarlo y
    // hardcodear en las funciones (o usar #define).
    // ------------------------------------------------------------
    inline static constant uint N         = N_0;
    inline static constant uint R_SQUARED = R_SQUARED_0;
    inline static constant uint N_PRIME   = N_PRIME_0;
    inline static constant uint R_SUB_N   = 0xFFFFFFFFu - N_0 + 1;

    // ------------------------------------------------------------
    // Función auxiliar para "montgofy" un uint normal
    // ------------------------------------------------------------
    inline uint mulMonty(uint lhs, uint rhs) const {
        ulong x = (ulong)lhs * (ulong)rhs;
        ulong t = (x * (ulong)N_PRIME) & 0xFFFFFFFFull;
        ulong u = t * (ulong)N;
        ulong x_sub_u = x - u;
        bool over = (x < u);

        uint res = (uint)(x_sub_u >> 32);
        if (over) {
            res += N;
        }
       // if (res >= N) {
       //     res -= N;
       // }
        return res;
    }
};

#endif // felt_u32_h
