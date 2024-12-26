// https://github.com/andrewmilson/ministark/blob/main/gpu-poly/src/metal/felt_u256.h.metal

#ifndef felt_u32_h
#define felt_u32_h

// #include "../unsigned_integer/u256.h.metal" // se deberia ir?

template <
    /* =N **/ uint N_0,
    /* =R_SQUARED **/ uint R_SQUARED_0,
    /* =N_PRIME **/ uint N_PRIME_0
>
class Fp32 {
public:
Fp32() = default;

    constexpr Fp32(uint v) : inner(v) {}

    constexpr explicit operator uint() const {
        return inner;
    }

    //============= Operator Overloads (mod N) =============//

     // 1) Addition
    constexpr Fp32 operator+(const Fp32& rhs) const {
        return Fp32(add(inner, rhs.inner));
    }

    // 2) Subtraction
    constexpr Fp32 operator-(const Fp32& rhs) const {
        return Fp32(sub(inner, rhs.inner));
    }

    // 3) Multiplication (Montgomery)
    Fp32 operator*(const Fp32& rhs) const {
        return Fp32(mul(inner, rhs.inner));
    }

    //============= Exponentiation and Inverse =============//


    Fp32 pow(uint exp) {
        // "Montgomery 1" = (1 * R^2) mod N
        Fp32 const ONE = Fp32(mul(1u, R_SQUARED));
        Fp32 result = ONE;

        while (exp > 0) {
            if (exp & 1) {
                result = result * *this;
            }
            exp >>= 1;
            *this = *this * *this;
        }
        return result;
    }
     // Inverse using Fermat's little theorem: a^(N-2) mod N
    Fp32 inverse() {
        return pow(N - 2u);
    }

// If inner == 0, result is 0; otherwise it's N - inner.
 Fp32 neg() const {
        return (inner == 0) ? Fp32(0) : Fp32(N - inner);
    }

private:
    uint inner; // The 32-bit Montgomery value

    // Compile-time constants
    constexpr static const uint N         = N_0;
    constexpr static const uint R_SQUARED = R_SQUARED_0;
    constexpr static const uint N_PRIME   = N_PRIME_0;
// (1 << 32) - N + 1, used in the add function for overflow correction
    constexpr static const uint R_SUB_N   = 0xFFFFFFFFu - N + 1;


    // Computes `lhs + rhs mod N`
    // Returns value in range [0,N)
inline uint add(uint lhs, uint rhs) const
{
    uint sum = lhs + rhs;
    // sum < lhs => 32-bit overflow occurred
    // sum >= N => sum >= modulus
    return sum
        - (uint)(sum >= N) * N
        + (uint)(sum < lhs) * R_SUB_N;
}

    // Computes `lhs - rhs mod N`
    // Assumes `rhs` value in range [0,N)
    inline uint sub(uint lhs, uint rhs) const {
        if (rhs <= lhs) {
            return lhs - rhs;
        } else {
            return N - (rhs - lhs);
        }
    }

// 3) Mul mod N (Montgomery multiplication)s
uint mul(uint lhs, uint rhs) const {
        // 64-bit product
        ulong x = (ulong)lhs * (ulong)rhs;

        // Partial reduction: (x * N_PRIME) mod 2^32
        ulong t = (x * (ulong)N_PRIME) & 0xFFFFFFFFull;

        // Multiply by modulus
        ulong u = t * (ulong)N;

        // Subtract in 64 bits and handle underflow
        bool underflow = (x < u);
        ulong x_sub_u = x - u;

        // Take the high 32 bits
        uint res = (uint)(x_sub_u >> 32);

        // Handle underflow and final reduction
        if (underflow) {
            res += N;
        }
        if (res >= N) {
            res -= N;
        }
        return res;
    }
};

#endif 
