#ifndef felt_u32_h
#define felt_u32_h
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
    constexpr Fp32 operator+(const thread Fp32& rhs) const {
        return Fp32(add(inner, rhs.inner));
    }
    constexpr Fp32 operator-(const thread Fp32& rhs) const {
        return Fp32(sub(inner, rhs.inner));
    }
    Fp32 operator*(const thread Fp32& rhs) const {
        return Fp32(mul(inner, rhs.inner));
    }
    //============= Exponentiation and Inverse =============//
    inline Fp32 pow(uint exp) {
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
    inline Fp32 pow_naive(uint exp) {
        uint i = 1u;
        Fp32 result =  *this;
        while (i < exp) {
            result = result * *this;
            i = i + 1u;
        }
        return result;
    }

    inline Fp32 inverse() {
        // TODO: Do we have to check that the input isnt zero?
        return pow_naive(N - 2u);
    }
    inline Fp32 neg() const {
        return (inner == 0) ? Fp32(0) : Fp32(N - inner);
    }
private:
    uint inner; // The 32-bit Montgomery value
    // Compile-time constants
    constexpr static const constant uint N         = N_0;
    constexpr static const constant uint R_SQUARED = R_SQUARED_0;
    constexpr static const constant uint N_PRIME   = N_PRIME_0;
    constexpr static const constant uint R_SUB_N   = 0xFFFFFFFFu - N + 1;
    // Computes `lhs + rhs mod N`
    // Returns value in range [0,N)
    inline uint add(uint lhs, uint rhs) const {
        uint sum = lhs + rhs;
        return sum
            - (uint)(sum >= N) * N;
            // TODO: Check if doing it as u256 is better:
            // + (uint)(sum < lhs) * R_SUB_N;
    }
    // Computes `lhs - rhs mod N`
    inline uint sub(uint lhs, uint rhs) const {
        return (rhs <= lhs) ? lhs - rhs : N - (rhs - lhs);
    }
    // Montgomery multiplication
    uint mul(uint lhs, uint rhs) const {
        ulong x = (ulong)lhs * (ulong)rhs;
        ulong t = (x * (ulong)N_PRIME) & 0xFFFFFFFFull; // 32 lsb of x * mu mod 2^64
        //debug_buffer[0] = t;
        ulong u = t * (ulong)N;
        ulong x_sub_u = x - u; // sub and wrap
        uint res = (uint)(x_sub_u >> 32);
        if (x < u) {
            res += N; // add and wrap
        }
        // // NOTE: Should we remove this if?
        //  if (res >= N) {
        //      res -= N;
        //  }
        return res;
    }
};
#endif
