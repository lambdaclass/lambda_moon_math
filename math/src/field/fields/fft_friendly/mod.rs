/// Implemenation of the Babybear Prime Field p = 2^31 - 2^27 + 1
pub mod babybear;
/// Implemenation of the quadratic extension of the babybear field
pub mod quadratic_babybear;
/// Implemenation of the quadric extension of the babybear field
pub mod quartic_babybear;
/// Implementation of the prime field used in [Stark101](https://starkware.co/stark-101/) tutorial, p = 3 * 2^30 + 1
pub mod stark_101_prime_field;
/// Implementation of two-adic prime field over 256 bit unsigned integers.
pub mod stark_252_prime_field;
/// Implemenation of the Goldilocks Prime Field p = 2^64 - 2^32 + 1
pub mod u64_goldilocks;
/// Implemenation of the Mersenne Prime field p = 2^31 - 1
pub mod u64_mersenne_montgomery_field;

/// Inmplementation of the Babybear Prime Field p = 2^31 - 2^27 + 1 using u32
pub mod babybear_u32;
pub mod quartic_babybear_u32;
