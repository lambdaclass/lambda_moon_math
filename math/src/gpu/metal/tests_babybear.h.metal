#pragma once

#include "shaders/field/fp_u32.h.metal"
// Kernel to test addition in the BabyBear field
[[kernel]] void test_add_babybear(
    device FpBabyBear* lhs  [[ buffer(0) ]],  // Input buffer for the first operand
    device FpBabyBear* rhs  [[ buffer(1) ]],  // Input buffer for the second operand
    device FpBabyBear* out  [[ buffer(2) ]],  // Output buffer for the result
    uint index              [[ thread_position_in_grid ]] // Thread index
) {
    // Perform addition for the current thread index
    out[index] = lhs[index] + rhs[index];
}
