// This is necessary because pragma once doesn't work as expected
// and some symbols are being duplicated.

// TODO: Investigate this issue, having .metal sources would be better
// than headers and a unique source.

#include "shaders/field/babybear31.h.metal"
#include "shaders/field/stark256.h.metal"

