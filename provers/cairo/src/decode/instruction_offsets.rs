use crate::Felt252;
use lambdaworks_math::field::{element::FieldElement, traits::IsField};

use super::instruction_flags::aux_get_last_nim_of_field_element;

const OFF_DST_OFF: u32 = 0;
const OFF_OP0_OFF: u32 = 16;
const OFF_OP1_OFF: u32 = 32;
const OFFX_MASK: u64 = 0xFFFF;

#[derive(Debug, PartialEq, Eq)]
pub struct InstructionOffsets {
    pub off_dst: i32,
    pub off_op0: i32,
    pub off_op1: i32,
}

impl InstructionOffsets {
    pub fn new(mem_value: &Felt252) -> Self {
        Self {
            off_dst: Self::decode_offset(mem_value, OFF_DST_OFF),
            off_op0: Self::decode_offset(mem_value, OFF_OP0_OFF),
            off_op1: Self::decode_offset(mem_value, OFF_OP1_OFF),
        }
    }

    pub fn decode_offset(mem_value: &Felt252, instruction_offset: u32) -> i32 {
        let offset = aux_get_last_nim_of_field_element(mem_value) >> instruction_offset & OFFX_MASK;
        let vectorized_offset = offset.to_le_bytes();
        let aux = [
            vectorized_offset[0],
            vectorized_offset[1].overflowing_sub(128).0,
        ];
        i32::from(i16::from_le_bytes(aux))
    }

    pub fn to_trace_representation(&self) -> (Felt252, Felt252, Felt252) {
        (
            to_unbiased_representation(self.off_dst),
            to_unbiased_representation(self.off_op0),
            to_unbiased_representation(self.off_op1),
        )
    }
}

/// Returns an unbiased representation of the number. This is applied to
/// instruction offsets as explained in section 9.4 of the Cairo whitepaper
/// to be in the range [0, 2^16). https://eprint.iacr.org/2021/1063.pdf
fn to_unbiased_representation<F: IsField>(n: i32) -> FieldElement<F> {
    let b15 = 2u64.pow(15u32);
    if n < 0 {
        FieldElement::<F>::from(b15 - n.unsigned_abs() as u64)
    } else {
        FieldElement::<F>::from(n as u64 + b15)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assert_opcode_flag_is_correct_1() {
        // Instruction A
        let value = Felt252::from(0x480680017fff8000);
        let instruction_offsets = InstructionOffsets::new(&value);

        assert_eq!(instruction_offsets.off_dst, 0);
        assert_eq!(instruction_offsets.off_op0, -1);
        assert_eq!(instruction_offsets.off_op1, 1);
    }

    #[test]
    fn assert_opcode_flag_is_correct_2() {
        // Instruction A
        let value = Felt252::from(0x208b7fff7fff7ffe);
        let instruction_offsets = InstructionOffsets::new(&value);

        assert_eq!(instruction_offsets.off_dst, -2);
        assert_eq!(instruction_offsets.off_op0, -1);
        assert_eq!(instruction_offsets.off_op1, -1);
    }

    #[test]
    fn assert_opcode_flag_is_correct_3() {
        // Instruction A
        let value = Felt252::from(0x48327ffc7ffa8000);
        let instruction_offsets = InstructionOffsets::new(&value);

        assert_eq!(instruction_offsets.off_dst, 0);
        assert_eq!(instruction_offsets.off_op0, -6);
        assert_eq!(instruction_offsets.off_op1, -4);
    }
}
