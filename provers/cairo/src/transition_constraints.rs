use crate::Felt252;
use lambdaworks_math::field::fields::fft_friendly::stark_252_prime_field::Stark252PrimeField;
use stark_platinum_prover::{
    constraints::transition::TransitionConstraint, frame::Frame, table::TableView,
};

#[derive(Clone)]
pub struct BitPrefixFlag;
impl BitPrefixFlag {
    pub fn new() -> Self {
        Self
    }
}
impl Default for BitPrefixFlag {
    fn default() -> Self {
        Self::new()
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for BitPrefixFlag {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        0
    }

    fn evaluate(
        &self,
        frame: &stark_platinum_prover::frame::Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let constraint_idx = self.constraint_idx();

        let current_flag = current_step.get_main_evaluation_element(0, 1);
        let next_flag = current_step.get_main_evaluation_element(1, 1);

        let one = Felt252::one();
        let two = Felt252::from(2);

        let bit = current_flag - two * next_flag;

        let res = bit * (bit - one);

        transition_evaluations[constraint_idx] = res;
    }

    fn exemptions_period(&self) -> Option<usize> {
        Some(16)
    }

    fn periodic_exemptions_offset(&self) -> Option<usize> {
        Some(15)
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct ZeroFlagConstraint;
impl Default for ZeroFlagConstraint {
    fn default() -> Self {
        Self::new()
    }
}

impl ZeroFlagConstraint {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for ZeroFlagConstraint {
    fn degree(&self) -> usize {
        1
    }

    fn constraint_idx(&self) -> usize {
        1
    }

    fn evaluate(
        &self,
        frame: &stark_platinum_prover::frame::Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let zero_flag = current_step.get_main_evaluation_element(15, 1);

        transition_evaluations[self.constraint_idx()] = *zero_flag;
    }

    fn period(&self) -> usize {
        16
    }
    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct FlagOp1BaseOp0BitConstraint;
impl Default for FlagOp1BaseOp0BitConstraint {
    fn default() -> Self {
        Self::new()
    }
}

impl FlagOp1BaseOp0BitConstraint {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for FlagOp1BaseOp0BitConstraint {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        22
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let one = Felt252::one();
        let two = Felt252::from(2);

        let f_op1_imm = current_step.get_main_evaluation_element(2, 1)
            - two * current_step.get_main_evaluation_element(3, 1);
        let f_op1_fp = current_step.get_main_evaluation_element(3, 1)
            - two * current_step.get_main_evaluation_element(4, 1);
        let f_op1_ap = current_step.get_main_evaluation_element(4, 1)
            - two * current_step.get_main_evaluation_element(5, 1);

        let f_op1_base_op0_bit = one - f_op1_imm - f_op1_fp - f_op1_ap;

        let res = f_op1_base_op0_bit * (f_op1_base_op0_bit - one);

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct FlagResOp1BitConstraint;
impl Default for FlagResOp1BitConstraint {
    fn default() -> Self {
        Self::new()
    }
}

impl FlagResOp1BitConstraint {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for FlagResOp1BitConstraint {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        23
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let one = Felt252::one();
        let two = Felt252::from(2);

        let f_res_add = current_step.get_main_evaluation_element(5, 1)
            - two * current_step.get_main_evaluation_element(6, 1);
        let f_res_mul = current_step.get_main_evaluation_element(6, 1)
            - two * current_step.get_main_evaluation_element(7, 1);
        let f_pc_jnz = current_step.get_main_evaluation_element(9, 1)
            - two * current_step.get_main_evaluation_element(10, 1);

        let f_res_op1_bit = one - f_res_add - f_res_mul - f_pc_jnz;

        let res = f_res_op1_bit * (f_res_op1_bit - one);

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct FlagPcUpdateRegularBit;
impl Default for FlagPcUpdateRegularBit {
    fn default() -> Self {
        Self::new()
    }
}

impl FlagPcUpdateRegularBit {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for FlagPcUpdateRegularBit {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        24
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let one = Felt252::one();
        let two = Felt252::from(2);

        let f_jump_abs = current_step.get_main_evaluation_element(7, 1)
            - two * current_step.get_main_evaluation_element(8, 1);
        let f_jump_rel = current_step.get_main_evaluation_element(8, 1)
            - two * current_step.get_main_evaluation_element(9, 1);
        let f_pc_jnz = current_step.get_main_evaluation_element(9, 1)
            - two * current_step.get_main_evaluation_element(10, 1);

        let flag_pc_update_regular_bit = one - f_jump_abs - f_jump_rel - f_pc_jnz;

        let res = flag_pc_update_regular_bit * (flag_pc_update_regular_bit - one);

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct FlagFpUpdateRegularBit;
impl Default for FlagFpUpdateRegularBit {
    fn default() -> Self {
        Self::new()
    }
}

impl FlagFpUpdateRegularBit {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for FlagFpUpdateRegularBit {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        25
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let one = Felt252::one();
        let two = Felt252::from(2);

        let f_opcode_call = current_step.get_main_evaluation_element(12, 1)
            - two * current_step.get_main_evaluation_element(13, 1);
        let f_opcode_ret = current_step.get_main_evaluation_element(13, 1)
            - two * current_step.get_main_evaluation_element(14, 1);

        let flag_fp_update_regular_bit = one - f_opcode_call - f_opcode_ret;

        let res = flag_fp_update_regular_bit * (flag_fp_update_regular_bit - one);

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct InstructionUnpacking;
impl Default for InstructionUnpacking {
    fn default() -> Self {
        Self::new()
    }
}

impl InstructionUnpacking {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for InstructionUnpacking {
    fn degree(&self) -> usize {
        1
    }

    fn constraint_idx(&self) -> usize {
        2
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let two = Felt252::from(2);
        let b16 = two.pow(16u32);
        let b32 = two.pow(32u32);
        let b48 = two.pow(48u32);

        // Named like this to match the Cairo whitepaper's notation.
        let f0_squiggle = current_step.get_main_evaluation_element(0, 1);

        let instruction = current_step.get_main_evaluation_element(1, 3);
        let off_dst = current_step.get_main_evaluation_element(0, 0);
        let off_op0 = current_step.get_main_evaluation_element(8, 0);
        let off_op1 = current_step.get_main_evaluation_element(4, 0);

        let res = off_dst + b16 * off_op0 + b32 * off_op1 + b48 * f0_squiggle - instruction;

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct CpuOpcodesCallOff0;
impl Default for CpuOpcodesCallOff0 {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuOpcodesCallOff0 {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuOpcodesCallOff0 {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        26
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);
        let two = Felt252::from(2);
        let b15 = two.pow(15u32);

        let f_opcode_call = current_step.get_main_evaluation_element(12, 1)
            - two * current_step.get_main_evaluation_element(13, 1);

        let off_dst = current_step.get_main_evaluation_element(0, 0);

        let res = f_opcode_call * (off_dst - b15);

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct CpuOpcodesCallOff1;
impl Default for CpuOpcodesCallOff1 {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuOpcodesCallOff1 {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuOpcodesCallOff1 {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        27
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let one = Felt252::one();
        let two = Felt252::from(2);
        let b15 = two.pow(15u32);

        let f_opcode_call = current_step.get_main_evaluation_element(12, 1)
            - two * current_step.get_main_evaluation_element(13, 1);

        let off_op0 = current_step.get_main_evaluation_element(8, 0);

        let res = f_opcode_call * (off_op0 - b15 - one);

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct CpuOpcodesCallFlags;
impl Default for CpuOpcodesCallFlags {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuOpcodesCallFlags {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuOpcodesCallFlags {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        28
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let one = Felt252::one();
        let two = Felt252::from(2);

        let f_opcode_call = current_step.get_main_evaluation_element(12, 1)
            - two * current_step.get_main_evaluation_element(13, 1);

        let bit_flag0 = current_step.get_main_evaluation_element(0, 1)
            - two * current_step.get_main_evaluation_element(1, 1);
        let bit_flag1 = current_step.get_main_evaluation_element(1, 1)
            - two * current_step.get_main_evaluation_element(2, 1);

        let res =
            f_opcode_call * (two * f_opcode_call + one + one - bit_flag0 - bit_flag1 - two - two);

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct CpuOpcodesRetOff0;
impl Default for CpuOpcodesRetOff0 {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuOpcodesRetOff0 {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuOpcodesRetOff0 {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        29
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let two = Felt252::from(2);
        let b15 = two.pow(15u32);

        let f_opcode_ret = current_step.get_main_evaluation_element(13, 1)
            - two * current_step.get_main_evaluation_element(14, 1);
        let off_dst = current_step.get_main_evaluation_element(0, 0);

        let res = f_opcode_ret * (off_dst + two - b15);

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct CpuOpcodesRetOff2;
impl Default for CpuOpcodesRetOff2 {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuOpcodesRetOff2 {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuOpcodesRetOff2 {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        30
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let one = Felt252::one();
        let two = Felt252::from(2);
        let b15 = two.pow(15u32);

        let f_opcode_ret = current_step.get_main_evaluation_element(13, 1)
            - two * current_step.get_main_evaluation_element(14, 1);
        let off_op1 = current_step.get_main_evaluation_element(4, 0);

        let res = f_opcode_ret * (off_op1 + one - b15);

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct CpuOpcodesRetFlags;
impl Default for CpuOpcodesRetFlags {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuOpcodesRetFlags {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuOpcodesRetFlags {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        31
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let one = Felt252::one();
        let two = Felt252::from(2);

        let f_opcode_ret = current_step.get_main_evaluation_element(13, 1)
            - two * current_step.get_main_evaluation_element(14, 1);
        let flag0 = current_step.get_main_evaluation_element(0, 1)
            - two * current_step.get_main_evaluation_element(1, 1);
        let flag3 = current_step.get_main_evaluation_element(3, 1)
            - two * current_step.get_main_evaluation_element(4, 1);
        let flag7 = current_step.get_main_evaluation_element(7, 1)
            - two * current_step.get_main_evaluation_element(8, 1);

        let f_res_add = current_step.get_main_evaluation_element(5, 1)
            - two * current_step.get_main_evaluation_element(6, 1);
        let f_res_mul = current_step.get_main_evaluation_element(6, 1)
            - two * current_step.get_main_evaluation_element(7, 1);
        let f_pc_jnz = current_step.get_main_evaluation_element(9, 1)
            - two * current_step.get_main_evaluation_element(10, 1);

        let f_res_op1_bit = one - f_res_add - f_res_mul - f_pc_jnz;

        let res = f_opcode_ret * (flag7 + flag0 + flag3 + f_res_op1_bit - two - two);

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct CpuOperandsMemDstAddr;
impl Default for CpuOperandsMemDstAddr {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuOperandsMemDstAddr {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuOperandsMemDstAddr {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        3
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let two = Felt252::from(2);
        let one = Felt252::one();
        let b15 = two.pow(15u32);

        let dst_fp = current_step.get_main_evaluation_element(0, 1)
            - two * current_step.get_main_evaluation_element(1, 1);
        let ap = current_step.get_main_evaluation_element(0, 5);
        let fp = current_step.get_main_evaluation_element(8, 5);
        let off_dst = current_step.get_main_evaluation_element(0, 0);
        let dst_addr = current_step.get_main_evaluation_element(8, 3);

        let res = dst_fp * fp + (one - dst_fp) * ap + (off_dst - b15) - dst_addr;

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct CpuOperandsMem0Addr;
impl Default for CpuOperandsMem0Addr {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuOperandsMem0Addr {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuOperandsMem0Addr {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        4
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let two = Felt252::from(2);
        let one = Felt252::one();
        let b15 = two.pow(15u32);

        let op0_fp = current_step.get_main_evaluation_element(1, 1)
            - two * current_step.get_main_evaluation_element(2, 1);

        let ap = current_step.get_main_evaluation_element(0, 5);
        let fp = current_step.get_main_evaluation_element(8, 5);

        let off_op0 = current_step.get_main_evaluation_element(8, 0);
        let op0_addr = current_step.get_main_evaluation_element(4, 3);

        let res = op0_fp * fp + (one - op0_fp) * ap + (off_op0 - b15) - op0_addr;

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct CpuOperandsMem1Addr;
impl Default for CpuOperandsMem1Addr {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuOperandsMem1Addr {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuOperandsMem1Addr {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        5
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let one = Felt252::one();
        let two = Felt252::from(2);
        let b15 = two.pow(15u32);

        let op1_val = current_step.get_main_evaluation_element(2, 1)
            - two * current_step.get_main_evaluation_element(3, 1);
        let op1_fp = current_step.get_main_evaluation_element(3, 1)
            - two * current_step.get_main_evaluation_element(4, 1);
        let op1_ap = current_step.get_main_evaluation_element(4, 1)
            - two * current_step.get_main_evaluation_element(5, 1);

        let op0 = current_step.get_main_evaluation_element(5, 3);
        let off_op1 = current_step.get_main_evaluation_element(4, 0);
        let op1_addr = current_step.get_main_evaluation_element(12, 3);

        let ap = current_step.get_main_evaluation_element(0, 5);
        let fp = current_step.get_main_evaluation_element(8, 5);
        let pc = current_step.get_main_evaluation_element(0, 3);

        let res = op1_val * pc
            + op1_ap * ap
            + op1_fp * fp
            + (one - op1_val - op1_ap - op1_fp) * op0
            + (off_op1 - b15)
            - op1_addr;

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

// cpu/update_registers/update_ap/ap_update
pub struct CpuUpdateRegistersApUpdate;
impl Default for CpuUpdateRegistersApUpdate {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuUpdateRegistersApUpdate {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuUpdateRegistersApUpdate {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        6
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);
        let next_step = frame.get_evaluation_step(1);

        let two = Felt252::from(2);

        let ap = current_step.get_main_evaluation_element(0, 5);
        let next_ap = next_step.get_main_evaluation_element(0, 5);
        let res = current_step.get_main_evaluation_element(12, 5);

        let ap_one = current_step.get_main_evaluation_element(11, 1)
            - two * current_step.get_main_evaluation_element(12, 1);
        let opc_call = current_step.get_main_evaluation_element(12, 1)
            - two * current_step.get_main_evaluation_element(13, 1);
        let ap_add = current_step.get_main_evaluation_element(10, 1)
            - two * current_step.get_main_evaluation_element(11, 1);

        let res = ap + ap_add * res + ap_one + opc_call * two - next_ap;

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        1
    }
}

pub struct CpuUpdateRegistersFpUpdate;
impl Default for CpuUpdateRegistersFpUpdate {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuUpdateRegistersFpUpdate {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuUpdateRegistersFpUpdate {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        7
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);
        let next_step = frame.get_evaluation_step(1);

        let one = Felt252::one();
        let two = Felt252::from(2);

        let ap = current_step.get_main_evaluation_element(0, 5);
        let fp = current_step.get_main_evaluation_element(8, 5);
        let next_fp = next_step.get_main_evaluation_element(8, 5);
        let dst = current_step.get_main_evaluation_element(9, 3);

        let opc_call = current_step.get_main_evaluation_element(12, 1)
            - two * current_step.get_main_evaluation_element(13, 1);
        let opc_ret = current_step.get_main_evaluation_element(13, 1)
            - two * current_step.get_main_evaluation_element(14, 1);

        let res = opc_ret * dst + opc_call * (ap + two) + (one - opc_ret - opc_call) * fp - next_fp;

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        1
    }
}

// cpu/update_registers/update_pc/pc_cond_negative:
pub struct CpuUpdateRegistersPcCondNegative;
impl Default for CpuUpdateRegistersPcCondNegative {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuUpdateRegistersPcCondNegative {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField>
    for CpuUpdateRegistersPcCondNegative
{
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        9
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);
        let next_step = frame.get_evaluation_step(1);

        let one = Felt252::one();
        let two = Felt252::from(2);

        let t0 = current_step.get_main_evaluation_element(2, 5);
        let pc = current_step.get_main_evaluation_element(0, 3);
        let next_pc = next_step.get_main_evaluation_element(0, 3);
        let op1 = current_step.get_main_evaluation_element(13, 3);

        let pc_jnz = current_step.get_main_evaluation_element(9, 1)
            - two * current_step.get_main_evaluation_element(10, 1);
        let pc_abs = current_step.get_main_evaluation_element(7, 1)
            - two * current_step.get_main_evaluation_element(8, 1);
        let pc_rel = current_step.get_main_evaluation_element(8, 1)
            - two * current_step.get_main_evaluation_element(9, 1);
        let res = current_step.get_main_evaluation_element(12, 5);

        let res = t0 * (next_pc - (pc + op1)) + (one - pc_jnz) * next_pc
            - ((one - pc_abs - pc_rel - pc_jnz) * (pc + frame_inst_size(current_step))
                + pc_abs * res
                + pc_rel * (pc + res));
        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        1
    }
}

pub struct CpuUpdateRegistersPcCondPositive;
impl Default for CpuUpdateRegistersPcCondPositive {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuUpdateRegistersPcCondPositive {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField>
    for CpuUpdateRegistersPcCondPositive
{
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        8
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);
        let next_step = frame.get_evaluation_step(1);

        let two = Felt252::from(2);

        let t1 = current_step.get_main_evaluation_element(10, 5);
        let pc_jnz = current_step.get_main_evaluation_element(9, 1)
            - two * current_step.get_main_evaluation_element(10, 1);
        let pc = current_step.get_main_evaluation_element(0, 3);
        let next_pc = next_step.get_main_evaluation_element(0, 3);

        let res = (t1 - pc_jnz) * (next_pc - (pc + frame_inst_size(current_step)));

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        1
    }
}

//cpu/update_registers/update_pc/tmp0
pub struct CpuUpdateRegistersUpdatePcTmp0;
impl Default for CpuUpdateRegistersUpdatePcTmp0 {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuUpdateRegistersUpdatePcTmp0 {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField>
    for CpuUpdateRegistersUpdatePcTmp0
{
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        10
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let two = Felt252::from(2);
        let dst = current_step.get_main_evaluation_element(9, 3);
        let t0 = current_step.get_main_evaluation_element(2, 5);
        let pc_jnz = current_step.get_main_evaluation_element(9, 1)
            - two * current_step.get_main_evaluation_element(10, 1);

        let res = pc_jnz * dst - t0;

        transition_evaluations[self.constraint_idx()] = res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct CpuUpdateRegistersUpdatePcTmp1;
impl Default for CpuUpdateRegistersUpdatePcTmp1 {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuUpdateRegistersUpdatePcTmp1 {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField>
    for CpuUpdateRegistersUpdatePcTmp1
{
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        11
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let t1 = current_step.get_main_evaluation_element(10, 5);
        let t0 = current_step.get_main_evaluation_element(2, 5);
        let res = current_step.get_main_evaluation_element(12, 5);

        let transition_res = t0 * res - t1;

        transition_evaluations[self.constraint_idx()] = transition_res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct CpuOperandsOpsMul;
impl Default for CpuOperandsOpsMul {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuOperandsOpsMul {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuOperandsOpsMul {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        12
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let mul = current_step.get_main_evaluation_element(4, 5);
        let op0 = current_step.get_main_evaluation_element(5, 3);
        let op1 = current_step.get_main_evaluation_element(13, 3);

        transition_evaluations[self.constraint_idx()] = mul - op0 * op1;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

// cpu/operands/res
pub struct CpuOperandsRes;
impl Default for CpuOperandsRes {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuOperandsRes {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuOperandsRes {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        13
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);
        let one = Felt252::one();
        let two = Felt252::from(2);

        let mul = current_step.get_main_evaluation_element(4, 5);
        let op0 = current_step.get_main_evaluation_element(5, 3);
        let op1 = current_step.get_main_evaluation_element(13, 3);
        let res = current_step.get_main_evaluation_element(12, 5);

        let res_add = current_step.get_main_evaluation_element(5, 1)
            - two * current_step.get_main_evaluation_element(6, 1);
        let res_mul = current_step.get_main_evaluation_element(6, 1)
            - two * current_step.get_main_evaluation_element(7, 1);
        let pc_jnz = current_step.get_main_evaluation_element(9, 1)
            - two * current_step.get_main_evaluation_element(10, 1);

        let transition_res =
            res_add * (op0 + op1) + res_mul * mul + (one - res_add - res_mul - pc_jnz) * op1
                - (one - pc_jnz) * res;

        transition_evaluations[self.constraint_idx()] = transition_res;
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

// cpu/opcodes/call/push_fp
pub struct CpuOpcodesCallPushFp;
impl Default for CpuOpcodesCallPushFp {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuOpcodesCallPushFp {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuOpcodesCallPushFp {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        14
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let two = Felt252::from(2);

        let opc_call = current_step.get_main_evaluation_element(12, 1)
            - two * current_step.get_main_evaluation_element(13, 1);

        let dst = current_step.get_main_evaluation_element(9, 3);
        let fp = current_step.get_main_evaluation_element(8, 5);

        transition_evaluations[self.constraint_idx()] = opc_call * (dst - fp);
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

pub struct CpuOpcodesCallPushPc;
impl Default for CpuOpcodesCallPushPc {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuOpcodesCallPushPc {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuOpcodesCallPushPc {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        15
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let two = Felt252::from(2);

        let opc_call = current_step.get_main_evaluation_element(12, 1)
            - two * current_step.get_main_evaluation_element(13, 1);

        let op0 = current_step.get_main_evaluation_element(5, 3);
        let pc = current_step.get_main_evaluation_element(0, 3);

        transition_evaluations[self.constraint_idx()] =
            opc_call * (op0 - (pc + frame_inst_size(current_step)));
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

// cpu/opcodes/assert_eq/assert_eq
pub struct CpuOpcodesAssertEq;
impl Default for CpuOpcodesAssertEq {
    fn default() -> Self {
        Self::new()
    }
}

impl CpuOpcodesAssertEq {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for CpuOpcodesAssertEq {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        16
    }

    fn period(&self) -> usize {
        16
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let two = Felt252::from(2);

        let opc_aeq = current_step.get_main_evaluation_element(14, 1)
            - two * current_step.get_main_evaluation_element(15, 1);
        let dst = current_step.get_main_evaluation_element(9, 3);
        let res = current_step.get_main_evaluation_element(12, 5);

        transition_evaluations[self.constraint_idx()] = opc_aeq * (dst - res)
    }

    fn end_exemptions(&self) -> usize {
        0
    }
}

// memory/diff_is_bit
pub struct MemoryDiffIsBit;
impl Default for MemoryDiffIsBit {
    fn default() -> Self {
        Self::new()
    }
}

impl MemoryDiffIsBit {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for MemoryDiffIsBit {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        17
    }

    fn period(&self) -> usize {
        2
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let one = Felt252::one();

        let mem_addr_sorted = current_step.get_main_evaluation_element(0, 4);
        let mem_addr_sorted_next = current_step.get_main_evaluation_element(2, 4);

        transition_evaluations[self.constraint_idx()] = (mem_addr_sorted - mem_addr_sorted_next)
            * (mem_addr_sorted_next - mem_addr_sorted - one);
    }

    fn end_exemptions(&self) -> usize {
        1
    }
}

// memory/is_func (single-valued)
pub struct MemoryIsFunc;
impl Default for MemoryIsFunc {
    fn default() -> Self {
        Self::new()
    }
}

impl MemoryIsFunc {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for MemoryIsFunc {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        18
    }

    fn period(&self) -> usize {
        2
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let one = Felt252::one();

        let mem_addr_sorted = current_step.get_main_evaluation_element(0, 4);
        let mem_addr_sorted_next = current_step.get_main_evaluation_element(2, 4);

        let mem_val_sorted = current_step.get_main_evaluation_element(1, 4);
        let mem_val_sorted_next = current_step.get_main_evaluation_element(3, 4);

        transition_evaluations[self.constraint_idx()] =
            (mem_val_sorted - mem_val_sorted_next) * (mem_addr_sorted_next - mem_addr_sorted - one);
    }

    fn end_exemptions(&self) -> usize {
        1
    }
}

// memory/multi_column_perm/perm/step0
pub struct MemoryMultiColumnPermStep0;
impl Default for MemoryMultiColumnPermStep0 {
    fn default() -> Self {
        Self::new()
    }
}

impl MemoryMultiColumnPermStep0 {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for MemoryMultiColumnPermStep0 {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        19
    }

    fn period(&self) -> usize {
        2
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let alpha = rap_challenges[0];
        let z = rap_challenges[1];

        let a1 = current_step.get_main_evaluation_element(2, 3);
        let v1 = current_step.get_main_evaluation_element(3, 3);

        let ap1 = current_step.get_main_evaluation_element(2, 4);
        let vp1 = current_step.get_main_evaluation_element(3, 4);
        let p0 = current_step.get_aux_evaluation_element(0, 1);
        let p1 = current_step.get_aux_evaluation_element(2, 1);

        transition_evaluations[self.constraint_idx()] =
            (z - (ap1 + alpha * vp1)) * p1 - (z - (a1 + alpha * v1)) * p0;
    }

    fn end_exemptions(&self) -> usize {
        1
    }
}

// rc16/diff_is_bit
pub struct Rc16DiffIsBit;
impl Default for Rc16DiffIsBit {
    fn default() -> Self {
        Self::new()
    }
}

impl Rc16DiffIsBit {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for Rc16DiffIsBit {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        20
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        _rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);
        let one = Felt252::one();

        let rc_col_1 = current_step.get_main_evaluation_element(0, 2);
        let rc_col_2 = current_step.get_main_evaluation_element(1, 2);

        transition_evaluations[self.constraint_idx()] =
            (rc_col_1 - rc_col_2) * (rc_col_2 - rc_col_1 - one);
    }

    fn end_exemptions(&self) -> usize {
        1
    }
}

// rc16/perm/step0
pub struct Rc16PermStep0;
impl Default for Rc16PermStep0 {
    fn default() -> Self {
        Self::new()
    }
}

impl Rc16PermStep0 {
    pub fn new() -> Self {
        Self
    }
}

impl TransitionConstraint<Stark252PrimeField, Stark252PrimeField> for Rc16PermStep0 {
    fn degree(&self) -> usize {
        2
    }

    fn constraint_idx(&self) -> usize {
        21
    }

    fn evaluate(
        &self,
        frame: &Frame<Stark252PrimeField, Stark252PrimeField>,
        transition_evaluations: &mut [Felt252],
        _periodic_values: &[Felt252],
        rap_challenges: &[Felt252],
    ) {
        let current_step = frame.get_evaluation_step(0);

        let z = rap_challenges[2];
        let a1 = current_step.get_main_evaluation_element(1, 0);
        let ap1 = current_step.get_main_evaluation_element(1, 2);

        let p0 = current_step.get_aux_evaluation_element(0, 0);
        let p1 = current_step.get_aux_evaluation_element(1, 0);

        transition_evaluations[self.constraint_idx()] = (z - ap1) * p1 - (z - a1) * p0;
    }

    fn end_exemptions(&self) -> usize {
        1
    }
}

fn frame_inst_size(step: &TableView<Stark252PrimeField, Stark252PrimeField>) -> Felt252 {
    let op1_val = step.get_main_evaluation_element(2, 1)
        - Felt252::from(2) * step.get_main_evaluation_element(3, 1);
    op1_val + Felt252::one()
}
