use std::{collections::VecDeque, iter};

use super::{
    cairo_mem::CairoMemory,
    decode::{
        instruction_flags::{
            aux_get_last_nim_of_field_element, ApUpdate, CairoInstructionFlags, CairoOpcode,
            DstReg, Op0Reg, Op1Src, PcUpdate, ResLogic,
        },
        instruction_offsets::InstructionOffsets,
    },
    register_states::RegisterStates,
};
use crate::layouts::plain::air::PublicInputs;
use cairo_vm::without_std::collections::HashMap;
use itertools::Itertools;
use lambdaworks_math::{
    field::fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
    unsigned_integer::element::UnsignedInteger,
};
use stark_platinum_prover::{fri::FieldElement, trace::TraceTable, Felt252};
pub type CairoTraceTable = TraceTable<Stark252PrimeField, Stark252PrimeField>;

// NOTE: This should be deleted and use CairoAIR::STEP_SIZE once it is set to 16
const CAIRO_STEP: usize = 16;

const PLAIN_LAYOUT_NUM_MAIN_COLUMNS: usize = 6;
const PLAIN_LAYOUT_NUM_AUX_COLUMNS: usize = 2;

/// Gets holes from the range-checked columns. These holes must be filled for the
/// permutation range-checks, as can be read in section 9.9 of the Cairo whitepaper.
/// Receives the trace and the indexes of the range-checked columns.
/// Outputs the holes that must be filled to make the range continuous and the extreme
/// values rc_min and rc_max, corresponding to the minimum and maximum values of the range.
/// NOTE: These extreme values should be received as public inputs in the future and not
/// calculated here.
fn get_rc_holes(sorted_rc_values: &[u16]) -> VecDeque<Felt252> {
    let mut rc_holes = VecDeque::new();

    let mut prev_rc_value = sorted_rc_values[0];

    for rc_value in sorted_rc_values.iter() {
        let rc_diff = rc_value - prev_rc_value;
        if rc_diff != 1 && rc_diff != 0 {
            let mut rc_hole = prev_rc_value + 1;

            while rc_hole < *rc_value {
                rc_holes.push_back(Felt252::from(rc_hole as u64));
                rc_hole += 1;
            }
        }
        prev_rc_value = *rc_value;
    }

    rc_holes
}

/// Get memory holes from accessed addresses. These memory holes appear
/// as a consequence of interaction with builtins.
/// Returns a vector of addresses that were not present in the input vector (holes)
///
/// # Arguments
///
/// * `sorted_addrs` - Vector of sorted memory addresses.
/// * `codelen` - the length of the Cairo program instructions.
fn get_memory_holes(
    sorted_addrs: &[Felt252],
    pub_memory: &HashMap<Felt252, Felt252>,
) -> VecDeque<Felt252> {
    let mut memory_holes = VecDeque::new();
    let mut prev_addr = &sorted_addrs[0];
    let one = Felt252::one();

    for addr in sorted_addrs.iter() {
        let addr_diff = addr - prev_addr;

        // If the candidate memory hole has an address belonging to the program segment (public
        // memory), that is not accounted here since public memory is added in a posterior step of
        // the protocol.
        if addr_diff != one && addr_diff != Felt252::zero() {
            let mut hole_addr = prev_addr + one;

            while hole_addr.representative() < addr.representative() {
                if !pub_memory.contains_key(&hole_addr) {
                    memory_holes.push_back(hole_addr);
                }
                hole_addr += one;
            }
        }
        prev_addr = addr;
    }

    let max_addr_plus_one = sorted_addrs.last().unwrap() + one;
    memory_holes.push_back(max_addr_plus_one);

    memory_holes
}

/// Receives the raw Cairo trace and memory as outputted from the Cairo VM and returns
/// the trace table used to Felt252ed the Cairo STARK prover.
/// The constraints of the Cairo AIR are defined over this trace rather than the raw trace
/// obtained from the Cairo VM, this is why this function is needed.
pub fn build_cairo_execution_trace(
    register_states: &RegisterStates,
    memory: &CairoMemory,
    public_input: &mut PublicInputs,
) -> CairoTraceTable {
    let num_steps = register_states.steps();

    // Instruction flags and offsets are decoded from the raw instructions and represented
    // by the CairoInstructionFlags and InstructionOffsets as an intermediate representation
    let (flags, biased_offsets): (Vec<CairoInstructionFlags>, Vec<InstructionOffsets>) =
        register_states
            .flags_and_offsets(memory)
            .unwrap()
            .into_iter()
            .unzip();

    // dst, op0, op1 and res are computed from flags and offsets
    let (dst_addrs, mut dsts): (Vec<Felt252>, Vec<Felt252>) =
        compute_dst(&flags, &biased_offsets, register_states, memory);
    let (op0_addrs, mut op0s): (Vec<Felt252>, Vec<Felt252>) =
        compute_op0(&flags, &biased_offsets, register_states, memory);
    let (op1_addrs, op1s): (Vec<Felt252>, Vec<Felt252>) =
        compute_op1(&flags, &biased_offsets, register_states, memory, &op0s);
    let mut res = compute_res(&flags, &op0s, &op1s, &dsts);

    // In some cases op0, dst or res may need to be updated from the already calculated values
    update_values(&flags, register_states, &mut op0s, &mut dsts, &mut res);

    // Flags and offsets are transformed to a bit representation. This is needed since
    // the flag constraints of the Cairo AIR are defined over bit representations of these
    let bit_prefix_flags: Vec<[Felt252; 16]> = flags
        .iter()
        .map(CairoInstructionFlags::to_trace_representation)
        .collect();
    let unbiased_offsets: Vec<(Felt252, Felt252, Felt252)> = biased_offsets
        .iter()
        .map(InstructionOffsets::to_trace_representation)
        .collect();

    // ap, fp, pc and instruction columns are computed
    let aps: Vec<Felt252> = register_states
        .rows
        .iter()
        .map(|t| Felt252::from(t.ap))
        .collect();

    let fps: Vec<Felt252> = register_states
        .rows
        .iter()
        .map(|t| Felt252::from(t.fp))
        .collect();

    let pcs: Vec<Felt252> = register_states
        .rows
        .iter()
        .map(|t| Felt252::from(t.pc))
        .collect();

    let instructions: Vec<Felt252> = register_states
        .rows
        .iter()
        .map(|t| *memory.get(&t.pc).unwrap())
        .collect();

    // t0, t1 and mul derived values are constructed. For details reFelt252r to
    // section 9.1 of the Cairo whitepaper
    let two = Felt252::from(2);
    let t0: Vec<Felt252> = bit_prefix_flags
        .iter()
        .zip(&dsts)
        .map(|(repr_flags, dst)| (repr_flags[9] - two * repr_flags[10]) * dst)
        .collect();
    let t1: Vec<Felt252> = t0.iter().zip(&res).map(|(t, r)| t * r).collect();
    let mul: Vec<Felt252> = op0s.iter().zip(&op1s).map(|(op0, op1)| op0 * op1).collect();

    let mut trace: CairoTraceTable = TraceTable::allocate_with_zeros(
        num_steps,
        PLAIN_LAYOUT_NUM_MAIN_COLUMNS,
        PLAIN_LAYOUT_NUM_AUX_COLUMNS,
        CAIRO_STEP,
    );

    let rc_values = set_rc_pool(&mut trace, unbiased_offsets);
    set_bit_prefix_flags(&mut trace, bit_prefix_flags);
    let mut sorted_addrs = set_mem_pool(
        &mut trace,
        pcs,
        instructions,
        op0_addrs,
        op0s,
        dst_addrs,
        dsts,
        op1_addrs,
        op1s,
    );
    set_update_pc(&mut trace, aps, t0, t1, mul, fps, res);

    // Sort values in rc pool
    let mut sorted_rc_value_representatives: Vec<u16> = rc_values
        .iter()
        .map(|x| x.representative().into())
        .collect();
    sorted_rc_value_representatives.sort();
    let rc_holes = get_rc_holes(&sorted_rc_value_representatives);
    let rc_max = Felt252::from(*(sorted_rc_value_representatives.last().unwrap()) as u64);
    finalize_rc_pool(&mut trace, rc_holes, rc_max);

    // Get all rc values.
    // NOTE: We are sorting these values again, once for finding rc holes and one for the sorted column construction.
    // This could be rethinked for better performance
    let mut sorted_rc_column = trace.get_column_main(0);
    sorted_rc_column.sort_by_key(|x| x.representative());
    set_sorted_rc_pool(&mut trace, sorted_rc_column);

    // Add memory holes
    sorted_addrs.sort_by_key(|x| x.representative());
    let memory_holes = get_memory_holes(&sorted_addrs, &public_input.public_memory);
    finalize_mem_pool(&mut trace, memory_holes);
    // Sort memory and insert to trace
    set_sorted_mem_pool(&mut trace, public_input.public_memory.clone());

    trace
}

/// Returns the vector of res values.
fn compute_res(
    flags: &[CairoInstructionFlags],
    op0s: &[Felt252],
    op1s: &[Felt252],
    dsts: &[Felt252],
) -> Vec<Felt252> {
    /*
    Cairo whitepaper, page 33 - https://eprint.iacr.org/2021/1063.pdf
    # Compute res.
    if pc_update == 4:
        if res_logic == 0 && opcode == 0 && ap_update != 1:
            res = Unused
        else:
            Undefined Behavior
    else if pc_update = 0, 1 or 2:
        switch res_logic:
            case 0: res = op1
            case 1: res = op0 + op1
            case 2: res = op0 * op1
            default: Undefined Behavior
    else: Undefined Behavior
    */
    flags
        .iter()
        .zip(op0s)
        .zip(op1s)
        .zip(dsts)
        .map(|(((f, op0), op1), dst)| {
            match f.pc_update {
                PcUpdate::Jnz => {
                    match (&f.res_logic, &f.opcode, &f.ap_update) {
                        (
                            ResLogic::Op1,
                            CairoOpcode::NOp,
                            ApUpdate::Regular | ApUpdate::Add1 | ApUpdate::Add2,
                        ) => {
                            // In a `jnz` instruction, res is not used, so it is used
                            // to hold the value v = dst^(-1) as an optimization.
                            // This is important for the calculation of the `t1` virtual column
                            // values later on.
                            // See section 9.5 of the Cairo whitepaper, page 53.
                            if dst == &Felt252::zero() {
                                *dst
                            } else {
                                dst.inv().unwrap()
                            }
                        }
                        _ => {
                            panic!("Undefined Behavior");
                        }
                    }
                }
                PcUpdate::Regular | PcUpdate::Jump | PcUpdate::JumpRel => match f.res_logic {
                    ResLogic::Op1 => *op1,
                    ResLogic::Add => op0 + op1,
                    ResLogic::Mul => op0 * op1,
                    ResLogic::Unconstrained => {
                        panic!("Undefined Behavior");
                    }
                },
            }
        })
        .collect()
}

/// Returns the vector of:
/// - dst_addrs
/// - dsts
fn compute_dst(
    flags: &[CairoInstructionFlags],
    offsets: &[InstructionOffsets],
    register_states: &RegisterStates,
    memory: &CairoMemory,
) -> (Vec<Felt252>, Vec<Felt252>) {
    /* Cairo whitepaper, page 33 - https://eprint.iacr.org/2021/1063.pdf

    # Compute dst
    if dst_reg == 0:
        dst = m(ap + offdst)
    else:
        dst = m(fp + offdst)
    */
    flags
        .iter()
        .zip(offsets)
        .zip(register_states.rows.iter())
        .map(|((f, o), t)| match f.dst_reg {
            DstReg::AP => {
                let addr = t.ap.checked_add_signed(o.off_dst.into()).unwrap();
                (Felt252::from(addr), *memory.get(&addr).unwrap())
            }
            DstReg::FP => {
                let addr = t.fp.checked_add_signed(o.off_dst.into()).unwrap();
                (Felt252::from(addr), *memory.get(&addr).unwrap())
            }
        })
        .unzip()
}

/// Returns the vector of:
/// - op0_addrs
/// - op0s
fn compute_op0(
    flags: &[CairoInstructionFlags],
    offsets: &[InstructionOffsets],
    register_states: &RegisterStates,
    memory: &CairoMemory,
) -> (Vec<Felt252>, Vec<Felt252>) {
    /* Cairo whitepaper, page 33 - https://eprint.iacr.org/2021/1063.pdf

    # Compute op0.
    if op0_reg == 0:
        op0 = m(ap + offop0)
    else:
        op0 = m(fp + offop0)
    */
    flags
        .iter()
        .zip(offsets)
        .zip(register_states.rows.iter())
        .map(|((f, o), t)| match f.op0_reg {
            Op0Reg::AP => {
                let addr = t.ap.checked_add_signed(o.off_op0.into()).unwrap();
                (Felt252::from(addr), *memory.get(&addr).unwrap())
            }
            Op0Reg::FP => {
                let addr = t.fp.checked_add_signed(o.off_op0.into()).unwrap();
                (Felt252::from(addr), *memory.get(&addr).unwrap())
            }
        })
        .unzip()
}

/// Returns the vector of:
/// - op1_addrs
/// - op1s
fn compute_op1(
    flags: &[CairoInstructionFlags],
    offsets: &[InstructionOffsets],
    register_states: &RegisterStates,
    memory: &CairoMemory,
    op0s: &[Felt252],
) -> (Vec<Felt252>, Vec<Felt252>) {
    /* Cairo whitepaper, page 33 - https://eprint.iacr.org/2021/1063.pdf
    # Compute op1 and instruction_size.
    switch op1_src:
        case 0:
            instruction_size = 1
            op1 = m(op0 + offop1)
        case 1:
            instruction_size = 2
            op1 = m(pc + offop1)
            # If offop1 = 1, we have op1 = immediate_value.
        case 2:
            instruction_size = 1
            op1 = m(fp + offop1)
        case 4:
            instruction_size = 1
            op1 = m(ap + offop1)
        default:
            Undefined Behavior
    */
    flags
        .iter()
        .zip(offsets)
        .zip(op0s)
        .zip(register_states.rows.iter())
        .map(|(((flag, offset), op0), trace_state)| match flag.op1_src {
            Op1Src::Op0 => {
                let addr = aux_get_last_nim_of_field_element(op0)
                    .checked_add_signed(offset.off_op1.into())
                    .unwrap();
                (Felt252::from(addr), *memory.get(&addr).unwrap())
            }
            Op1Src::Imm => {
                let pc = trace_state.pc;
                let addr = pc.checked_add_signed(offset.off_op1.into()).unwrap();
                (Felt252::from(addr), *memory.get(&addr).unwrap())
            }
            Op1Src::AP => {
                let ap = trace_state.ap;
                let addr = ap.checked_add_signed(offset.off_op1.into()).unwrap();
                (Felt252::from(addr), *memory.get(&addr).unwrap())
            }
            Op1Src::FP => {
                let fp = trace_state.fp;
                let addr = fp.checked_add_signed(offset.off_op1.into()).unwrap();
                (Felt252::from(addr), *memory.get(&addr).unwrap())
            }
        })
        .unzip()
}

/// Depending on the instruction opcodes, some values should be updated.
/// This function updates op0s, dst, res in place when the conditions hold.
fn update_values(
    flags: &[CairoInstructionFlags],
    register_states: &RegisterStates,
    op0s: &mut [Felt252],
    dst: &mut [Felt252],
    res: &mut [Felt252],
) {
    for (i, f) in flags.iter().enumerate() {
        if f.opcode == CairoOpcode::Call {
            let instruction_size = if flags[i].op1_src == Op1Src::Imm {
                2
            } else {
                1
            };
            op0s[i] = (register_states.rows[i].pc + instruction_size).into();
            dst[i] = register_states.rows[i].fp.into();
        } else if f.opcode == CairoOpcode::AssertEq {
            res[i] = dst[i];
        }
    }
}

// NOTE: Leaving this function despite not being used anywhere. It could be useful once
// we implement layouts with the range-check builtin.
#[allow(dead_code)]
fn decompose_rc_values_into_trace_columns(rc_values: &[&Felt252]) -> [Vec<Felt252>; 8] {
    let mask = UnsignedInteger::from_hex("FFFF").unwrap();
    let mut rc_base_types: Vec<UnsignedInteger<4>> =
        rc_values.iter().map(|x| x.representative()).collect();

    let mut decomposition_columns: Vec<Vec<Felt252>> = Vec::new();

    for _ in 0..8 {
        decomposition_columns.push(
            rc_base_types
                .iter()
                .map(|&x| Felt252::from(&(x & mask)))
                .collect(),
        );

        rc_base_types = rc_base_types.iter().map(|&x| x >> 16).collect();
    }

    // This can't fail since we have 8 pushes
    decomposition_columns.try_into().unwrap()
}

// Column 1
fn set_bit_prefix_flags(trace: &mut CairoTraceTable, bit_prefix_flags: Vec<[Felt252; 16]>) {
    for (step_idx, flags) in bit_prefix_flags.into_iter().enumerate() {
        for (flag_idx, flag) in flags.into_iter().enumerate() {
            trace.set_main(flag_idx + CAIRO_STEP * step_idx, 1, flag);
        }
    }
}

// Column 0
fn set_rc_pool(
    trace: &mut CairoTraceTable,
    offsets: Vec<(Felt252, Felt252, Felt252)>,
) -> Vec<Felt252> {
    // NOTE: We should check that these offsets correspond to the off0, off1 and off2.
    const OFF_DST_OFFSET: usize = 0;
    const OFF_OP0_OFFSET: usize = 8;
    const OFF_OP1_OFFSET: usize = 4;

    let mut rc_values = Vec::new();
    for (step_idx, (off_dst, off_op0, off_op1)) in offsets.into_iter().enumerate() {
        trace.set_main(OFF_DST_OFFSET + CAIRO_STEP * step_idx, 0, off_dst);
        trace.set_main(OFF_OP0_OFFSET + CAIRO_STEP * step_idx, 0, off_op0);
        trace.set_main(OFF_OP1_OFFSET + CAIRO_STEP * step_idx, 0, off_op1);

        rc_values.push(off_dst);
        rc_values.push(off_op0);
        rc_values.push(off_op1);
    }

    rc_values
}

// Column 3
#[allow(clippy::too_many_arguments)]
fn set_mem_pool(
    trace: &mut CairoTraceTable,
    pcs: Vec<Felt252>,
    instructions: Vec<Felt252>,
    op0_addrs: Vec<Felt252>,
    op0_vals: Vec<Felt252>,
    dst_addrs: Vec<Felt252>,
    dst_vals: Vec<Felt252>,
    op1_addrs: Vec<Felt252>,
    op1_vals: Vec<Felt252>,
) -> Vec<Felt252> {
    const PC_OFFSET: usize = 0;
    const INST_OFFSET: usize = 1;
    const OP0_ADDR_OFFSET: usize = 4;
    const OP0_VAL_OFFSET: usize = 5;
    const DST_ADDR_OFFSET: usize = 8;
    const DST_VAL_OFFSET: usize = 9;
    const OP1_ADDR_OFFSET: usize = 12;
    const OP1_VAL_OFFSET: usize = 13;

    let mut addrs: Vec<Felt252> = Vec::new();
    for (step_idx, (pc, inst, op0_addr, op0_val, dst_addr, dst_val, op1_addr, op1_val)) in
        itertools::izip!(
            pcs,
            instructions,
            op0_addrs,
            op0_vals,
            dst_addrs,
            dst_vals,
            op1_addrs,
            op1_vals
        )
        .enumerate()
    {
        trace.set_main(PC_OFFSET + CAIRO_STEP * step_idx, 3, pc);
        trace.set_main(INST_OFFSET + CAIRO_STEP * step_idx, 3, inst);
        trace.set_main(OP0_ADDR_OFFSET + CAIRO_STEP * step_idx, 3, op0_addr);
        trace.set_main(OP0_VAL_OFFSET + CAIRO_STEP * step_idx, 3, op0_val);
        trace.set_main(DST_ADDR_OFFSET + CAIRO_STEP * step_idx, 3, dst_addr);
        trace.set_main(DST_VAL_OFFSET + CAIRO_STEP * step_idx, 3, dst_val);
        trace.set_main(OP1_ADDR_OFFSET + CAIRO_STEP * step_idx, 3, op1_addr);
        trace.set_main(OP1_VAL_OFFSET + CAIRO_STEP * step_idx, 3, op1_val);

        addrs.push(pc);
        addrs.push(op0_addr);
        addrs.push(dst_addr);
        addrs.push(op1_addr);
    }

    addrs
}

// Column 5
fn set_update_pc(
    trace: &mut CairoTraceTable,
    aps: Vec<Felt252>,
    t0s: Vec<Felt252>,
    t1s: Vec<Felt252>,
    mul: Vec<Felt252>,
    fps: Vec<Felt252>,
    res: Vec<Felt252>,
) {
    const AP_OFFSET: usize = 0;
    const TMP0_OFFSET: usize = 2;
    const OPS_MUL_OFFSET: usize = 4;
    const FP_OFFSET: usize = 8;
    const TMP1_OFFSET: usize = 10;
    const RES_OFFSET: usize = 12;

    for (step_idx, (ap, tmp0, m, fp, tmp1, res)) in
        itertools::izip!(aps, t0s, mul, fps, t1s, res).enumerate()
    {
        trace.set_main(AP_OFFSET + CAIRO_STEP * step_idx, 5, ap);
        trace.set_main(TMP0_OFFSET + CAIRO_STEP * step_idx, 5, tmp0);
        trace.set_main(OPS_MUL_OFFSET + CAIRO_STEP * step_idx, 5, m);
        trace.set_main(FP_OFFSET + CAIRO_STEP * step_idx, 5, fp);
        trace.set_main(TMP1_OFFSET + CAIRO_STEP * step_idx, 5, tmp1);
        trace.set_main(RES_OFFSET + CAIRO_STEP * step_idx, 5, res);
    }
}

fn finalize_mem_pool(trace: &mut CairoTraceTable, memory_holes: VecDeque<Felt252>) {
    const MEM_POOL_UNUSED_ADDR_OFFSET: usize = 6;
    const MEM_POOL_UNUSED_VALUE_OFFSET: usize = 7;
    const MEM_POOL_UNUSED_CELL_STEP: usize = 8;

    let mut memory_holes = memory_holes;

    let last_hole_addr = memory_holes.pop_back().unwrap();

    for step_idx in 0..trace.num_steps() {
        if let Some(hole_addr) = memory_holes.pop_front() {
            trace.set_main(
                MEM_POOL_UNUSED_ADDR_OFFSET + CAIRO_STEP * step_idx,
                3,
                hole_addr,
            );
            trace.set_main(
                MEM_POOL_UNUSED_VALUE_OFFSET + CAIRO_STEP * step_idx,
                3,
                Felt252::zero(),
            );
        } else {
            trace.set_main(
                MEM_POOL_UNUSED_ADDR_OFFSET + CAIRO_STEP * step_idx,
                3,
                last_hole_addr,
            );
            trace.set_main(
                MEM_POOL_UNUSED_VALUE_OFFSET + CAIRO_STEP * step_idx,
                3,
                Felt252::zero(),
            );
        }

        if let Some(hole_addr) = memory_holes.pop_front() {
            trace.set_main(
                MEM_POOL_UNUSED_ADDR_OFFSET + MEM_POOL_UNUSED_CELL_STEP + CAIRO_STEP * step_idx,
                3,
                hole_addr,
            );
            trace.set_main(
                MEM_POOL_UNUSED_VALUE_OFFSET + MEM_POOL_UNUSED_CELL_STEP + CAIRO_STEP * step_idx,
                3,
                Felt252::zero(),
            );
        } else {
            trace.set_main(
                MEM_POOL_UNUSED_ADDR_OFFSET + MEM_POOL_UNUSED_CELL_STEP + CAIRO_STEP * step_idx,
                3,
                last_hole_addr,
            );
            trace.set_main(
                MEM_POOL_UNUSED_VALUE_OFFSET + MEM_POOL_UNUSED_CELL_STEP + CAIRO_STEP * step_idx,
                3,
                Felt252::zero(),
            );
        }

        assert!(memory_holes.is_empty());
    }
}

fn set_sorted_rc_pool(trace: &mut CairoTraceTable, sorted_rc_column: Vec<Felt252>) {
    for (row_idx, rc_value) in sorted_rc_column.into_iter().enumerate() {
        trace.set_main(row_idx, 2, rc_value);
    }
}

fn finalize_rc_pool(trace: &mut CairoTraceTable, rc_holes: VecDeque<Felt252>, rc_max: Felt252) {
    let mut rc_holes = rc_holes;

    let reserved_cell_idxs = [4, 8];
    for step_idx in 0..trace.num_steps() {
        for step_cell_idx in 1..CAIRO_STEP {
            if reserved_cell_idxs.contains(&step_cell_idx) {
                continue;
            };
            if let Some(rc_hole) = rc_holes.pop_front() {
                trace.set_main(step_idx * CAIRO_STEP + step_cell_idx, 0, rc_hole);
            } else {
                trace.set_main(step_idx * CAIRO_STEP + step_cell_idx, 0, rc_max);
            }
        }
    }

    assert!(rc_holes.is_empty());
}

fn set_sorted_mem_pool(trace: &mut CairoTraceTable, pub_memory: HashMap<Felt252, Felt252>) {
    const PUB_MEMORY_ADDR_OFFSET: usize = 2;
    const PUB_MEMORY_VALUE_OFFSET: usize = 3;
    const PUB_MEMORY_STEP: usize = 8;

    assert!(2 * trace.num_steps() >= pub_memory.len());

    let mut mem_pool = trace.get_column_main(3);
    let first_pub_memory_addr = Felt252::one();
    let first_pub_memory_value = *pub_memory.get(&first_pub_memory_addr).unwrap();
    let first_pub_memory_entry_padding_len = 2 * trace.num_steps() - pub_memory.len();

    let padding_num_steps = first_pub_memory_entry_padding_len.div_ceil(2);
    let mut first_addr_padding =
        iter::repeat(first_pub_memory_addr).take(first_pub_memory_entry_padding_len);
    let mut padding_step_flag = 0;

    // this loop is for padding with (first_pub_addr, first_pub_value)
    for step_idx in 0..padding_num_steps {
        if let Some(first_addr) = first_addr_padding.next() {
            mem_pool[PUB_MEMORY_ADDR_OFFSET + CAIRO_STEP * step_idx] = first_addr;
            mem_pool[PUB_MEMORY_VALUE_OFFSET + CAIRO_STEP * step_idx] = first_pub_memory_value;
        } else {
            padding_step_flag = 0;
            break;
        }
        if let Some(first_addr) = first_addr_padding.next() {
            mem_pool[PUB_MEMORY_STEP + PUB_MEMORY_ADDR_OFFSET + CAIRO_STEP * step_idx] = first_addr;
            mem_pool[PUB_MEMORY_STEP + PUB_MEMORY_VALUE_OFFSET + CAIRO_STEP * step_idx] =
                first_pub_memory_value;
        } else {
            padding_step_flag = 1;
            break;
        }
    }

    let mut pub_memory_iter = pub_memory.iter();
    if padding_step_flag == 1 {
        let (pub_memory_addr, pub_memory_value) = pub_memory_iter.next().unwrap();
        mem_pool[PUB_MEMORY_STEP + PUB_MEMORY_ADDR_OFFSET + CAIRO_STEP * (padding_num_steps - 1)] =
            *pub_memory_addr;
        mem_pool
            [PUB_MEMORY_STEP + PUB_MEMORY_VALUE_OFFSET + CAIRO_STEP * (padding_num_steps - 1)] =
            *pub_memory_value;
    }

    for step_idx in padding_num_steps..trace.num_steps() {
        let (pub_memory_addr, pub_memory_value) = pub_memory_iter.next().unwrap();
        mem_pool[PUB_MEMORY_ADDR_OFFSET + CAIRO_STEP * step_idx] = *pub_memory_addr;
        mem_pool[PUB_MEMORY_VALUE_OFFSET + CAIRO_STEP * step_idx] = *pub_memory_value;

        let (pub_memory_addr, pub_memory_value) = pub_memory_iter.next().unwrap();
        mem_pool[PUB_MEMORY_STEP + PUB_MEMORY_ADDR_OFFSET + CAIRO_STEP * step_idx] =
            *pub_memory_addr;
        mem_pool[PUB_MEMORY_STEP + PUB_MEMORY_VALUE_OFFSET + CAIRO_STEP * step_idx] =
            *pub_memory_value;
    }

    let mut addrs = Vec::with_capacity(trace.num_rows() / 2);
    let mut values = Vec::with_capacity(trace.num_rows() / 2);
    let mut sorted_addrs = Vec::with_capacity(trace.num_rows() / 2);
    let mut sorted_values = Vec::with_capacity(trace.num_rows() / 2);
    for (idx, mem_cell) in mem_pool.into_iter().enumerate() {
        if idx % 2 == 0 {
            let addr = mem_cell;
            addrs.push(addr);
        } else {
            let value = mem_cell;
            values.push(value);
        }
    }

    let mut sorted_addr_idxs: Vec<usize> = (0..addrs.len()).collect();
    sorted_addr_idxs.sort_by_key(|&idx| addrs[idx]);
    for idx in sorted_addr_idxs.iter() {
        sorted_addrs.push(addrs[*idx]);
        sorted_values.push(values[*idx]);
    }

    let mut sorted_addrs_iter = sorted_addrs.into_iter();
    let mut sorted_values_iter = sorted_values.into_iter();
    for row_idx in 0..trace.num_rows() {
        if row_idx % 2 == 0 {
            let addr = sorted_addrs_iter.next().unwrap();
            trace.set_main(row_idx, 4, addr);
        } else {
            let value = sorted_values_iter.next().unwrap();
            trace.set_main(row_idx, 4, value);
        }
    }
}

pub(crate) fn set_rc_permutation_column(trace: &mut CairoTraceTable, z: &Felt252) {
    let mut denominator_evaluations = trace
        .get_column_main(2)
        .iter()
        .map(|a_prime| z - a_prime)
        .collect_vec();
    FieldElement::inplace_batch_inverse(&mut denominator_evaluations).unwrap();

    let rc_cumulative_procuts = trace
        .get_column_main(0)
        .iter()
        .zip(&denominator_evaluations)
        .scan(Felt252::one(), |product, (num_i, den_i)| {
            let ret = *product;
            *product = ret * (z - num_i) * den_i;
            Some(*product)
        })
        .collect_vec();

    for (i, rc_perm_i) in rc_cumulative_procuts.into_iter().enumerate() {
        trace.set_aux(i, 0, rc_perm_i)
    }
}

pub(crate) fn set_mem_permutation_column(
    trace: &mut CairoTraceTable,
    alpha_mem: &Felt252,
    z_mem: &Felt252,
) {
    let sorted_mem_pool = trace.get_column_main(4);
    let sorted_addrs = sorted_mem_pool.iter().step_by(2).collect_vec();
    let sorted_values = sorted_mem_pool[1..].iter().step_by(2).collect_vec();

    let mut denominator = std::iter::zip(sorted_addrs, sorted_values)
        .map(|(ap, vp)| z_mem - (ap + alpha_mem * vp))
        .collect_vec();
    FieldElement::inplace_batch_inverse(&mut denominator).unwrap();

    let mem_pool = trace.get_column_main(3);
    let addrs = mem_pool.iter().step_by(2).collect_vec();
    let values = mem_pool[1..].iter().step_by(2).collect_vec();

    let mem_cumulative_products = itertools::izip!(addrs, values, denominator)
        .scan(Felt252::one(), |product, (a_i, v_i, den_i)| {
            let ret = *product;
            *product = ret * ((z_mem - (a_i + alpha_mem * v_i)) * den_i);
            Some(*product)
        })
        .collect_vec();

    for (i, row_idx) in (0..trace.num_rows()).step_by(2).enumerate() {
        let mem_cumul_prod = mem_cumulative_products[i];
        trace.set_aux(row_idx, 1, mem_cumul_prod);
    }
}

#[cfg(test)]
mod test {

    use super::*;
    use crate::{
        cairo_layout::CairoLayout, runner::run::run_program, tests::utils::cairo0_program_path,
    };

    #[test]
    fn test_rc_decompose() {
        let fifteen = Felt252::from_hex("000F000F000F000F000F000F000F000F").unwrap();
        let sixteen = Felt252::from_hex("00100010001000100010001000100010").unwrap();
        let one_two_three = Felt252::from_hex("00010002000300040005000600070008").unwrap();

        let decomposition_columns =
            decompose_rc_values_into_trace_columns(&[&fifteen, &sixteen, &one_two_three]);

        for row in &decomposition_columns {
            assert_eq!(row[0], Felt252::from_hex("F").unwrap());
            assert_eq!(row[1], Felt252::from_hex("10").unwrap());
        }

        assert_eq!(decomposition_columns[0][2], Felt252::from_hex("8").unwrap());
        assert_eq!(decomposition_columns[1][2], Felt252::from_hex("7").unwrap());
        assert_eq!(decomposition_columns[2][2], Felt252::from_hex("6").unwrap());
        assert_eq!(decomposition_columns[3][2], Felt252::from_hex("5").unwrap());
        assert_eq!(decomposition_columns[4][2], Felt252::from_hex("4").unwrap());
        assert_eq!(decomposition_columns[5][2], Felt252::from_hex("3").unwrap());
        assert_eq!(decomposition_columns[6][2], Felt252::from_hex("2").unwrap());
        assert_eq!(decomposition_columns[7][2], Felt252::from_hex("1").unwrap());
    }

    #[test]
    fn test_get_memory_holes_empty_pub_memory() {
        // We construct a sorted addresses list [1, 2, 3, 6, 7, 8, 9, 13, 14, 15], and
        // an empty public memory. This way, any holes present between
        // the min and max addresses should be returned by the function.
        // NOTE: The memory hole at address 16 will also be returned because the max addr
        // +1 is considered a memory hole too.
        let mut addrs: Vec<Felt252> = (1..4).map(Felt252::from).collect();
        let addrs_extension: Vec<Felt252> = (6..10).map(Felt252::from).collect();
        addrs.extend_from_slice(&addrs_extension);
        let addrs_extension: Vec<Felt252> = (13..16).map(Felt252::from).collect();
        addrs.extend_from_slice(&addrs_extension);
        let pub_memory = HashMap::new();

        let expected_memory_holes = VecDeque::from([
            Felt252::from(4),
            Felt252::from(5),
            Felt252::from(10),
            Felt252::from(11),
            Felt252::from(12),
            Felt252::from(16),
        ]);
        let calculated_memory_holes = get_memory_holes(&addrs, &pub_memory);

        assert_eq!(expected_memory_holes, calculated_memory_holes);
    }

    #[test]
    fn test_get_memory_holes_inside_program_section() {
        // We construct a sorted addresses list [1, 2, 3, 8, 9] and we
        // set public memory from address 1 to 9. Since all the holes will be inside the
        // program segment (meaning from addresses 1 to 9), the function
        // should not return any of them.
        let mut addrs: Vec<Felt252> = (1..4).map(Felt252::from).collect();
        let addrs_extension: Vec<Felt252> = (8..10).map(Felt252::from).collect();
        addrs.extend_from_slice(&addrs_extension);

        let mut pub_memory = HashMap::new();
        for addr in 1..=9 {
            let addr = Felt252::from(addr);
            pub_memory.insert(addr, Felt252::zero());
        }

        let calculated_memory_holes = get_memory_holes(&addrs, &pub_memory);

        // max_addr + 1 (10, in this case) is always returned by the get_memory_holes function
        let expected_memory_holes = VecDeque::from([Felt252::from(10)]);

        assert_eq!(expected_memory_holes, calculated_memory_holes);
    }

    #[test]
    fn test_get_memory_holes_outside_program_section() {
        // We construct a sorted addresses list [1, 2, 3, 8, 9] and we
        // set public memory from addresses 1 to 6. The holes found inside the program section,
        // i.e. in the address range between 1 to 6, should not be returned.
        // So addresses 4, 5 and 6 will no be returned, only address 7.
        let mut addrs: Vec<Felt252> = (1..4).map(Felt252::from).collect();
        let addrs_extension: Vec<Felt252> = (8..10).map(Felt252::from).collect();
        addrs.extend_from_slice(&addrs_extension);

        let mut pub_memory = HashMap::new();
        for addr in 0..=6 {
            let addr = Felt252::from(addr);
            pub_memory.insert(addr, Felt252::zero());
        }

        let calculated_memory_holes = get_memory_holes(&addrs, &pub_memory);
        let expected_memory_holes = VecDeque::from([Felt252::from(7), Felt252::from(10)]);

        assert_eq!(expected_memory_holes, calculated_memory_holes);
    }

    #[test]
    fn set_rc_pool_works() {
        let program_content = std::fs::read(cairo0_program_path("fibonacci_stone.json")).unwrap();
        let mut trace: CairoTraceTable = TraceTable::allocate_with_zeros(128, 6, 2, 16);
        let (register_states, memory, _) =
            run_program(None, CairoLayout::Plain, &program_content).unwrap();

        let (_, biased_offsets): (Vec<CairoInstructionFlags>, Vec<InstructionOffsets>) =
            register_states
                .flags_and_offsets(&memory)
                .unwrap()
                .into_iter()
                .unzip();

        let unbiased_offsets: Vec<(Felt252, Felt252, Felt252)> = biased_offsets
            .iter()
            .map(InstructionOffsets::to_trace_representation)
            .collect();

        let rc_values = set_rc_pool(&mut trace, unbiased_offsets);
        let mut sorted_rc_values: Vec<u16> = rc_values
            .iter()
            .map(|x| x.representative().into())
            .collect();
        sorted_rc_values.sort();
        let rc_holes = get_rc_holes(&sorted_rc_values);
        let rc_max = Felt252::from(*(sorted_rc_values.last().unwrap()) as u64);
        finalize_rc_pool(&mut trace, rc_holes, rc_max);

        let mut sorted_rc_column = trace.get_column_main(0);
        sorted_rc_column.sort_by_key(|x| x.representative());
        set_sorted_rc_pool(&mut trace, sorted_rc_column);

        trace.main_table.columns()[2]
            .iter()
            .enumerate()
            .for_each(|(i, v)| println!("SORTED RC VAL {} - {}", i, v));
    }

    #[test]
    fn set_update_pc_works() {
        let program_content = std::fs::read(cairo0_program_path("fibonacci_stone.json")).unwrap();
        let mut trace: CairoTraceTable = TraceTable::allocate_with_zeros(128, 6, 2, 16);
        let (register_states, memory, _) =
            run_program(None, CairoLayout::Plain, &program_content).unwrap();

        let (flags, biased_offsets): (Vec<CairoInstructionFlags>, Vec<InstructionOffsets>) =
            register_states
                .flags_and_offsets(&memory)
                .unwrap()
                .into_iter()
                .unzip();

        // dst, op0, op1 and res are computed from flags and offsets
        let (_dst_addrs, mut dsts): (Vec<Felt252>, Vec<Felt252>) =
            compute_dst(&flags, &biased_offsets, &register_states, &memory);
        let (_op0_addrs, mut op0s): (Vec<Felt252>, Vec<Felt252>) =
            compute_op0(&flags, &biased_offsets, &register_states, &memory);
        let (_op1_addrs, op1s): (Vec<Felt252>, Vec<Felt252>) =
            compute_op1(&flags, &biased_offsets, &register_states, &memory, &op0s);
        let mut res = compute_res(&flags, &op0s, &op1s, &dsts);

        update_values(&flags, &register_states, &mut op0s, &mut dsts, &mut res);

        let aps: Vec<Felt252> = register_states
            .rows
            .iter()
            .map(|t| Felt252::from(t.ap))
            .collect();
        let fps: Vec<Felt252> = register_states
            .rows
            .iter()
            .map(|t| Felt252::from(t.fp))
            .collect();

        let trace_repr_flags: Vec<[Felt252; 16]> = flags
            .iter()
            .map(CairoInstructionFlags::to_trace_representation)
            .collect();

        let two = Felt252::from(2);
        let t0: Vec<Felt252> = trace_repr_flags
            .iter()
            .zip(&dsts)
            .map(|(repr_flags, dst)| (repr_flags[9] - two * repr_flags[10]) * dst)
            .collect();
        let t1: Vec<Felt252> = t0.iter().zip(&res).map(|(t, r)| t * r).collect();
        let mul: Vec<Felt252> = op0s.iter().zip(&op1s).map(|(op0, op1)| op0 * op1).collect();

        set_update_pc(&mut trace, aps, t0, t1, mul, fps, res);

        trace.main_table.columns()[5][0..50]
            .iter()
            .enumerate()
            .for_each(|(i, v)| println!("ROW {} - VALUE: {}", i, v));
    }

    #[test]
    fn set_mem_pool_works() {
        let program_content = std::fs::read(cairo0_program_path("fibonacci_stone.json")).unwrap();
        let mut trace: CairoTraceTable = TraceTable::allocate_with_zeros(128, 6, 2, 16);
        let (register_states, memory, pub_inputs) =
            run_program(None, CairoLayout::Plain, &program_content).unwrap();

        let (flags, biased_offsets): (Vec<CairoInstructionFlags>, Vec<InstructionOffsets>) =
            register_states
                .flags_and_offsets(&memory)
                .unwrap()
                .into_iter()
                .unzip();

        // dst, op0, op1 and res are computed from flags and offsets
        let (dst_addrs, mut dsts): (Vec<Felt252>, Vec<Felt252>) =
            compute_dst(&flags, &biased_offsets, &register_states, &memory);
        let (op0_addrs, mut op0s): (Vec<Felt252>, Vec<Felt252>) =
            compute_op0(&flags, &biased_offsets, &register_states, &memory);
        let (op1_addrs, op1s): (Vec<Felt252>, Vec<Felt252>) =
            compute_op1(&flags, &biased_offsets, &register_states, &memory, &op0s);
        let mut res = compute_res(&flags, &op0s, &op1s, &dsts);

        update_values(&flags, &register_states, &mut op0s, &mut dsts, &mut res);

        let pcs: Vec<Felt252> = register_states
            .rows
            .iter()
            .map(|t| Felt252::from(t.pc))
            .collect();
        let instructions: Vec<Felt252> = register_states
            .rows
            .iter()
            .map(|t| *memory.get(&t.pc).unwrap())
            .collect();

        let mut sorted_addrs = set_mem_pool(
            &mut trace,
            pcs,
            instructions,
            op0_addrs,
            op0s,
            dst_addrs,
            dsts,
            op1_addrs,
            op1s,
        );

        sorted_addrs.sort_by_key(|x| x.representative());
        let memory_holes = get_memory_holes(&sorted_addrs, &pub_inputs.public_memory);
        finalize_mem_pool(&mut trace, memory_holes);

        set_sorted_mem_pool(&mut trace, pub_inputs.public_memory);

        let z = Felt252::from_hex_unchecked(
            "0x6896a2e62f03d4d1f625efb97468ef93f31105bb51a83d550bca6fdebd035de",
        );
        let alpha = Felt252::from_hex_unchecked(
            "0x64de8f5be59594e112d438c13ec4916e138b013e7d388b681c11b03ede7962e",
        );
        set_mem_permutation_column(&mut trace, &alpha, &z);

        trace.aux_table.columns()[1]
            .iter()
            .enumerate()
            .for_each(|(i, v)| println!("ROW {} - MEM CUMUL PROD: {}", i, v));
    }

    #[test]
    fn set_bit_prefix_flags_works() {
        let program_content = std::fs::read(cairo0_program_path("fibonacci_stone.json")).unwrap();
        let mut trace: CairoTraceTable = TraceTable::allocate_with_zeros(128, 6, 2, 16);
        let (register_states, memory, _) =
            run_program(None, CairoLayout::Plain, &program_content).unwrap();

        let (flags, _biased_offsets): (Vec<CairoInstructionFlags>, Vec<InstructionOffsets>) =
            register_states
                .flags_and_offsets(&memory)
                .unwrap()
                .into_iter()
                .unzip();

        let bit_prefix_flags: Vec<[Felt252; 16]> = flags
            .iter()
            .map(CairoInstructionFlags::to_trace_representation)
            .collect();

        set_bit_prefix_flags(&mut trace, bit_prefix_flags);

        trace.main_table.columns()[1][0..50]
            .iter()
            .enumerate()
            .for_each(|(i, v)| println!("ROW {} - VALUE: {}", i, v));
    }

    #[test]
    fn set_rc_permutation_col_works() {
        let program_content = std::fs::read(cairo0_program_path("fibonacci_stone.json")).unwrap();
        let mut trace: CairoTraceTable = TraceTable::allocate_with_zeros(128, 6, 2, 16);
        let (register_states, memory, _) =
            run_program(None, CairoLayout::Plain, &program_content).unwrap();

        let (_, biased_offsets): (Vec<CairoInstructionFlags>, Vec<InstructionOffsets>) =
            register_states
                .flags_and_offsets(&memory)
                .unwrap()
                .into_iter()
                .unzip();

        let unbiased_offsets: Vec<(Felt252, Felt252, Felt252)> = biased_offsets
            .iter()
            .map(InstructionOffsets::to_trace_representation)
            .collect();

        let rc_values = set_rc_pool(&mut trace, unbiased_offsets);
        let mut sorted_rc_values: Vec<u16> = rc_values
            .iter()
            .map(|x| x.representative().into())
            .collect();
        sorted_rc_values.sort();

        let rc_holes = get_rc_holes(&sorted_rc_values);
        let rc_max = Felt252::from(*(sorted_rc_values.last().unwrap()) as u64);
        finalize_rc_pool(&mut trace, rc_holes, rc_max);

        let mut sorted_rc_column = trace.get_column_main(0);
        sorted_rc_column.sort_by_key(|x| x.representative());
        set_sorted_rc_pool(&mut trace, sorted_rc_column);

        let z = Felt252::from_hex_unchecked(
            "0x221ee7f99bdf1f11e16445f06fd90f413146e1764a1d16d46525148456cc3eb",
        );

        set_rc_permutation_column(&mut trace, &z);

        trace.aux_table.columns()[0][..20]
            .iter()
            .enumerate()
            .for_each(|(i, v)| println!("RC PERMUTATION ARG {} - {}", i, v));
    }
}
