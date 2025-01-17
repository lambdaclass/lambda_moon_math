use lambdaworks_math::{
    fft::cpu::roots_of_unity::get_powers_of_primitive_root_coset,
    field::{element::FieldElement, traits::IsFFTField},
};

use super::traits::AIR;

pub struct Domain<F: IsFFTField> {
    pub(crate) root_order: u32,
    pub(crate) lde_roots_of_unity_coset: Vec<FieldElement<F>>,
    pub(crate) trace_primitive_root: FieldElement<F>,
    pub(crate) trace_roots_of_unity: Vec<FieldElement<F>>,
    pub(crate) coset_offset: FieldElement<F>,
    pub(crate) blowup_factor: usize,
    pub(crate) interpolation_domain_size: usize,
}

impl<F: IsFFTField> Domain<F> {
    pub fn new<A>(air: &A) -> Self
    where
        A: AIR<Field = F>,
    {
        // Initial definitions
        let blowup_factor = air.options().blowup_factor as usize;
        let coset_offset = FieldElement::from(air.options().coset_offset);
        let interpolation_domain_size = air.trace_length();
        let root_order = air.trace_length().trailing_zeros();
        // * Generate Coset
        let trace_primitive_root = F::get_primitive_root_of_unity(root_order as u64).unwrap();
        let trace_roots_of_unity = get_powers_of_primitive_root_coset(
            root_order as u64,
            interpolation_domain_size,
            &FieldElement::one(),
        )
        .unwrap();

        let lde_root_order = (air.trace_length() * blowup_factor).trailing_zeros();
        let lde_roots_of_unity_coset = get_powers_of_primitive_root_coset(
            lde_root_order as u64,
            air.trace_length() * blowup_factor,
            &coset_offset,
        )
        .unwrap();

        Self {
            root_order,
            lde_roots_of_unity_coset,
            trace_primitive_root,
            trace_roots_of_unity,
            blowup_factor,
            coset_offset,
            interpolation_domain_size,
        }
    }
}
