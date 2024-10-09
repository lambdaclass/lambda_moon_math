use crate::{
    fft::cpu::roots_of_unity::get_powers_of_primitive_root,
    field::{
        element::FieldElement,
        traits::{IsFFTField, RootsConfig},
    },
};
use alloc::vec::Vec;

/// Calculates the (non-unitary) Discrete Fourier Transform of `input` via the DFT matrix.
pub fn naive_matrix_dft_test<F: IsFFTField>(input: &[FieldElement<F>]) -> Vec<FieldElement<F>> {
    let n = input.len();
    assert!(n.is_power_of_two());
    let order = n.trailing_zeros();

    let twiddles =
        get_powers_of_primitive_root::<F>(order.into(), n, RootsConfig::Natural).unwrap();

    let mut output = Vec::with_capacity(n);
    for row in 0..n {
        let mut sum = FieldElement::zero();

        for (col, element) in input.iter().enumerate() {
            let i = (row * col) % n; // w^i = w^(i mod n)
            sum += element.clone() * twiddles[i].clone();
        }

        output.push(sum);
    }

    output
}

#[cfg(test)]
mod fft_helpers_test {
    use crate::{
        fft::{
            cpu::roots_of_unity::get_powers_of_primitive_root, test_helpers::naive_matrix_dft_test,
        },
        field::{
            element::FieldElement, test_fields::u64_test_field::U64TestField, traits::RootsConfig,
        },
        polynomial::Polynomial,
    };
    use proptest::{collection, prelude::*};

    type F = U64TestField;
    type FE = FieldElement<F>;

    prop_compose! {
        fn powers_of_two(max_exp: u8)(exp in 1..max_exp) -> usize { 1 << exp }
        // max_exp cannot be multiple of the bits that represent a usize, generally 64 or 32.
        // also it can't exceed the test field's two-adicity.
    }
    prop_compose! {
        fn field_element()(num in any::<u64>().prop_filter("Avoid null coefficients", |x| x != &0)) -> FE {
            FE::from(num)
        }
    }
    prop_compose! {
        fn field_vec(max_exp: u8)(vec in collection::vec(field_element(), 2..1<<max_exp).prop_filter("Avoid polynomials of size not power of two", |vec| vec.len().is_power_of_two())) -> Vec<FE> {
            vec
        }
    }

    proptest! {
        // Property-based test that ensures dft() gives the same result as a naive polynomial evaluation.
        #[test]
        fn test_dft_same_as_eval(coeffs in field_vec(8)) {
            let dft = naive_matrix_dft_test(&coeffs);

            let poly = Polynomial::new(&coeffs);
            let order = coeffs.len().trailing_zeros();
            let twiddles = get_powers_of_primitive_root(order.into(), coeffs.len(), RootsConfig::Natural).unwrap();
            let evals: Vec<FE> = twiddles.iter().map(|x| poly.evaluate(x)).collect();

            prop_assert_eq!(evals, dft);
        }
    }
}
