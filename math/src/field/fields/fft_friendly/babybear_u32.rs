use crate::field::{
    fields::u32_montgomery_backend_prime_field::U32MontgomeryBackendPrimeField, traits::IsFFTField,
};

// Babybear Prime p = 2^31 - 2^27 + 1 = 0x78000001 = 2013265921
pub type Babybear31PrimeField = U32MontgomeryBackendPrimeField<2013265921>;

// p = 2^31 - 2^27 + 1 = 2^27 * (2^4-1) + 1, then
// there is a gruop in the field of order 2^27.
// Since we want to have margin to be able to define a bigger group (blow-up group),
// we define TWO_ADICITY as 24 (so the blow-up factor can be 2^3 = 8).
// A two-adic primitive root of unity is 21^(2^24) because
// 21^(2^24)=1 mod 2013265921.
// In the future we should allow this with metal and cuda feature, and just dispatch it to the CPU until the implementation is done
#[cfg(any(not(feature = "metal"), not(feature = "cuda")))]
impl IsFFTField for Babybear31PrimeField {
    const TWO_ADICITY: u64 = 24;

    const TWO_ADIC_PRIMITVE_ROOT_OF_UNITY: Self::BaseType = 21;

    fn field_name() -> &'static str {
        "babybear31"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    mod test_babybear_31_ops {
        use super::*;
        use crate::{
            errors::CreationError,
            field::{element::FieldElement, errors::FieldError, traits::IsPrimeField},
            traits::ByteConversion,
        };
        type FE = FieldElement<Babybear31PrimeField>;

        #[test]
        fn two_plus_one_is_three() {
            let a = FE::from(2);
            let b = FE::one();
            let res = FE::from(3);

            assert_eq!(a + b, res)
        }

        #[test]
        fn one_minus_two_is_minus_one() {
            let a = FE::from(2);
            let b = FE::one();
            let res = FE::from(2013265920);
            assert_eq!(b - a, res)
        }

        #[test]
        fn mul_by_zero_is_zero() {
            let a = FE::from(2);
            let b = FE::zero();
            assert_eq!(a * b, b)
        }

        #[test]
        fn neg_zero_is_zero() {
            let zero = FE::from(0);

            assert_eq!(-&zero, zero);
        }

        #[test]
        fn doubling() {
            assert_eq!(FE::from(2).double(), FE::from(2) + FE::from(2),);
        }

        const ORDER: usize = 2013265921;

        #[test]
        fn order_is_0() {
            assert_eq!(FE::from((ORDER - 1) as u64) + FE::from(1), FE::from(0));
        }

        #[test]
        fn when_comparing_13_and_13_they_are_equal() {
            let a: FE = FE::from(13);
            let b: FE = FE::from(13);
            assert_eq!(a, b);
        }

        #[test]
        fn when_comparing_13_and_8_they_are_different() {
            let a: FE = FE::from(13);
            let b: FE = FE::from(8);
            assert_ne!(a, b);
        }

        #[test]
        fn mul_neutral_element() {
            let a: FE = FE::from(1);
            let b: FE = FE::from(2);
            assert_eq!(a * b, FE::from(2));
        }

        #[test]
        fn mul_2_3_is_6() {
            let a: FE = FE::from(2);
            let b: FE = FE::from(3);
            assert_eq!(a * b, FE::from(6));
        }

        #[test]
        fn mul_order_minus_1() {
            let a: FE = FE::from((ORDER - 1) as u64);
            let b: FE = FE::from((ORDER - 1) as u64);
            assert_eq!(a * b, FE::from(1));
        }

        #[test]
        fn inv_0_error() {
            let result = FE::from(0).inv();
            assert!(matches!(result, Err(FieldError::InvZeroError)))
        }

        #[test]
        fn inv_2_mul_2_is_1() {
            let a: FE = FE::from(2);
            assert_eq!(a * a.inv().unwrap(), FE::from(1));
        }

        #[test]
        fn square_2_is_4() {
            assert_eq!(FE::from(2).square(), FE::from(4))
        }

        #[test]
        fn pow_2_3_is_8() {
            assert_eq!(FE::from(2).pow(3_u64), FE::from(8))
        }

        #[test]
        fn pow_p_minus_1() {
            assert_eq!(FE::from(2).pow(ORDER - 1), FE::from(1))
        }

        #[test]
        fn div_1() {
            assert_eq!(FE::from(2) / FE::from(1), FE::from(2))
        }

        #[test]
        fn div_4_2() {
            assert_eq!(FE::from(4) / FE::from(2), FE::from(2))
        }

        #[test]
        fn two_plus_its_additive_inv_is_0() {
            let two = FE::from(2);

            assert_eq!(two + (-&two), FE::from(0))
        }

        #[test]
        fn four_minus_three_is_1() {
            let four = FE::from(4);
            let three = FE::from(3);

            assert_eq!(four - three, FE::from(1))
        }

        #[test]
        fn zero_minus_1_is_order_minus_1() {
            let zero = FE::from(0);
            let one = FE::from(1);

            assert_eq!(zero - one, FE::from((ORDER - 1) as u64))
        }

        #[test]
        fn babybear_uses_31_bits() {
            assert_eq!(Babybear31PrimeField::field_bit_size(), 31);
        }

        #[test]
        fn montgomery_backend_prime_field_compute_mu_parameter() {
            let mu_expected: u32 = 2281701377;
            assert_eq!(Babybear31PrimeField::MU, mu_expected);
        }

        #[test]
        fn montgomery_backend_prime_field_compute_r2_parameter() {
            let r2_expected: u32 = 1172168163;
            assert_eq!(Babybear31PrimeField::R2, r2_expected);
        }

        #[test]
        #[cfg(feature = "alloc")]
        fn from_hex_bigger_than_u64_returns_error() {
            let x = FE::from_hex("5f103b0bd4397d4df560eb559f38353f80eeb6");
            assert!(matches!(x, Err(CreationError::InvalidHexString)))
        }

        #[test]
        #[cfg(feature = "alloc")]
        fn to_bytes_from_bytes_be_is_the_identity() {
            let x = FE::from_hex("5f103b").unwrap();
            assert_eq!(FE::from_bytes_be(&x.to_bytes_be()).unwrap(), x);
        }

        #[test]
        #[cfg(feature = "alloc")]
        fn from_bytes_to_bytes_be_is_the_identity() {
            let bytes = [0, 0, 0, 1];
            assert_eq!(FE::from_bytes_be(&bytes).unwrap().to_bytes_be(), bytes);
        }

        #[test]
        #[cfg(feature = "alloc")]
        fn to_bytes_from_bytes_le_is_the_identity() {
            let x = FE::from_hex("5f103b").unwrap();
            assert_eq!(FE::from_bytes_le(&x.to_bytes_le()).unwrap(), x);
        }

        #[test]
        #[cfg(feature = "alloc")]
        fn from_bytes_to_bytes_le_is_the_identity_4_bytes() {
            let bytes = [1, 0, 0, 0];
            assert_eq!(FE::from_bytes_le(&bytes).unwrap().to_bytes_le(), bytes);
        }

        #[test]
        #[cfg(feature = "alloc")]
        fn byte_serialization_for_a_number_matches_with_byte_conversion_implementation_le() {
            let element = FE::from_hex("0123456701234567").unwrap();
            let bytes = element.to_bytes_le();
            let expected_bytes: [u8; 4] = ByteConversion::to_bytes_le(&element).try_into().unwrap();
            assert_eq!(bytes, expected_bytes);
        }

        #[test]
        #[cfg(feature = "alloc")]
        fn byte_serialization_for_a_number_matches_with_byte_conversion_implementation_be() {
            let element = FE::from_hex("0123456701234567").unwrap();
            let bytes = element.to_bytes_be();
            let expected_bytes: [u8; 4] = ByteConversion::to_bytes_be(&element).try_into().unwrap();
            assert_eq!(bytes, expected_bytes);
        }

        #[test]
        fn byte_serialization_and_deserialization_works_le() {
            let element = FE::from_hex("0x7654321076543210").unwrap();
            let bytes = element.to_bytes_le();
            let from_bytes = FE::from_bytes_le(&bytes).unwrap();
            assert_eq!(element, from_bytes);
        }

        #[test]
        fn byte_serialization_and_deserialization_works_be() {
            let element = FE::from_hex("7654321076543210").unwrap();
            let bytes = element.to_bytes_be();
            let from_bytes = FE::from_bytes_be(&bytes).unwrap();
            assert_eq!(element, from_bytes);
        }
    }

    #[cfg(all(feature = "std", not(feature = "instruments")))]
    mod test_babybear_31_fft {
        use super::*;
        #[cfg(not(any(feature = "metal", feature = "cuda")))]
        use crate::fft::cpu::roots_of_unity::{
            get_powers_of_primitive_root, get_powers_of_primitive_root_coset,
        };
        use crate::field::element::FieldElement;
        #[cfg(not(any(feature = "metal", feature = "cuda")))]
        use crate::field::traits::{IsFFTField, RootsConfig};
        use crate::polynomial::Polynomial;
        use proptest::{collection, prelude::*, std_facade::Vec};

        #[cfg(not(any(feature = "metal", feature = "cuda")))]
        fn gen_fft_and_naive_evaluation<F: IsFFTField>(
            poly: Polynomial<FieldElement<F>>,
        ) -> (Vec<FieldElement<F>>, Vec<FieldElement<F>>) {
            let len = poly.coeff_len().next_power_of_two();
            let order = len.trailing_zeros();
            let twiddles =
                get_powers_of_primitive_root(order.into(), len, RootsConfig::Natural).unwrap();

            let fft_eval = Polynomial::evaluate_fft::<F>(&poly, 1, None).unwrap();
            let naive_eval = poly.evaluate_slice(&twiddles);

            (fft_eval, naive_eval)
        }

        #[cfg(not(any(feature = "metal", feature = "cuda")))]
        fn gen_fft_coset_and_naive_evaluation<F: IsFFTField>(
            poly: Polynomial<FieldElement<F>>,
            offset: FieldElement<F>,
            blowup_factor: usize,
        ) -> (Vec<FieldElement<F>>, Vec<FieldElement<F>>) {
            let len = poly.coeff_len().next_power_of_two();
            let order = (len * blowup_factor).trailing_zeros();
            let twiddles =
                get_powers_of_primitive_root_coset(order.into(), len * blowup_factor, &offset)
                    .unwrap();

            let fft_eval =
                Polynomial::evaluate_offset_fft::<F>(&poly, blowup_factor, None, &offset).unwrap();
            let naive_eval = poly.evaluate_slice(&twiddles);

            (fft_eval, naive_eval)
        }

        #[cfg(not(any(feature = "metal", feature = "cuda")))]
        fn gen_fft_and_naive_interpolate<F: IsFFTField>(
            fft_evals: &[FieldElement<F>],
        ) -> (Polynomial<FieldElement<F>>, Polynomial<FieldElement<F>>) {
            let order = fft_evals.len().trailing_zeros() as u64;
            let twiddles =
                get_powers_of_primitive_root(order, 1 << order, RootsConfig::Natural).unwrap();

            let naive_poly = Polynomial::interpolate(&twiddles, fft_evals).unwrap();
            let fft_poly = Polynomial::interpolate_fft::<F>(fft_evals).unwrap();

            (fft_poly, naive_poly)
        }

        #[cfg(not(any(feature = "metal", feature = "cuda")))]
        fn gen_fft_and_naive_coset_interpolate<F: IsFFTField>(
            fft_evals: &[FieldElement<F>],
            offset: &FieldElement<F>,
        ) -> (Polynomial<FieldElement<F>>, Polynomial<FieldElement<F>>) {
            let order = fft_evals.len().trailing_zeros() as u64;
            let twiddles = get_powers_of_primitive_root_coset(order, 1 << order, offset).unwrap();

            let naive_poly = Polynomial::interpolate(&twiddles, fft_evals).unwrap();
            let fft_poly = Polynomial::interpolate_offset_fft(fft_evals, offset).unwrap();

            (fft_poly, naive_poly)
        }

        #[cfg(not(any(feature = "metal", feature = "cuda")))]
        fn gen_fft_interpolate_and_evaluate<F: IsFFTField>(
            poly: Polynomial<FieldElement<F>>,
        ) -> (Polynomial<FieldElement<F>>, Polynomial<FieldElement<F>>) {
            let eval = Polynomial::evaluate_fft::<F>(&poly, 1, None).unwrap();
            let new_poly = Polynomial::interpolate_fft::<F>(&eval).unwrap();

            (poly, new_poly)
        }

        prop_compose! {
            fn powers_of_two(max_exp: u8)(exp in 1..max_exp) -> usize { 1 << exp }
            // max_exp cannot be multiple of the bits that represent a usize, generally 64 or 32.
            // also it can't exceed the test field's two-adicity.
        }
        prop_compose! {
            fn field_element()(num in any::<u64>().prop_filter("Avoid null coefficients", |x| x != &0)) -> FieldElement<Babybear31PrimeField> {
                FieldElement::<Babybear31PrimeField>::from(num)
            }
        }
        prop_compose! {
            fn offset()(num in any::<u64>(), factor in any::<u64>()) -> FieldElement<Babybear31PrimeField> { FieldElement::<Babybear31PrimeField>::from(num).pow(factor) }
        }
        prop_compose! {
            fn field_vec(max_exp: u8)(vec in collection::vec(field_element(), 0..1 << max_exp)) -> Vec<FieldElement<Babybear31PrimeField>> {
                vec
            }
        }
        prop_compose! {
            fn non_power_of_two_sized_field_vec(max_exp: u8)(vec in collection::vec(field_element(), 2..1<<max_exp).prop_filter("Avoid polynomials of size power of two", |vec| !vec.len().is_power_of_two())) -> Vec<FieldElement<Babybear31PrimeField>> {
                vec
            }
        }
        prop_compose! {
            fn poly(max_exp: u8)(coeffs in field_vec(max_exp)) -> Polynomial<FieldElement<Babybear31PrimeField>> {
                Polynomial::new(&coeffs)
            }
        }
        prop_compose! {
            fn poly_with_non_power_of_two_coeffs(max_exp: u8)(coeffs in non_power_of_two_sized_field_vec(max_exp)) -> Polynomial<FieldElement<Babybear31PrimeField>> {
                Polynomial::new(&coeffs)
            }
        }

        proptest! {
            // Property-based test that ensures FFT eval. gives same result as a naive polynomial evaluation.
            #[test]
            #[cfg(not(any(feature = "metal",feature = "cuda")))]
            fn test_fft_matches_naive_evaluation(poly in poly(8)) {
                let (fft_eval, naive_eval) = gen_fft_and_naive_evaluation(poly);
                prop_assert_eq!(fft_eval, naive_eval);
            }

            // Property-based test that ensures FFT eval. with coset gives same result as a naive polynomial evaluation.
            #[test]
            #[cfg(not(any(feature = "metal",feature = "cuda")))]
            fn test_fft_coset_matches_naive_evaluation(poly in poly(4), offset in offset(), blowup_factor in powers_of_two(4)) {
                let (fft_eval, naive_eval) = gen_fft_coset_and_naive_evaluation(poly, offset, blowup_factor);
                prop_assert_eq!(fft_eval, naive_eval);
            }

            // #[cfg(not(any(feature = "metal"),not(feature = "cuda")))]
            // Property-based test that ensures FFT interpolation is the same as naive..
            #[test]
            #[cfg(not(any(feature = "metal",feature = "cuda")))]
            fn test_fft_interpolate_matches_naive(fft_evals in field_vec(4)
                                                           .prop_filter("Avoid polynomials of size not power of two",
                                                                        |evals| evals.len().is_power_of_two())) {
                let (fft_poly, naive_poly) = gen_fft_and_naive_interpolate(&fft_evals);
                prop_assert_eq!(fft_poly, naive_poly);
            }

            // Property-based test that ensures FFT interpolation with an offset is the same as naive.
            #[test]
            #[cfg(not(any(feature = "metal",feature = "cuda")))]
            fn test_fft_interpolate_coset_matches_naive(offset in offset(), fft_evals in field_vec(4)
                                                           .prop_filter("Avoid polynomials of size not power of two",
                                                                        |evals| evals.len().is_power_of_two())) {
                let (fft_poly, naive_poly) = gen_fft_and_naive_coset_interpolate(&fft_evals, &offset);
                prop_assert_eq!(fft_poly, naive_poly);
            }

            // Property-based test that ensures interpolation is the inverse operation of evaluation.
            #[test]
            #[cfg(not(any(feature = "metal",feature = "cuda")))]
            fn test_fft_interpolate_is_inverse_of_evaluate(
                poly in poly(4).prop_filter("Avoid non pows of two", |poly| poly.coeff_len().is_power_of_two())) {
                let (poly, new_poly) = gen_fft_interpolate_and_evaluate(poly);
                prop_assert_eq!(poly, new_poly);
            }
        }
    }
}
