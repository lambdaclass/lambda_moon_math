extern crate alloc;
use super::point::CirclePoint;
#[cfg(feature = "alloc")]
use super::{
    cfft::{cfft, icfft, order_cfft_result_naive, order_icfft_input_naive},
    cosets::Coset,
    twiddles::{get_twiddles, TwiddlesConfig},
};
use crate::{
    fft::cpu::bit_reversing::in_place_bit_reverse_permute,
    field::{element::FieldElement, fields::mersenne31::field::Mersenne31Field},
};
use alloc::vec::Vec;

/// Given the 2^n coefficients of a two-variables polynomial of degree 2^n - 1 in the basis {1, y, x, xy, 2xˆ2 -1, 2xˆ2y-y, 2xˆ3-x, 2xˆ3y-xy,...}
/// returns the evaluation of the polynomial on the points of the standard coset of size 2^n.
/// Note that coeff has to be a vector with length a power of two 2^n.
#[cfg(feature = "alloc")]
pub fn evaluate_cfft(
    coeff: Vec<FieldElement<Mersenne31Field>>,
) -> Vec<FieldElement<Mersenne31Field>> {
    let mut coeff = coeff;

    // We get the twiddles for the Evaluation.
    let domain_log_2_size: u32 = coeff.len().trailing_zeros();
    let coset = Coset::new_standard(domain_log_2_size);
    let config = TwiddlesConfig::Evaluation;
    let twiddles = get_twiddles(coset, config);

    // For our algorithm to work, we must give as input the coefficients in bit reverse order.
    in_place_bit_reverse_permute::<FieldElement<Mersenne31Field>>(&mut coeff);
    cfft(&mut coeff, twiddles);

    // The cfft returns the evaluations in a certain order, so we permute them to get the natural order.
    order_cfft_result_naive(&coeff)
}

/// Interpolates the 2^n evaluations of a two-variables polynomial of degree 2^n - 1 on the points of the standard coset of size 2^n.
/// As a result we obtain the coefficients of the polynomial in the basis: {1, y, x, xy, 2xˆ2 -1, 2xˆ2y-y, 2xˆ3-x, 2xˆ3y-xy,...}
/// Note that eval has to be a vector of length a power of two 2^n.
/// If the vector of evaluations is empty, it returns an empty vector.
#[cfg(feature = "alloc")]
pub fn interpolate_cfft(
    eval: Vec<FieldElement<Mersenne31Field>>,
) -> Vec<FieldElement<Mersenne31Field>> {
    let mut eval = eval;

    if eval.is_empty() {
        let poly: Vec<FieldElement<Mersenne31Field>> = Vec::new();
        return poly;
    }

    // We get the twiddles for the interpolation.
    let domain_log_2_size: u32 = eval.len().trailing_zeros();
    let coset = Coset::new_standard(domain_log_2_size);
    let config = TwiddlesConfig::Interpolation;
    let twiddles = get_twiddles(coset, config);

    // For our algorithm to work, we must give as input the evaluations ordered in a certain way.
    let mut eval_ordered = order_icfft_input_naive(&mut eval);
    icfft(&mut eval_ordered, twiddles);

    // The icfft returns the polynomial coefficients in bit reverse order. So we premute it to get the natural order.
    in_place_bit_reverse_permute::<FieldElement<Mersenne31Field>>(&mut eval_ordered);

    // The icfft returns all the coefficients multiplied by 2^n, the length of the evaluations.
    // So we multiply every element that outputs the icfft by the inverse of 2^n to get the actual coefficients.
    // Note that this `unwrap` will never panic because eval.len() != 0.
    let factor = (FieldElement::<Mersenne31Field>::from(eval.len() as u64))
        .inv()
        .unwrap();
    eval_ordered.iter().map(|coef| coef * factor).collect()
}

/// Note: This implementation uses a straightforward approach and is intended for testing purposes only.
pub fn evaluate_point(
    coef: &Vec<FieldElement<Mersenne31Field>>,
    point: CirclePoint<Mersenne31Field>,
) -> FieldElement<Mersenne31Field> {
    let order = coef.len();
    assert!(
        order.is_power_of_two(),
        "Coefficient length must be a power of 2"
    );

    let v_len = order.trailing_zeros() as usize;

    let mut a = point.x;
    let mut v = Vec::with_capacity(v_len);
    v.push(FieldElement::one());
    v.push(point.x);
    for _ in 2..v_len {
        a = a.square().double() - FieldElement::one();
        v.push(a);
    }

    let mut result = FieldElement::zero();

    for i in 0..order {
        let mut term = coef[i];

        if i % 2 == 1 {
            term = term * point.y;
        }

        let mut idx = i / 2;
        let mut pos = 0;
        while idx > 0 {
            if idx % 2 == 1 {
                term = term * v[pos + 1];
            }
            idx /= 2;
            pos += 1;
        }

        result = result + term;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circle::cosets::Coset;
    type FE = FieldElement<Mersenne31Field>;
    use alloc::vec;

    /// Naive evaluation of a polynomial of degree 3.
    fn evaluate_poly_4(coef: &[FE; 4], x: FE, y: FE) -> FE {
        coef[0] + coef[1] * y + coef[2] * x + coef[3] * x * y
    }

    /// Naive evaluation of a polynomial of degree 7.
    fn evaluate_poly_8(coef: &[FE; 8], x: FE, y: FE) -> FE {
        coef[0]
            + coef[1] * y
            + coef[2] * x
            + coef[3] * x * y
            + coef[4] * (x.square().double() - FE::one())
            + coef[5] * (x.square().double() - FE::one()) * y
            + coef[6] * ((x.square() * x).double() - x)
            + coef[7] * ((x.square() * x).double() - x) * y
    }

    /// Naive evaluation of a polynomial of degree 15.
    fn evaluate_poly_16(coef: &[FE; 16], x: FE, y: FE) -> FE {
        let mut a = x;
        let mut v = Vec::new();
        v.push(FE::one());
        v.push(x);
        for _ in 2..4 {
            a = a.square().double() - FE::one();
            v.push(a);
        }

        coef[0] * v[0]
            + coef[1] * y * v[0]
            + coef[2] * v[1]
            + coef[3] * y * v[1]
            + coef[4] * v[2]
            + coef[5] * y * v[2]
            + coef[6] * v[1] * v[2]
            + coef[7] * y * v[1] * v[2]
            + coef[8] * v[3]
            + coef[9] * y * v[3]
            + coef[10] * v[1] * v[3]
            + coef[11] * y * v[1] * v[3]
            + coef[12] * v[2] * v[3]
            + coef[13] * y * v[2] * v[3]
            + coef[14] * v[1] * v[2] * v[3]
            + coef[15] * y * v[1] * v[2] * v[3]
    }

    #[test]
    /// cfft evaluation equals naive evaluation.
    fn cfft_evaluation_4_points() {
        // We define the coefficients of a polynomial of degree 3.
        let input = [FE::from(1), FE::from(2), FE::from(3), FE::from(4)];

        // We create the coset points and evaluate the polynomial with the naive function.
        let coset = Coset::new_standard(2);
        let points = Coset::get_coset_points(&coset);
        let mut expected_result: Vec<FE> = Vec::new();
        for point in points {
            let point_eval = evaluate_poly_4(&input, point.x, point.y);
            expected_result.push(point_eval);
        }

        let input_vec = input.to_vec();
        // We evaluate the polynomial using now the cfft.
        let result = evaluate_cfft(input_vec);
        let slice_result: &[FE] = &result;

        assert_eq!(slice_result, expected_result);
    }

    #[test]
    /// cfft evaluation equals naive evaluation.
    fn cfft_evaluation_8_points() {
        // We define the coefficients of a polynomial of degree 7.
        let input = [
            FE::from(1),
            FE::from(2),
            FE::from(3),
            FE::from(4),
            FE::from(5),
            FE::from(6),
            FE::from(7),
            FE::from(8),
        ];

        // We create the coset points and evaluate them without the fft.
        let coset = Coset::new_standard(3);
        let points = Coset::get_coset_points(&coset);
        let mut expected_result: Vec<FE> = Vec::new();
        for point in points {
            let point_eval = evaluate_poly_8(&input, point.x, point.y);
            expected_result.push(point_eval);
        }

        // We evaluate the polynomial using now the cfft.
        let result = evaluate_cfft(input.to_vec());
        let slice_result: &[FE] = &result;

        assert_eq!(slice_result, expected_result);
    }

    #[test]
    /// cfft evaluation equals naive evaluation.
    fn cfft_evaluation_16_points() {
        // We define the coefficients of a polynomial of degree 15.
        let input = [
            FE::from(1),
            FE::from(2),
            FE::from(3),
            FE::from(4),
            FE::from(5),
            FE::from(6),
            FE::from(7),
            FE::from(8),
            FE::from(9),
            FE::from(10),
            FE::from(11),
            FE::from(12),
            FE::from(13),
            FE::from(14),
            FE::from(15),
            FE::from(16),
        ];

        // We create the coset points and evaluate them without the fft.
        let coset = Coset::new_standard(4);
        let points = Coset::get_coset_points(&coset);
        let mut expected_result: Vec<FE> = Vec::new();
        for point in points {
            let point_eval = evaluate_poly_16(&input, point.x, point.y);
            expected_result.push(point_eval);
        }

        // We evaluate the polynomial using now the cfft.
        let result = evaluate_cfft(input.to_vec());
        let slice_result: &[FE] = &result;

        assert_eq!(slice_result, expected_result);
    }

    #[test]
    fn evaluate_and_interpolate_8_points_is_identity() {
        // We define the 8 coefficients of a polynomial of degree 7.
        let coeff = vec![
            FE::from(1),
            FE::from(2),
            FE::from(3),
            FE::from(4),
            FE::from(5),
            FE::from(6),
            FE::from(7),
            FE::from(8),
        ];
        let evals = evaluate_cfft(coeff.clone());
        let new_coeff = interpolate_cfft(evals);

        assert_eq!(coeff, new_coeff);
    }

    #[test]
    fn evaluate_and_interpolate_8_other_points() {
        let coeff = vec![
            FE::from(2147483650),
            FE::from(147483647),
            FE::from(2147483700),
            FE::from(2147483647),
            FE::from(3147483647),
            FE::from(4147483647),
            FE::from(2147483640),
            FE::from(5147483647),
        ];
        let evals = evaluate_cfft(coeff.clone());
        let new_coeff = interpolate_cfft(evals);

        assert_eq!(coeff, new_coeff);
    }

    #[test]
    fn evaluate_and_interpolate_32_points() {
        // We define 32 coefficients of a polynomial of degree 31.
        let coeff = vec![
            FE::from(1),
            FE::from(2),
            FE::from(3),
            FE::from(4),
            FE::from(5),
            FE::from(6),
            FE::from(7),
            FE::from(8),
            FE::from(9),
            FE::from(10),
            FE::from(11),
            FE::from(12),
            FE::from(13),
            FE::from(14),
            FE::from(15),
            FE::from(16),
            FE::from(17),
            FE::from(18),
            FE::from(19),
            FE::from(20),
            FE::from(21),
            FE::from(22),
            FE::from(23),
            FE::from(24),
            FE::from(25),
            FE::from(26),
            FE::from(27),
            FE::from(28),
            FE::from(29),
            FE::from(30),
            FE::from(31),
            FE::from(32),
        ];
        let evals = evaluate_cfft(coeff.clone());

        let coset = Coset::new_standard(5);
        let coset_points = Coset::get_coset_points(&coset);

        assert_eq!(evals[4], evaluate_point(&coeff, coset_points[4].clone()));

        let new_coeff = interpolate_cfft(evals);

        assert_eq!(coeff, new_coeff);
    }
}
