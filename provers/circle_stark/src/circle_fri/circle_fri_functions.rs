use lambdaworks_math::circle::cosets::Coset;
use lambdaworks_math::circle::twiddles::{get_twiddles, TwiddlesConfig};
use lambdaworks_math::field::fields::mersenne31::field::Mersenne31Field;
use lambdaworks_math::field::{
    element::FieldElement,
    traits::{IsFFTField, IsField, IsSubFieldOf},
};

type F = Mersenne31Field;
type FE = FieldElement<F>;

/// First folding step in the Circle FRI protocol, converting a bivariate polynomial over the circle
/// into a univariate polynomial.
pub fn fold_circle_to_univariate(evaluations: &[FE], beta: &FE, domain: &Coset) -> Vec<FE> {
    let n = evaluations.len() / 2;
    let twiddles = get_twiddles(domain.clone(), TwiddlesConfig::Evaluation);

    // Use the corresponding twiddles for folding
    evaluations
        .iter()
        .take(n)
        .zip(evaluations.iter().skip(n))
        .enumerate()
        .map(|(i, (f_p, f_neg_p))| {
            let y_twiddle = &twiddles[0][i];
            let sum = f_p + f_neg_p;
            let diff = (f_p - f_neg_p) * y_twiddle;
            // Ensure the constant is of type FieldElement<F>
            let inv_two = FieldElement::<F>::from(2u64).inv().unwrap();
            (sum + beta * diff) * inv_two
        })
        .collect()
}

/// Subsequent folding steps in the Circle FRI protocol using univariate polynomials.
pub fn fold_univariate(
    evaluations: &[FE],
    beta: &FieldElement<F>,
    domain: &Coset,
    layer_index: usize,
) -> Vec<FE> {
    let n = evaluations.len() / 2;
    let twiddles = get_twiddles(domain.clone(), TwiddlesConfig::Evaluation);
    let twiddle_layer = &twiddles[layer_index];

    evaluations
        .iter()
        .take(n)
        .zip(evaluations.iter().skip(n))
        .enumerate()
        .map(|(i, (f_x, f_neg_x))| {
            let twiddle = &twiddle_layer[i];
            let sum = f_x + f_neg_x;
            let diff = (f_x - f_neg_x) * twiddle;
            // Ensure the constant is of type FieldElement<F>
            let inv_two = FieldElement::<F>::from(2u64).inv().unwrap();
            (sum + beta * diff) * inv_two
        })
        .collect()
}

/*
#[cfg(test)]
mod tests {
    use super::{fold_circle_to_univariate, fold_univariate};
    use lambdaworks_math::{
        circle::cosets::Coset,
        fft::cpu::roots_of_unity::get_powers_of_primitive_root_coset,
        field::{
            element::FieldElement,
            fields::mersenne31::field::Mersenne31Field,
            traits::{IsFFTField, IsField},
        },
    };

    type F = Mersenne31Field;
    type FE = FieldElement<F>;

    #[test]
    fn test_fold_circle_to_univariate() {
        // Define the field
        let two = FE::from(2u64);
        let three = FE::from(3u64);
        let five = FE::from(5u64);
        let seven = FE::from(7u64);
        let eleven = FE::from(11u64);
        let thirteen = FE::from(13u64);

        // Create sample evaluations (size must be a power of two)
        let evaluations = vec![two, three, five, seven, eleven, thirteen, two, three];

        // Define beta challenge
        let beta = FE::from(4u64);

        // Define initial domain (coset)
        let coset_shift = FE::from(1u64); // For simplicity, using 1 as shift
        let log_2_size = 3; // Since evaluations.len() = 8, log_2(8) = 3
        let domain = Coset {
            shift: coset_shift.clone(),
            log_2_size,
        };

        // Perform the folding
        let folded_evaluations = fold_circle_to_univariate(&evaluations, &beta, &domain);

        // Expected output (you need to compute this based on your implementation)
        // For illustration purposes, let's assume some expected values
        let expected_folded_evaluations = vec![
            FE::from(5u64),  // Placeholder values
            FE::from(10u64), // You should compute the actual expected values
            FE::from(15u64), // based on the mathematical operations in your function
            FE::from(20u64),
        ];

        // Verify the output
        assert_eq!(folded_evaluations, expected_folded_evaluations);
    }
}

*/
