use lambdaworks_math::circle::cosets::Coset;
use lambdaworks_math::circle::twiddles::{get_twiddles, TwiddlesConfig};
use lambdaworks_math::field::fields::mersenne31::field::Mersenne31Field;
use lambdaworks_math::field::{
    element::FieldElement,
    traits::{IsField, IsSubFieldOf},
};

/// First folding step in the Circle FRI protocol, converting a bivariate polynomial over the circle
/// into a univariate polynomial.
pub fn fold_circle_to_univariate<F>(
    evaluations: &[FieldElement<F>],
    beta: &FieldElement<F>,
    domain: &Coset,
) -> Vec<FieldElement<F>>
where
    F: IsField + IsSubFieldOf<Mersenne31Field>,
    Mersenne31Field: IsSubFieldOf<F>,
{
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
pub fn fold_univariate<F>(
    evaluations: &[FieldElement<F>],
    beta: &FieldElement<F>,
    domain: &Coset,
    layer_index: usize,
) -> Vec<FieldElement<F>>
where
    F: IsField + IsSubFieldOf<Mersenne31Field>,
    Mersenne31Field: IsSubFieldOf<F>,
{
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
