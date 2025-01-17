use super::Polynomial;
use lambdaworks_math::{
    field::{element::FieldElement, traits::IsField},
    polynomial,
};

pub fn fold_polynomial<F>(
    poly: &Polynomial<FieldElement<F>>,
    beta: &FieldElement<F>,
) -> Polynomial<FieldElement<F>>
where
    F: IsField,
{
    let coef = poly.coefficients();
    let even_coef: Vec<FieldElement<F>> = coef.iter().step_by(2).cloned().collect();

    // odd coeficients of poly are multiplied by beta
    let odd_coef_mul_beta: Vec<FieldElement<F>> = coef
        .iter()
        .skip(1)
        .step_by(2)
        .map(|v| (v.clone()) * beta)
        .collect();

    let (even_poly, odd_poly) = polynomial::pad_with_zero_coefficients(
        &Polynomial::new(&even_coef),
        &Polynomial::new(&odd_coef_mul_beta),
    );
    even_poly + odd_poly
}

#[cfg(test)]
mod tests {
    use super::fold_polynomial;
    use lambdaworks_math::field::element::FieldElement;
    use lambdaworks_math::field::fields::u64_prime_field::U64PrimeField;
    const MODULUS: u64 = 293;
    type FE = FieldElement<U64PrimeField<MODULUS>>;
    use lambdaworks_math::polynomial::Polynomial;

    #[test]
    fn test_fold() {
        let p0 = Polynomial::new(&[
            FE::new(3),
            FE::new(1),
            FE::new(2),
            FE::new(7),
            FE::new(3),
            FE::new(5),
        ]);
        let beta = FE::new(4);
        let p1 = fold_polynomial(&p0, &beta);
        assert_eq!(
            p1,
            Polynomial::new(&[FE::new(7), FE::new(30), FE::new(23),])
        );

        let gamma = FE::new(3);
        let p2 = fold_polynomial(&p1, &gamma);
        assert_eq!(p2, Polynomial::new(&[FE::new(97), FE::new(23),]));

        let delta = FE::new(2);
        let p3 = fold_polynomial(&p2, &delta);
        assert_eq!(p3, Polynomial::new(&[FE::new(143)]));
        assert_eq!(p3.degree(), 0);
    }
}
