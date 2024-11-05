use itertools::izip;
use lambdaworks_math::field::fields::mersenne31::extensions::Degree4ExtensionField;
use lambdaworks_math::{
    circle::{
        domain::CircleDomain,
        point::CirclePoint,
        polynomial::{evaluate_cfft, evaluate_point},
    },
    field::{
        element::FieldElement,
        //extensions::binomial::BinomialExtension,
        fields::mersenne31::field::Mersenne31Field,
        traits::{IsField, IsSubFieldOf},
    },
};
type F = Mersenne31Field;
type EF = Degree4ExtensionField;

/// Computes numerator and denominator of the "vanishing part" of the DEEP quotient
/// Section 6, Remark 21 of Circle Starks
/// Re(1/v_gamma) + alpha^L Im(1/v_gamma)
pub fn deep_quotient_vanishing_part(
    x: CirclePoint<F>,
    zeta: CirclePoint<EF>,
    alpha_pow_width: FieldElement<EF>,
) -> (FieldElement<EF>, FieldElement<EF>) {
    let (re_v_zeta, im_v_zeta) = v_p(x, zeta);
    let re_v_zeta_clone = re_v_zeta.clone();
    (
        re_v_zeta - alpha_pow_width * &im_v_zeta,
        re_v_zeta_clone.square() + im_v_zeta.clone().square(),
    )
}
// should move to point.rs
pub fn v_p(point: CirclePoint<F>, at: CirclePoint<EF>) -> (FieldElement<EF>, FieldElement<EF>) {
    // Caculate at infinity?
    let diff = CirclePoint {
        x: -at.x + point.x.to_extension(),
        y: -at.y + point.y.to_extension(),
    };

    (FieldElement::<EF>::one() - diff.x, -diff.y)
}

pub fn deep_quotient_reduce_row(
    alpha: FieldElement<EF>,
    x: CirclePoint<F>,
    zeta: CirclePoint<EF>,
    ps_at_x: &[FieldElement<F>],
    ps_at_zeta: &[FieldElement<EF>],
) -> FieldElement<EF> {
    let (vp_num, vp_denom) = deep_quotient_vanishing_part(x, zeta, alpha.pow(ps_at_x.len() as u64));

    // Generate an iterator of powers of `alpha` up to `ps_at_x.len()`
    let powers_iter = powers(alpha, ps_at_x.len());

    // Compute terms for the dot product
    let terms =
        izip!(ps_at_x, ps_at_zeta).map(|(p_at_x, p_at_zeta)| -p_at_zeta + p_at_x.to_extension());

    // Calculate the dot product using the iterator of powers
    let powers_vec: Vec<_> = powers_iter.collect();
    let terms_vec: Vec<_> = terms.collect();
    let dot_product = dot_product(&powers_vec, &terms_vec);

    (vp_num / vp_denom) * dot_product
}

pub fn deep_quotient_reduce(
    coeff: Vec<FieldElement<F>>, // Coeficientes de entrada
    alpha: FieldElement<EF>,
    zeta: CirclePoint<EF>,
    ps_at_zeta: &[FieldElement<EF>],
) -> Vec<FieldElement<EF>> {
    // Realiza la FFT sobre los coeficientes para obtener evaluaciones en el coset.
    let mut evaluations = evaluate_cfft(coeff);

    // Calcula el factor `alpha_pow_width` para normalizar las evaluaciones.
    let alpha_pow_width = alpha.pow(evaluations.len() as u64);

    // Calcula la parte que se anula en cada evaluación
    let (vp_nums, vp_denoms): (Vec<_>, Vec<_>) = evaluations
        .iter()
        .map(|eval| deep_quotient_vanishing_part(*eval, zeta, alpha_pow_width.clone()))
        .unzip();

    let vp_denom_invs = batch_multiplicative_inverse(&vp_denoms);

    // Calcula la reducción en `ps_at_zeta`
    let alpha_powers: Vec<_> = powers(alpha.clone(), ps_at_zeta.len()).collect();
    let ps_at_zeta_cloned: Vec<_> = ps_at_zeta.iter().cloned().collect();
    let alpha_reduced_ps_at_zeta = dot_product(&alpha_powers, &ps_at_zeta_cloned);

    // Reduce cada evaluación de acuerdo al deep quotient y almacena los resultados
    evaluations
        .iter_mut()
        .enumerate()
        .map(|(i, eval)| {
            let alpha_powers_eval: Vec<_> = powers(alpha.clone(), evaluations.len()).collect();
            let eval_extension = eval.to_extension();
            let reduced_eval = dot_product(&alpha_powers_eval, &eval_extension);
            (vp_nums[i] / vp_denom_invs[i]) * (reduced_eval - alpha_reduced_ps_at_zeta.clone())
        })
        .collect()
}

/// Extrae el múltiplo del polinomio de anulación del dominio original
/// Ver Sección 4.3, Lemma 6: < v_n, f > = 0 para cualquier f en espacio FFT
pub fn extract_lambda(lde: &mut [FieldElement<EF>], log_blowup: usize) -> FieldElement<EF> {
    let log_lde_size = lde.len().trailing_zeros() as usize;

    // v_n es constante en cosets del mismo tamaño que el dominio original
    let v_d_init = CircleDomain::standard(log_lde_size)
        .points()
        .take(1 << log_blowup)
        .map(|p| p.v_n(log_lde_size - log_blowup))
        .collect::<Vec<_>>();

    // Los valores únicos se repiten sobre el resto del dominio como
    // 0 1 2 .. n-1 n n n-1 .. 1 0 0 1 ..
    let v_d: Vec<_> = v_d_init
        .iter()
        .chain(v_d_init.iter().rev())
        .cycle()
        .take(lde.len())
        .cloned()
        .collect();

    // < v_d, v_d >
    let v_d_2 = FieldElement::from(2).pow(log_lde_size as u64 - 1);
    let v_d_2_inv = v_d_2.inv().unwrap();

    // Calcula lambda
    let lambda = dot_product(&lde, &v_d) * v_d_2_inv;

    // Actualiza lde
    for (y, v_x) in lde.iter_mut().zip(v_d.iter()) {
        *y = y.clone() - lambda.clone() * v_x;
    }

    lambda
}

// aux functions

fn dot_product(a: &[FieldElement<EF>], b: &[FieldElement<EF>]) -> FieldElement<EF> {
    let mut sum = FieldElement::zero();
    for i in 0..a.len() {
        sum = sum + a[i].clone() * b[i].clone();
    }
    sum
}
/// Returns an iterator that generates successive powers of a `FieldElement`.
/// Starts at 1 and multiplies repeatedly by `base`.
pub fn powers<F: IsField>(base: FieldElement<F>, count: usize) -> Powers<F> {
    Powers {
        base,
        current: FieldElement::one(),
        count,
        index: 0,
    }
}

/// Structure to iterate over the powers of a `FieldElement`.
pub struct Powers<F: IsField> {
    base: FieldElement<F>,
    current: FieldElement<F>,
    count: usize,
    index: usize,
}

impl<F: IsField> Iterator for Powers<F> {
    type Item = FieldElement<F>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.count {
            let result = self.current.clone();
            self.current = self.current.clone() * self.base.clone();
            self.index += 1;
            Some(result)
        } else {
            None
        }
    }
}

fn batch_multiplicative_inverse<F: IsField>(elements: &[FieldElement<F>]) -> Vec<FieldElement<F>> {
    if elements.is_empty() {
        return vec![];
    }

    let mut acc = FieldElement::one();
    let mut products: Vec<_> = elements
        .iter()
        .map(|x| {
            acc = acc.clone() * x.clone();
            acc.clone()
        })
        .collect();

    let mut inv = acc.inv().unwrap();
    for i in (0..elements.len()).rev() {
        let new_inv = inv.clone() * elements[i].clone();
        products[i] = inv;
        inv = new_inv;
    }
    products
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{thread_rng, Rng};

    #[test]
    fn test_deep_quotient_reduce() {
        let domain = CircleDomain::standard(5);
        let values = vec![
            vec![FieldElement::from(1); 8],
            vec![FieldElement::from(2); 8],
        ];
        let evals = CircleEvaluations::new(domain, values);

        let alpha: FieldElement<EF> = FieldElement::random(&mut thread_rng());
        let zeta = CirclePoint::random(&mut thread_rng());

        let ps_at_zeta = evals.evaluate_at_point(&zeta);
        let reduced = deep_quotient_reduce(&evals, alpha, zeta, &ps_at_zeta);

        assert_eq!(reduced.len(), evals.values().len());
    }

    #[test]
    fn test_extract_lambda() {
        let log_n = 5;
        let log_blowup = 1;

        let mut lde = vec![FieldElement::zero(); 1 << (log_n + log_blowup)];
        // Fill lde with some test values
        for i in 0..lde.len() {
            lde[i] = FieldElement::from(i as u64);
        }

        let lambda = extract_lambda(&mut lde, log_blowup);
        assert!(!lambda.is_zero());
    }
}
