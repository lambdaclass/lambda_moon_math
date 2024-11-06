use itertools::izip;
use lambdaworks_math::circle::cosets::Coset;
use lambdaworks_math::field::fields::mersenne31::extensions::Degree4ExtensionField;
use lambdaworks_math::{
    circle::{
        //  domain::CircleDomain,
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
///
/// Returns the numerator and denominator as a tuple
/// (Re(vγ) - μL · Im(vγ) , Re(vγ)² + Im(vγ)²)
// There is another part of the vanishing part that is not calculated here
// (ḡ - v̄γ) is calulated in the deep_quotient_reduce_row function
pub fn deep_quotient_vanishing_part(
    x: CirclePoint<F>,
    zeta: CirclePoint<EF>,
    alpha_pow_width: FieldElement<EF>,
) -> (FieldElement<EF>, FieldElement<EF>) {
    let (re_v_zeta, im_v_zeta) = v_p(x, zeta);
    let numerator = &re_v_zeta - alpha_pow_width * &im_v_zeta;
    let denominator = re_v_zeta.square() + im_v_zeta.square();
    (numerator, denominator)
}
// should move to point.rs
pub fn v_p(point: CirclePoint<F>, at: CirclePoint<EF>) -> (FieldElement<EF>, FieldElement<EF>) {
    // vP(x,y) = v(Px·x + Py·y, -Py·x + Px·y)
    //= 1 - ((Px·x + Py·y) + i·(-Py·x + Px·y))
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
    ps_at_x: &[FieldElement<F>],     // [p₁(x), p₂(x), p₃(x).....]
    ps_at_zeta: &[FieldElement<EF>], // [p₁(ζ), p₂(ζ), p₃(ζ).....]
) -> FieldElement<EF> {
    let (vp_numerator, vp_denomnator) =
        deep_quotient_vanishing_part(x, zeta, alpha.pow(ps_at_x.len() as u64));

    // Generate an iterator of powers of `alpha` up to `ps_at_x.len()`
    // Given alpha = 2 and ps_at_x.len() = 3, the powers generated are [2^0, 2^1, 2^2] = [1, 2, 4].

    let powers_iter = powers(alpha, ps_at_x.len());
    // is this the same as alpha.powers in plonky3

    // Compute terms for the dot product
    let terms =
        izip!(ps_at_x, ps_at_zeta).map(|(p_at_x, p_at_zeta)| -p_at_zeta + p_at_x.to_extension());

    // Calculate the dot product using the iterator of powers
    let powers_vec: Vec<_> = powers_iter.collect();
    let terms_vec: Vec<_> = terms.collect();
    let dot_product = dot_product(&powers_vec, &terms_vec);

    (vp_numerator / vp_denomnator) * dot_product

    // Vanishing Part               Differences with powers of alpha
    // ┌────────────────────────────┐  ┌─────────────────────────────────────┐
    //   Re(vγ) - μL · Im(vγ)         Σ (pᵢ(x) - pᵢ(ζ))·αⁱ
    //   ───────────────────── ·
    //   Re(vγ)² + Im(vγ)²
}

pub fn deep_quotient_reduce(
    coeff: Vec<FieldElement<F>>, // Coeficientes de entrada
    alpha: FieldElement<EF>,
    zeta: CirclePoint<EF>,
    ps_at_zeta: &[FieldElement<EF>],
    log_n: usize,
) -> Vec<FieldElement<EF>> {
    // Realiza la FFT sobre los coeficientes para obtener evaluaciones en el coset.
    let mut evaluations = evaluate_cfft(coeff);
    // plonky3's has the evaluations so i need to evaluate and check if I
    // need to permform a permutation to get a correct evaluation
    // TO DO: ask Nicole and Colo about the permutation

    // Calcula el factor `alpha_pow_width` para normalizar las evaluaciones.
    let alpha_pow_width = alpha.pow(evaluations.len() as u64);

    // Calcula la parte que se anula en cada evaluación
    // Generate or retrieve your domain points
    let coset = Coset::new_standard(log_n as u32);
    let points = Coset::get_coset_points(&coset);

    // Apply the permutation to reorder the points if necessary
    let permuted_points: Vec<_> = (0..points.len())
        .map(|i| points[cfft_permute_index(i, log_n)].clone())
        .collect();

    // Compute the vanishing part for each permuted point
    let (vp_nums, vp_denoms): (Vec<_>, Vec<_>) = permuted_points
        .iter()
        .map(|point| {
            deep_quotient_vanishing_part(point.clone(), zeta.clone(), alpha_pow_width.clone())
        })
        .unzip();
    // This will invert the denominators using batch inversion which is optimized for this case
    let vp_denom_invs = batch_multiplicative_inverse(&vp_denoms);

    // need to make this more efficient and clear
    let alpha_powers: Vec<_> = powers(alpha.clone(), ps_at_zeta.len()).collect(); //
    let ps_at_zeta_cloned: Vec<_> = ps_at_zeta.iter().cloned().collect();
    let alpha_reduced_ps_at_zeta = dot_product(&alpha_powers, &ps_at_zeta_cloned);
    // need to make this more efficient and clear
    let eval_len = evaluations.len();
    let alpha_powers_eval: Vec<_> = powers(alpha.clone(), eval_len).collect();
    evaluations
        .iter_mut()
        .enumerate()
        .map(|(i, eval)| {
            let eval_extension = eval.to_extension();
            let eval_extension_converted: FieldElement<Degree4ExtensionField> =
                eval_extension.into();
            let reduced_eval = dot_product(&alpha_powers_eval, &[eval_extension_converted]);
            (vp_nums[i].clone() / vp_denom_invs[i].clone())
                * (reduced_eval - alpha_reduced_ps_at_zeta.clone())
        })
        .collect()
}

/// Extrae el múltiplo del polinomio de anulación del dominio original
/// Ver Sección 4.3, Lemma 6: < v_n, f > = 0 para cualquier f en espacio FFT
pub fn extract_lambda(lde: &mut [FieldElement<EF>], log_blowup: usize) -> FieldElement<EF> {
    let log_lde_size = lde.len().trailing_zeros() as usize;

    // v_n es constante en cosets del mismo tamaño que el dominio original
    let coset = Coset::new_standard(log_lde_size as u32);
    let v_d_init = Coset::get_coset_points(&coset)
        .into_iter()
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

    // Convert v_d to extension field elements
    let v_d_ext: Vec<FieldElement<EF>> = v_d.iter().map(|x| x.to_extension()).collect();

    // Calcula lambda
    let lambda = dot_product(&lde, &v_d_ext) * v_d_2_inv;

    // Actualiza lde
    for (y, v_x) in lde.iter_mut().zip(v_d_ext.iter()) {
        *y = y.clone() - lambda.clone() * v_x.clone();
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
pub fn cfft_permute_index(index: usize, log_n: usize) -> usize {
    let mut result = 0;
    for i in 0..log_n {
        result |= ((index >> i) & 1) << (log_n - 1 - i);
    }
    result
}
/// Evaluates the vanishing polynomial for the given CirclePoint with respect to `log_n`.
/// This function follows the iterative structure described and performs operations on the `x` coordinate.
pub fn v_n(mut point: CirclePoint<F>, log_n: usize) -> FieldElement<F> {
    // Applies the polynomial transformation iteratively based on log_n.
    for _ in 0..(log_n - 1) {
        point.x = point.x.square().double() - FieldElement::one();
    }
    point.x
}
/// Structure to iterate over the powers of a `FieldElement`.
// Need to se where it can be defined
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
