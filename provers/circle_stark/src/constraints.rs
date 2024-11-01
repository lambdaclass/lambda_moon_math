use lambdaworks_math::{
    circle::point::CirclePoint,
    field::{element::FieldElement, fields::mersenne31::field::Mersenne31Field},
};

/// Evaluate the vanishing polynomial of the standard coset of size 2^log_2_size in a point.
/// The vanishing polynomial of a coset is the polynomial that takes the value zero when evaluated
/// in all the points of the coset.
/// We are using that if we take a point in g_{2n} + <g_n> and double it n-1 times, then
/// we'll get the point (0, 1) or (0, -1); so its coordinate x is always 0.
pub fn evaluate_vanishing_poly(
    log_2_size: u32,
    point: &CirclePoint<Mersenne31Field>,
) -> FieldElement<Mersenne31Field> {
    let mut x = point.x;
    for _ in 1..log_2_size {
        x = x.square().double() - FieldElement::one();
    }
    x
}

// Evaluate the polynomial that vanishes at a specific point in the domain at an arbitrary point.
// This use the "tangent" line to the domain in the vanish point.
// Check: https://vitalik.eth.limo/general/2024/07/23/circlestarks.html for details.
pub fn evaluate_single_point_zerofier(
    vanish_point: CirclePoint<Mersenne31Field>,
    eval_point: &CirclePoint<Mersenne31Field>,
) -> FieldElement<Mersenne31Field> {
    (eval_point + vanish_point.conjugate()).x - FieldElement::<Mersenne31Field>::one()
}

#[cfg(test)]
mod tests {
    use super::*;
    use lambdaworks_math::circle::cosets::Coset;

    type FE = FieldElement<Mersenne31Field>;

    #[test]
    fn vanishing_poly_vanishes_in_coset() {
        let log_2_size = 3;
        let coset = Coset::new_standard(log_2_size);
        let points = Coset::get_coset_points(&coset);
        for point in points {
            assert_eq!(evaluate_vanishing_poly(log_2_size, &point), FE::zero());
        }
    }
    #[test]
    fn vanishing_poly_doesnt_vanishes_outside_coset() {
        let log_2_size = 3;
        let coset = Coset::new_standard(log_2_size + 1);
        let points = Coset::get_coset_points(&coset);
        for point in points {
            assert_ne!(evaluate_vanishing_poly(log_2_size, &point), FE::zero());
        }
    }

    #[test]
    fn single_point_zerofier_vanishes_only_in_vanish_point() {
        let vanish_point = CirclePoint::GENERATOR;
        let eval_point = &vanish_point * 3;
        assert_eq!(
            evaluate_single_point_zerofier(vanish_point.clone(), &vanish_point),
            FE::zero()
        );
        assert_ne!(
            evaluate_single_point_zerofier(vanish_point.clone(), &eval_point),
            FE::zero()
        );
    }
}
