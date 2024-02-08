use super::curve::MILLER_LOOP_CONSTANT;
use super::default_types::FrField;
use super::{
    curve::BLS12381Curve,
    field_extension::{BLS12381PrimeField, Degree12ExtensionField, Degree2ExtensionField},
    twist::BLS12381TwistCurve,
};
use crate::{cyclic_group::IsGroup, elliptic_curve::traits::IsPairing, errors::PairingError};

use crate::{
    elliptic_curve::short_weierstrass::{
        curves::bls12_381::field_extension::{Degree6ExtensionField, LevelTwoResidue},
        point::ShortWeierstrassProjectivePoint,
        traits::IsShortWeierstrass,
    },
    field::{element::FieldElement, extensions::cubic::HasCubicNonResidue},
    unsigned_integer::element::{UnsignedInteger, U256},
};

pub const SUBGROUP_ORDER: U256 =
    U256::from_hex_unchecked("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001");

#[derive(Clone)]
pub struct BLS12381AtePairing;

impl IsPairing for BLS12381AtePairing {
    type G1Point = ShortWeierstrassProjectivePoint<BLS12381Curve>;
    type G2Point = ShortWeierstrassProjectivePoint<BLS12381TwistCurve>;
    type OutputField = Degree12ExtensionField;
    type BaseField = FrField;

    /// Compute the product of the ate pairings for a list of point pairs.
    fn compute_batch(
        pairs: &[(&Self::G1Point, &Self::G2Point)],
    ) -> Result<FieldElement<Self::OutputField>, PairingError> {
        let mut result = FieldElement::one();
        for (p, q) in pairs {
            if !p.is_in_subgroup() || !q.is_in_subgroup() {
                return Err(PairingError::PointNotInSubgroup);
            }
            if !p.is_neutral_element() && !q.is_neutral_element() {
                let p = p.to_affine();
                let q = q.to_affine();
                result *= miller(&q, &p);
            }
        }
        Ok(final_exponentiation(&result))
    }
}

fn double_accumulate_line(
    t: &mut ShortWeierstrassProjectivePoint<BLS12381TwistCurve>,
    p: &ShortWeierstrassProjectivePoint<BLS12381Curve>,
    accumulator: &mut FieldElement<Degree12ExtensionField>,
) {
    let [x1, y1, z1] = t.coordinates();
    let [px, py, _] = p.coordinates();
    let residue = LevelTwoResidue::residue();
    let two_inv = FieldElement::<Degree2ExtensionField>::new_base("d0088f51cbff34d258dd3db21a5d66bb23ba5c279c2895fb39869507b587b120f55ffff58a9ffffdcff7fffffffd556");
    let three = FieldElement::<BLS12381PrimeField>::from(3);

    let a = &two_inv * x1 * y1;
    let b = y1.square();
    let c = z1.square();
    let d = &three * &c;
    let e = BLS12381TwistCurve::b() * d;
    let f = &three * &e;
    let g = two_inv * (&b + &f);
    let h = (y1 + z1).square() - (&b + &c);

    let x3 = &a * (&b - &f);
    let y3 = g.square() - (&three * e.square());
    let z3 = &b * &h;

    let [h0, h1] = h.value();
    let x1_sq_3 = three * x1.square();
    let [x1_sq_30, x1_sq_31] = x1_sq_3.value();

    t.0.value = [x3, y3, z3];

    // (a0 + a2w2 + a4w4 + a1w + a3w3 + a5w5) * (b0 + b2 w2 + b3 w3) =
    // (a0b0 + r (a3b3 + a4b2)) w0 + (a1b0 + r (a4b3 + a5b2)) w
    // (a2b0 + r  a5b3 + a0b2 ) w2 + (a3b0 + a0b3 + a1b2    ) w3
    // (a4b0 +    a1b3 + a2b2 ) w4 + (a5b0 + a2b3 + a3b2    ) w5
    let accumulator_sq = accumulator.square();
    let [x, y] = accumulator_sq.value();
    let [a0, a2, a4] = x.value();
    let [a1, a3, a5] = y.value();
    let b0 = e - b;
    let b2 = FieldElement::new([x1_sq_30 * px, x1_sq_31 * px]);
    let b3 = FieldElement::<Degree2ExtensionField>::new([-h0 * py, -h1 * py]);
    *accumulator = FieldElement::new([
        FieldElement::new([
            a0 * &b0 + &residue * (a3 * &b3 + a4 * &b2), // w0
            a2 * &b0 + &residue * a5 * &b3 + a0 * &b2,   // w2
            a4 * &b0 + a1 * &b3 + a2 * &b2,              // w4
        ]),
        FieldElement::new([
            a1 * &b0 + &residue * (a4 * &b3 + a5 * &b2), // w1
            a3 * &b0 + a0 * &b3 + a1 * &b2,              // w3
            a5 * &b0 + a2 * &b3 + a3 * &b2,              // w5
        ]),
    ]);
}

fn add_accumulate_line(
    t: &mut ShortWeierstrassProjectivePoint<BLS12381TwistCurve>,
    q: &ShortWeierstrassProjectivePoint<BLS12381TwistCurve>,
    p: &ShortWeierstrassProjectivePoint<BLS12381Curve>,
    accumulator: &mut FieldElement<Degree12ExtensionField>,
) {
    let [x1, y1, z1] = t.coordinates();
    let [x2, y2, _] = q.coordinates();
    let [px, py, _] = p.coordinates();
    let residue = LevelTwoResidue::residue();

    let a = y2 * z1;
    let b = x2 * z1;
    let theta = y1 - a;
    let lambda = x1 - b;
    let c = theta.square();
    let d = lambda.square();
    let e = &lambda * &d;
    let f = z1 * c;
    let g = x1 * d;
    let h = &e + f - FieldElement::<BLS12381PrimeField>::from(2) * &g;
    let i = y1 * &e;

    let x3 = &lambda * &h;
    let y3 = &theta * (g - h) - i;
    let z3 = z1 * e;

    t.0.value = [x3, y3, z3];

    let [lambda0, lambda1] = lambda.value();
    let [theta0, theta1] = theta.value();

    let [x, y] = accumulator.value();
    let [a0, a2, a4] = x.value();
    let [a1, a3, a5] = y.value();
    let b0 = -lambda.clone() * y2 + theta.clone() * x2;
    let b2 = FieldElement::new([-theta0 * px, -theta1 * px]);
    let b3 = FieldElement::<Degree2ExtensionField>::new([lambda0 * py, lambda1 * py]);
    *accumulator = FieldElement::new([
        FieldElement::new([
            a0 * &b0 + &residue * (a3 * &b3 + a4 * &b2), // w0
            a2 * &b0 + &residue * a5 * &b3 + a0 * &b2,   // w2
            a4 * &b0 + a1 * &b3 + a2 * &b2,              // w4
        ]),
        FieldElement::new([
            a1 * &b0 + &residue * (a4 * &b3 + a5 * &b2), // w1
            a3 * &b0 + a0 * &b3 + a1 * &b2,              // w3
            a5 * &b0 + a2 * &b3 + a3 * &b2,              // w5
        ]),
    ]);
}
/// Implements the miller loop for the ate pairing of the BLS12 381 curve.
/// Based on algorithm 9.2, page 212 of the book
/// "Topics in computational number theory" by W. Bons and K. Lenstra
#[allow(unused)]
fn miller(
    q: &ShortWeierstrassProjectivePoint<BLS12381TwistCurve>,
    p: &ShortWeierstrassProjectivePoint<BLS12381Curve>,
) -> FieldElement<Degree12ExtensionField> {
    let mut r = q.clone();
    let mut f = FieldElement::<Degree12ExtensionField>::one();
    let mut miller_loop_constant = MILLER_LOOP_CONSTANT;
    let mut miller_loop_constant_bits: alloc::vec::Vec<bool> = alloc::vec![];

    while miller_loop_constant > 0 {
        miller_loop_constant_bits.insert(0, (miller_loop_constant & 1) == 1);
        miller_loop_constant >>= 1;
    }

    for bit in miller_loop_constant_bits[1..].iter() {
        double_accumulate_line(&mut r, p, &mut f);
        if *bit {
            add_accumulate_line(&mut r, q, p, &mut f);
        }
    }
    f.inv().unwrap()
}

/// Auxiliary function for the final exponentiation of the ate pairing.
fn frobenius_square(
    f: &FieldElement<Degree12ExtensionField>,
) -> FieldElement<Degree12ExtensionField> {
    let [a, b] = f.value();
    let w_raised_to_p_squared_minus_one = FieldElement::<Degree6ExtensionField>::new_base("1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad");
    let omega_3 = FieldElement::<Degree2ExtensionField>::new_base("1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac");
    let omega_3_squared = FieldElement::<Degree2ExtensionField>::new_base(
        "5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe",
    );

    let [a0, a1, a2] = a.value();
    let [b0, b1, b2] = b.value();

    let f0 = FieldElement::new([a0.clone(), a1 * &omega_3, a2 * &omega_3_squared]);
    let f1 = FieldElement::new([b0.clone(), b1 * omega_3, b2 * omega_3_squared]);

    FieldElement::new([f0, w_raised_to_p_squared_minus_one * f1])
}

// To understand more about how to reduce the final exponentiation
// read "Efficient Final Exponentiation via Cyclotomic Structure for
// Pairings over Families of Elliptic Curves" (https://eprint.iacr.org/2020/875.pdf)
//
// TODO: implement optimizations for the hard part of the final exponentiation.
#[allow(unused)]
fn final_exponentiation(
    base: &FieldElement<Degree12ExtensionField>,
) -> FieldElement<Degree12ExtensionField> {
    const PHI_DIVIDED_BY_R: UnsignedInteger<20> = UnsignedInteger::from_hex_unchecked("f686b3d807d01c0bd38c3195c899ed3cde88eeb996ca394506632528d6a9a2f230063cf081517f68f7764c28b6f8ae5a72bce8d63cb9f827eca0ba621315b2076995003fc77a17988f8761bdc51dc2378b9039096d1b767f17fcbde783765915c97f36c6f18212ed0b283ed237db421d160aeb6a1e79983774940996754c8c71a2629b0dea236905ce937335d5b68fa9912aae208ccf1e516c3f438e3ba79");

    let f1 = base.conjugate() * base.inv().unwrap();
    let f2 = frobenius_square(&f1) * f1;
    f2.pow(PHI_DIVIDED_BY_R)
}

#[cfg(test)]
mod tests {
    use crate::{
        cyclic_group::IsGroup, elliptic_curve::traits::IsEllipticCurve,
        unsigned_integer::element::U384,
    };

    use super::*;

    #[test]
    fn test_double_accumulate_line_doubles_point_correctly() {
        let g1 = BLS12381Curve::generator();
        let g2 = BLS12381TwistCurve::generator();
        let mut r = g2.clone();
        let mut f = FieldElement::one();
        double_accumulate_line(&mut r, &g1, &mut f);
        assert_eq!(r, g2.operate_with(&g2));
    }

    #[test]
    fn test_add_accumulate_line_adds_points_correctly() {
        let g1 = BLS12381Curve::generator();
        let g = BLS12381TwistCurve::generator();
        let a: u64 = 12;
        let b: u64 = 23;
        let g2 = g.operate_with_self(a).to_affine();
        let g3 = g.operate_with_self(b).to_affine();
        let expected = g.operate_with_self(a + b);
        let mut r = g2;
        let mut f = FieldElement::one();
        add_accumulate_line(&mut r, &g3, &g1, &mut f);
        assert_eq!(r, expected);
    }

    #[test]
    fn batch_ate_pairing_bilinearity() {
        let p = BLS12381Curve::generator();
        let q = BLS12381TwistCurve::generator();
        let a = U384::from_u64(11);
        let b = U384::from_u64(93);

        let result = BLS12381AtePairing::compute_batch(&[
            (
                &p.operate_with_self(a).to_affine(),
                &q.operate_with_self(b).to_affine(),
            ),
            (
                &p.operate_with_self(a * b).to_affine(),
                &q.neg().to_affine(),
            ),
        ])
        .unwrap();
        assert_eq!(result, FieldElement::one());
    }

    #[test]
    fn ate_pairing_returns_one_when_one_element_is_the_neutral_element() {
        let p = BLS12381Curve::generator().to_affine();
        let q = ShortWeierstrassProjectivePoint::neutral_element();
        let result = BLS12381AtePairing::compute_batch(&[(&p.to_affine(), &q)]).unwrap();
        assert_eq!(result, FieldElement::one());

        let p = ShortWeierstrassProjectivePoint::neutral_element();
        let q = BLS12381TwistCurve::generator();
        let result = BLS12381AtePairing::compute_batch(&[(&p, &q.to_affine())]).unwrap();
        assert_eq!(result, FieldElement::one());
    }

    #[test]
    fn ate_pairing_errors_when_one_element_is_not_in_subgroup() {
        let p = ShortWeierstrassProjectivePoint::new([
            FieldElement::one(),
            FieldElement::one(),
            FieldElement::one(),
        ]);
        let q = ShortWeierstrassProjectivePoint::neutral_element();
        let result = BLS12381AtePairing::compute_batch(&[(&p.to_affine(), &q)]);
        assert!(result.is_err())
    }
}
