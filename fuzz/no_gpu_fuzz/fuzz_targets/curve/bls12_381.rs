#![no_main]

use lambdaworks_math::{
    cyclic_group::IsGroup,
    elliptic_curve::{
        short_weierstrass::{
            curves::bls12_381::{
                curve::BLS12381Curve, field_extension::Degree12ExtensionField,
                pairing::BLS12381AtePairing, twist::BLS12381TwistCurve,
            },
            point::ShortWeierstrassProjectivePoint,
        },
        traits::{IsEllipticCurve, IsPairing},
    },
    field::element::FieldElement,
    unsigned_integer::element::U384,
};
use libfuzzer_sys::fuzz_target;

type LambdaG1 = ShortWeierstrassProjectivePoint<BLS12381Curve>;
type LambdaG2 = ShortWeierstrassProjectivePoint<BLS12381TwistCurve>;

//TODO: Derive arbitrary for Affine and Projective or change this to use &[u8] as input to cover more cases.
//TODO: Use more advanced options to generate values over curve specifically given most inputs will fail curve check.
//TODO: Investigate normalization of projective coordinates in arkworks to allow for differential fuzzing.
fuzz_target!(|values: (u64, u64)| {
    let (a_val, b_val) = values;

    let a_g1 = BLS12381Curve::generator().operate_with_self(a_val);
    let b_g1 = BLS12381Curve::generator().operate_with_self(b_val);

    let a_g2 = BLS12381TwistCurve::generator().operate_with_self(a_val);
    let b_g2 = BLS12381TwistCurve::generator().operate_with_self(b_val);

    // ***AXIOM SOUNDNESS***
    let g1_zero = LambdaG1::neutral_element();

    let g2_zero = LambdaG2::neutral_element();

    // G1
    // -O = O
    assert_eq!(g1_zero.neg(), g1_zero, "Neutral mul element a failed");

    // P + O = O
    assert_eq!(
        a_g1.operate_with(&g1_zero),
        a_g1,
        "Neutral operate_with element a failed"
    );
    assert_eq!(
        b_g1.operate_with(&g1_zero),
        b_g1,
        "Neutral operate_with element b failed"
    );

    // P + Q = Q + P
    assert_eq!(
        a_g1.operate_with(&b_g1),
        b_g1.operate_with(&a_g1),
        "Commutative add property failed"
    );

    // (P + Q) + R = Q + (P + R)
    let c_g1 = a_g1.operate_with(&b_g1);
    assert_eq!(
        (a_g1.operate_with(&b_g1)).operate_with(&c_g1),
        a_g1.operate_with(&b_g1.operate_with(&c_g1)),
        "Associative operate_with property failed"
    );

    // P + -P = O
    assert_eq!(
        a_g1.operate_with(&a_g1.neg()),
        g1_zero,
        "Inverse add a failed"
    );
    assert_eq!(
        b_g1.operate_with(&b_g1.neg()),
        g1_zero,
        "Inverse add b failed"
    );

    // G2
    // -O = O
    assert_eq!(g2_zero.neg(), g2_zero, "Neutral mul element a failed");

    // P + O = O
    assert_eq!(
        a_g2.operate_with(&g2_zero),
        a_g2,
        "Neutral operate_with element a failed"
    );
    assert_eq!(
        b_g2.operate_with(&g2_zero),
        b_g2,
        "Neutral operate_with element b failed"
    );

    // P + Q = Q + P
    assert_eq!(
        a_g2.operate_with(&b_g2),
        b_g2.operate_with(&a_g2),
        "Commutative add property failed"
    );

    // (P + Q) + R = Q + (P + R)
    let c_g2 = a_g2.operate_with(&b_g2);
    assert_eq!(
        (a_g2.operate_with(&b_g2)).operate_with(&c_g2),
        a_g2.operate_with(&b_g2.operate_with(&c_g2)),
        "Associative operate_with property failed"
    );

    // P + -P = O
    assert_eq!(
        a_g2.operate_with(&a_g2.neg()),
        g2_zero,
        "Inverse add a failed"
    );
    assert_eq!(
        b_g2.operate_with(&b_g2.neg()),
        g2_zero,
        "Inverse add b failed"
    );

    // Pairing Bilinearity
    let a = U384::from_u64(a_val);
    let b = U384::from_u64(b_val);
    let result = BLS12381AtePairing::compute_batch(&[
        (
            &a_g1.operate_with_self(a).to_affine(),
            &a_g2.operate_with_self(b).to_affine(),
        ),
        (
            &a_g1.operate_with_self(a * b).to_affine(),
            &a_g2.neg().to_affine(),
        ),
    ]);
    assert_eq!(result, Ok(FieldElement::<Degree12ExtensionField>::one()));

    // Ate Pairing returns one with one element is neutral element
    let result =
        BLS12381AtePairing::compute_batch(&[(&a_g1.to_affine(), &LambdaG2::neutral_element())]);
    assert_eq!(result, Ok(FieldElement::<Degree12ExtensionField>::one()));

    let result =
        BLS12381AtePairing::compute_batch(&[(&LambdaG1::neutral_element(), &a_g2.to_affine())]);
    assert_eq!(result, Ok(FieldElement::<Degree12ExtensionField>::one()));
});
