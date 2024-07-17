#![no_main]

use libfuzzer_sys::fuzz_target;
use lambdaworks_math::{
    cyclic_group::IsGroup,
    elliptic_curve::{
        traits::{IsEllipticCurve, IsPairing},
        short_weierstrass::{
            curves::grumpkin::curve::GrumpkinCurve,
            point::ShortWeierstrassProjectivePoint,
        }
    },
    field::element::FieldElement,
};

type LambdaG1 = ShortWeierstrassProjectivePoint<GrumpkinCurve>;

//TODO: derive arbitrary for Affine and Projective or change this to use &[u8] as input to cover more cases
fuzz_target!(|values: (u64, u64)| {
    let (a_val, b_val) = values;

    let a_g1 = GrumpkinCurve::generator().operate_with_self(a_val);
    let b_g1 = GrumpkinCurve::generator().operate_with_self(b_val);

    // ***AXIOM SOUNDNESS***
    let g1_zero = LambdaG1::neutral_element();

    // -O = O
    assert_eq!(g1_zero.neg(), g1_zero, "Neutral mul element a failed");

    // P * O = O
    assert_eq!(a_g1.operate_with(&g1_zero), a_g1, "Neutral operate_with element a failed");
    assert_eq!(b_g1.operate_with(&g1_zero), b_g1, "Neutral operate_with element b failed");

    // P * Q = Q * P
    assert_eq!(a_g1.operate_with(&b_g1), b_g1.operate_with(&a_g1), "Commutative add property failed");

    // (P * Q) * R = Q * (P * R)
    let c_g1 = a_g1.operate_with(&b_g1);
    assert_eq!((a_g1.operate_with(&b_g1)).operate_with(&c_g1), a_g1.operate_with(&b_g1.operate_with(&c_g1)), "Associative operate_with property failed");

    // P * -P = O
    assert_eq!(a_g1.operate_with(&a_g1.neg()), g1_zero, "Inverse add a failed");
    assert_eq!(b_g1.operate_with(&b_g1.neg()), g1_zero, "Inverse add b failed");
});
