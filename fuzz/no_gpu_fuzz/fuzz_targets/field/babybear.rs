#![no_main]

use lambdaworks_math::field::{
    element::FieldElement,
    fields::u32_montgomery_backend_prime_field::U32MontgomeryBackendPrimeField,
};
use libfuzzer_sys::fuzz_target;
use p3_baby_bear::BabyBear;
use p3_field::{Field, FieldAlgebra, PrimeField32};

pub type U32Babybear31PrimeField = U32MontgomeryBackendPrimeField<2013265921>;
pub type F = FieldElement<U32Babybear31PrimeField>;

fuzz_target!(|values: (u32, u32)| {
    // Note: we filter values outside of order as it triggers an assert within plonky3 disallowing values n >= Self::Order
    let (value_u32_a, value_u32_b) = values;

    if value_u32_a >= 2013265921 || value_u32_b >= 2013265921 {
        return;
    }
    let a = F::from(value_u32_a as u64);
    let b = F::from(value_u32_b as u64);

    // Note: if we parse using from_canonical_u32 fails due to check that n < Self::Order
    let a_expected = BabyBear::from_canonical_u32(value_u32_a);
    let b_expected = BabyBear::from_canonical_u32(value_u32_b);

    let add_u32 = &a + &b;
    let addition = a_expected + b_expected;
    assert_eq!(add_u32.representative(), addition.as_canonical_u32());

    let sub_u32 = &a - &b;
    let substraction = a_expected - b_expected;
    assert_eq!(sub_u32.representative(), substraction.as_canonical_u32());

    let mul_u32 = &a * &b;
    let multiplication = a_expected * b_expected;
    assert_eq!(mul_u32.representative(), multiplication.as_canonical_u32());

    // Axioms soundness
    let one = F::one();
    let zero = F::zero();

    assert_eq!(&a + &zero, a, "Neutral add element a failed");
    assert_eq!(&b + &zero, b, "Neutral mul element b failed");
    assert_eq!(&a * &one, a, "Neutral add element a failed");
    assert_eq!(&b * &one, b, "Neutral mul element b failed");

    assert_eq!(&a + &b, &b + &a, "Commutative add property failed");
    assert_eq!(&a * &b, &b * &a, "Commutative mul property failed");

    let c = &a * &b;
    assert_eq!(
        (&a + &b) + &c,
        &a + (&b + &c),
        "Associative add property failed"
    );
    assert_eq!(
        (&a * &b) * &c,
        &a * (&b * &c),
        "Associative mul property failed"
    );

    assert_eq!(
        &a * (&b + &c),
        &a * &b + &a * &c,
        "Distributive property failed"
    );

    assert_eq!(&a - &a, zero, "Inverse add a failed");
    assert_eq!(&b - &b, zero, "Inverse add b failed");

    if a != zero {
        assert_eq!(&a * a.inv().unwrap(), one, "Inverse mul a failed");
    }
    if b != zero {
        assert_eq!(&b * b.inv().unwrap(), one, "Inverse mul b failed");
    }
});
