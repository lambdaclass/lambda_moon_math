#![no_main]

use lambdaworks_math::field::{
    element::FieldElement,
    fields::mersenne31::field::{Mersenne31Field, MERSENNE_31_PRIME_FIELD_ORDER},
};
use libfuzzer_sys::fuzz_target;
use p3_field::{AbstractField, Field, PrimeField32, PrimeField64};
use p3_mersenne_31::Mersenne31;

fuzz_target!(|values: (u32, u32)| {
    // Note: we filter values outside of order as it triggers an assert within plonky3 disallowing values n >= Self::Order
    if values.0 >= MERSENNE_31_PRIME_FIELD_ORDER || values.1 >= MERSENNE_31_PRIME_FIELD_ORDER {
        return;
    }

    let (value_u32_a, value_u32_b) = values;

    let a = FieldElement::<Mersenne31Field>::from(value_u32_a as u64);
    let b = FieldElement::<Mersenne31Field>::from(value_u32_b as u64);

    // Note: if we parse using from_canonical_u32 fails due to check that n < Self::Order
    let a_expected = Mersenne31::from_canonical_u32(value_u32_a);
    let b_expected = Mersenne31::from_canonical_u32(value_u32_b);

    let add_u32 = &a + &b;
    let addition = a_expected + b_expected;

    assert_eq!(add_u32.representative(), addition.as_canonical_u32());

    let sub_u32 = &a - &b;
    let substraction = a_expected - b_expected;
    assert_eq!(sub_u32.representative(), substraction.as_canonical_u32());

    let mul_u32 = &a * &b;
    let multiplication = a_expected * b_expected;
    assert_eq!(mul_u32.representative(), multiplication.as_canonical_u32());

    let pow = &a.pow(b.representative());
    let expected_pow = a_expected.exp_u64(b_expected.as_canonical_u64());
    assert_eq!(pow.representative(), expected_pow.as_canonical_u32());

    if value_u32_b != 0 && b.inv().is_ok() && b_expected.try_inverse().is_some() {
        let div = &a / &b;
        assert_eq!(&div * &b, a.clone());
        let expected_div = a_expected / b_expected;
        assert_eq!(div.representative(), expected_div.as_canonical_u32());
    }

    for n in [&a, &b] {
        match n.sqrt() {
            Some((fst_sqrt, snd_sqrt)) => {
                assert_eq!(
                    fst_sqrt.square(),
                    snd_sqrt.square(),
                    "Squared roots don't match each other"
                );
                assert_eq!(
                    n,
                    &fst_sqrt.square(),
                    "Squared roots don't match original number"
                );
            }
            None => {}
        };
    }

    // Axioms soundness

    let one = FieldElement::<Mersenne31Field>::one();
    let zero = FieldElement::<Mersenne31Field>::zero();

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
