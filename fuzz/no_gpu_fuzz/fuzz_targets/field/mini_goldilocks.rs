#![no_main]

use libfuzzer_sys::fuzz_target;
use lambdaworks_math::field::{
    element::FieldElement, 
    fields::{
        u64_prime_field::U64FieldElement,
        u64_goldilocks_field::Goldilocks64Field,
    },
};
use p3_goldilocks::Goldilocks;
use p3_field::{Field, PrimeField64, AbstractField};

fuzz_target!(|values: (u64, u64)| {

    let (value_u64_a, value_u64_b) = values;

    let a =  FieldElement::<Goldilocks64Field>::from(value_u64_a);
    let b =  FieldElement::<Goldilocks64Field>::from(value_u64_b);

    let a_expected = Goldilocks::from_canonical_u64(value_u64_a);
    let b_expected = Goldilocks::from_canonical_u64(value_u64_b);

    let add_u64 = &a + &b;
    let addition = a_expected + b_expected;

    assert_eq!(add_u64.representative(), addition.as_canonical_u64());

    let sub_u64 = &a - &b;
    let substraction = a_expected - b_expected;
    assert_eq!(sub_u64.representative(), substraction.as_canonical_u64());

    let mul_u64 = &a * &b;
    let multiplication = a_expected  * b_expected;
    assert_eq!(mul_u64.representative(), multiplication.as_canonical_u64());

    let pow = &a.pow(b.representative());
    let expected_pow = a_expected.exp_u64(b_expected.as_canonical_u64());
    assert_eq!(pow.representative(), expected_pow.as_canonical_u64());

    if value_u64_b != 0 && b.inv().is_ok() && b_expected.try_inverse().is_some() { 

        let div = &a / &b; 
        assert_eq!(&div * &b, a.clone());
        let expected_div = a_expected / b_expected;
        assert_eq!(div.representative(), expected_div.as_canonical_u64());
    }

    for n in [&a, &b] {
        match n.sqrt() {
            Some((fst_sqrt, snd_sqrt)) => {
                assert_eq!(fst_sqrt.square(), snd_sqrt.square(), "Squared roots don't match each other");
                assert_eq!(n, &fst_sqrt.square(), "Squared roots don't match original number");
            }
            None => {}
        };
    }

    // Axioms soundness

    let one = FieldElement::<Goldilocks64Field>::one();
    let zero = FieldElement::<Goldilocks64Field>::zero();

    assert_eq!(&a + &zero, a, "Neutral add element a failed");
    assert_eq!(&b + &zero, b, "Neutral mul element b failed");
    assert_eq!(&a * &one, a, "Neutral add element a failed");
    assert_eq!(&b * &one, b, "Neutral mul element b failed");

    assert_eq!(&a + &b, &b + &a, "Commutative add property failed");
    assert_eq!(&a * &b, &b * &a, "Commutative mul property failed");

    let c = &a * &b;
    assert_eq!((&a + &b) + &c, &a + (&b + &c), "Associative add property failed");
    assert_eq!((&a * &b) * &c, &a * (&b * &c), "Associative mul property failed");

    assert_eq!(&a * (&b + &c), &a * &b + &a * &c, "Distributive property failed");

    assert_eq!(&a - &a, zero, "Inverse add a failed");
    assert_eq!(&b - &b, zero, "Inverse add b failed");

    if a != zero {
        assert_eq!(&a * a.inv().unwrap(), one, "Inverse mul a failed");
    }
    if b != zero {
        assert_eq!(&b * b.inv().unwrap(), one, "Inverse mul b failed");
    }
});
