#![no_main]

use ibig::{modular::ModuloRing, UBig};
use lambdaworks_math::field::{
    element::FieldElement, fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
};
use lambdaworks_math::unsigned_integer::element::UnsignedInteger;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|values: (u64, u64)| {
    let (value_u64_a, value_u64_b) = values;
    let prime = UBig::from_str_radix(
        "800000000000011000000000000000000000000000000000000000000000001",
        16,
    )
    .unwrap();
    let ring = ModuloRing::new(&prime);

    let value_a = UnsignedInteger::from(value_u64_a);
    let value_b = UnsignedInteger::from(value_u64_b);

    let a = FieldElement::<Stark252PrimeField>::from_raw(value_a);
    let b = FieldElement::<Stark252PrimeField>::from_raw(value_b);

    let _a_expected = ring.from(value_u64_a);
    let _b_expected = ring.from(value_u64_b);

    let _add_u64 = &a + &b;

    let _sub_u64 = &a - &b;

    let _mul_u64 = &a * &b;

    let _pow = &a.pow(b.representative());

    if value_u64_b != 0 {
        let div = &a / &b;
        assert_eq!(&div * &b, a.clone());
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

    let one = FieldElement::<Stark252PrimeField>::one();
    let zero = FieldElement::<Stark252PrimeField>::zero();

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
