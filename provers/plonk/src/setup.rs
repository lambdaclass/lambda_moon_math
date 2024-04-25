use std::collections::HashMap;

use crate::constraint_system::{get_permutation, ConstraintSystem, Variable};
use crate::test_utils::utils::{generate_domain, generate_permutation_coefficients};
use lambdaworks_crypto::commitments::traits::IsCommitmentScheme;
use lambdaworks_crypto::fiat_shamir::{
    default_transcript::DefaultTranscript, is_transcript::IsTranscript,
};
use lambdaworks_math::field::traits::IsFFTField;
use lambdaworks_math::field::{element::FieldElement, traits::IsField};
use lambdaworks_math::polynomial::Polynomial;
use lambdaworks_math::traits::{AsBytes, ByteConversion};

//TODO: feature gate
use lambdaworks_math::gpu::icicle::IcicleFFT;

// TODO: implement getters
pub struct Witness<F: IsField> {
    pub a: Vec<FieldElement<F>>,
    pub b: Vec<FieldElement<F>>,
    pub c: Vec<FieldElement<F>>,
}

impl<F: IsField> Witness<F> {
    pub fn new(values: HashMap<Variable, FieldElement<F>>, system: &ConstraintSystem<F>) -> Self {
        let (lro, _) = system.to_matrices();
        let abc: Vec<_> = lro.iter().map(|v| values[v].clone()).collect();
        let n = lro.len() / 3;

        Self {
            a: abc[..n].to_vec(),
            b: abc[n..2 * n].to_vec(),
            c: abc[2 * n..].to_vec(),
        }
    }
}

// TODO: implement getters
#[derive(Clone)]
pub struct CommonPreprocessedInput<F: IsField + IcicleFFT>
where
    FieldElement<F>: ByteConversion,
{
    pub n: usize,
    /// Number of constraints
    pub domain: Vec<FieldElement<F>>,
    pub omega: FieldElement<F>,
    pub k1: FieldElement<F>,

    pub ql: Polynomial<FieldElement<F>>,
    pub qr: Polynomial<FieldElement<F>>,
    pub qo: Polynomial<FieldElement<F>>,
    pub qm: Polynomial<FieldElement<F>>,
    pub qc: Polynomial<FieldElement<F>>,

    pub s1: Polynomial<FieldElement<F>>,
    pub s2: Polynomial<FieldElement<F>>,
    pub s3: Polynomial<FieldElement<F>>,

    pub s1_lagrange: Vec<FieldElement<F>>,
    pub s2_lagrange: Vec<FieldElement<F>>,
    pub s3_lagrange: Vec<FieldElement<F>>,
}

impl<F: IsFFTField + IcicleFFT> CommonPreprocessedInput<F>
where
    FieldElement<F>: ByteConversion,
{
    pub fn from_constraint_system(
        system: &ConstraintSystem<F>,
        order_r_minus_1_root_unity: &FieldElement<F>,
    ) -> Self {
        let (lro, q) = system.to_matrices();
        let n = lro.len() / 3;
        let omega = F::get_primitive_root_of_unity(n.trailing_zeros() as u64).unwrap();
        let domain = generate_domain(&omega, n);

        let m = q.len() / 5;
        let ql: Vec<_> = q[..m].to_vec();
        let qr: Vec<_> = q[m..2 * m].to_vec();
        let qm: Vec<_> = q[2 * m..3 * m].to_vec();
        let qo: Vec<_> = q[3 * m..4 * m].to_vec();
        let qc: Vec<_> = q[4 * m..].to_vec();

        let permutation = get_permutation(&lro);
        let permuted =
            generate_permutation_coefficients(&omega, n, &permutation, order_r_minus_1_root_unity);

        let s1_lagrange: Vec<_> = permuted[..n].to_vec();
        let s2_lagrange: Vec<_> = permuted[n..2 * n].to_vec();
        let s3_lagrange: Vec<_> = permuted[2 * n..].to_vec();

        Self {
            domain,
            n,
            omega,
            k1: order_r_minus_1_root_unity.clone(),
            ql: Polynomial::interpolate_fft::<F>(&ql).unwrap(), // TODO: Remove unwraps
            qr: Polynomial::interpolate_fft::<F>(&qr).unwrap(),
            qo: Polynomial::interpolate_fft::<F>(&qo).unwrap(),
            qm: Polynomial::interpolate_fft::<F>(&qm).unwrap(),
            qc: Polynomial::interpolate_fft::<F>(&qc).unwrap(),
            s1: Polynomial::interpolate_fft::<F>(&s1_lagrange).unwrap(),
            s2: Polynomial::interpolate_fft::<F>(&s2_lagrange).unwrap(),
            s3: Polynomial::interpolate_fft::<F>(&s3_lagrange).unwrap(),
            s1_lagrange,
            s2_lagrange,
            s3_lagrange,
        }
    }
}

pub struct VerificationKey<G1Point> {
    pub qm_1: G1Point,
    pub ql_1: G1Point,
    pub qr_1: G1Point,
    pub qo_1: G1Point,
    pub qc_1: G1Point,

    pub s1_1: G1Point,
    pub s2_1: G1Point,
    pub s3_1: G1Point,
}

pub fn setup<F: IsField + IcicleFFT, CS: IsCommitmentScheme<F>>(
    common_input: &CommonPreprocessedInput<F>,
    commitment_scheme: &CS,
) -> VerificationKey<CS::Commitment>
where
    FieldElement<F>: ByteConversion,
{
    VerificationKey {
        qm_1: commitment_scheme.commit(&common_input.qm),
        ql_1: commitment_scheme.commit(&common_input.ql),
        qr_1: commitment_scheme.commit(&common_input.qr),
        qo_1: commitment_scheme.commit(&common_input.qo),
        qc_1: commitment_scheme.commit(&common_input.qc),

        s1_1: commitment_scheme.commit(&common_input.s1),
        s2_1: commitment_scheme.commit(&common_input.s2),
        s3_1: commitment_scheme.commit(&common_input.s3),
    }
}

pub fn new_strong_fiat_shamir_transcript<F, CS>(
    vk: &VerificationKey<CS::Commitment>,
    public_input: &[FieldElement<F>],
) -> DefaultTranscript<F>
where
    F: IsField,
    FieldElement<F>: ByteConversion,
    CS: IsCommitmentScheme<F>,
    CS::Commitment: AsBytes,
{
    let mut transcript = DefaultTranscript::default();

    transcript.append_bytes(&vk.s1_1.as_bytes());
    transcript.append_bytes(&vk.s2_1.as_bytes());
    transcript.append_bytes(&vk.s3_1.as_bytes());
    transcript.append_bytes(&vk.ql_1.as_bytes());
    transcript.append_bytes(&vk.qr_1.as_bytes());
    transcript.append_bytes(&vk.qm_1.as_bytes());
    transcript.append_bytes(&vk.qo_1.as_bytes());
    transcript.append_bytes(&vk.qc_1.as_bytes());

    for value in public_input.iter() {
        transcript.append_field_element(value);
    }

    transcript
}

#[cfg(test)]
mod tests {
    use lambdaworks_math::elliptic_curve::short_weierstrass::curves::bls12_381::default_types::FrField;
    use lambdaworks_math::elliptic_curve::{
        short_weierstrass::curves::bls12_381::curve::BLS12381Curve, traits::IsEllipticCurve,
    };

    use super::*;
    use crate::test_utils::circuit_1::test_common_preprocessed_input_1;
    use crate::test_utils::utils::{test_srs, FpElement, KZG};

    #[test]
    fn setup_works_for_simple_circuit() {
        let common_input = test_common_preprocessed_input_1();
        let srs = test_srs(common_input.n);
        let kzg = KZG::new(srs);

        let vk = setup::<FrField, KZG>(&common_input, &kzg);

        let expected_ql = BLS12381Curve::create_point_from_affine(
            FpElement::from_hex_unchecked("1492341357755e31a6306abf3237f84f707ded7cb526b8ffd40901746234ef27f12bc91ef638e4977563db208b765f12"),
            FpElement::from_hex_unchecked("ec3ff8288ea339010658334f494a614f7470c19a08d53a9cf5718e0613bb65d2cdbc1df374057d9b45c35cf1f1b5b72"),
        ).unwrap();
        let expected_qr = BLS12381Curve::create_point_from_affine(
            FpElement::from_hex_unchecked("107ab09b6b8c6fc55087aeb8045e17a6d016bdacbc64476264328e71f3e85a4eacaee34ee963e9c9249b6b1bc9653674"),
            FpElement::from_hex_unchecked("f98e3fe5a53545b67a51da7e7a6cedc51af467abdefd644113fb97edf339aeaa5e2f6a5713725ec76754510b76a10be"),
        ).unwrap();
        let expected_qo = BLS12381Curve::create_point_from_affine(
            FpElement::from_hex_unchecked("9fd00baa112a0064ce5c3c2d243e657b25df8a2f237b91eec27e83157f6ca896a2401d07ec7d7d097d2f2a344e2018f"),
            FpElement::from_hex_unchecked("15922cfa65972d80823c6bb9aeb0637c864b636267bfee2818413e9cdc5f7948575c4ce097bb8b9db8087c4ed5056592"),
        ).unwrap();
        let expected_qm = BLS12381Curve::create_point_from_affine(
            FpElement::from_hex_unchecked("9fd00baa112a0064ce5c3c2d243e657b25df8a2f237b91eec27e83157f6ca896a2401d07ec7d7d097d2f2a344e2018f"),
            FpElement::from_hex_unchecked("46ee4efd3e8b919c8df3bfc949b495ade2be8228bc524974eef94041a517cdbc74fb31e1998746201f683b12afa4519"),
        ).unwrap();

        let expected_s1 = BLS12381Curve::create_point_from_affine(
            FpElement::from_hex_unchecked("187ee12de08728650d18912aa9fe54863922a9eeb37e34ff43041f1d039f00028ad2cdf23705e6f6ab7ea9406535c1b0"),
            FpElement::from_hex_unchecked("4f29051990de0d12b38493992845d9abcb48ef18239eca8b8228618c78ec371d39917bc0d45cf6dc4f79bd64baa9ae2")
        ).unwrap();
        let expected_s2 = BLS12381Curve::create_point_from_affine(
            FpElement::from_hex_unchecked("167c0384025887c01ea704234e813842a4acef7d765c3a94a5442ca685b4fc1d1b425ba7786a7413bd4a7d6a1eb5a35a"),
            FpElement::from_hex_unchecked("12b644100c5d00af27c121806c4779f88e840ff3fdac44124b8175a303d586c4d910486f909b37dda1505c485f053da1")
        ).unwrap();
        let expected_s3 = BLS12381Curve::create_point_from_affine(
            FpElement::from_hex_unchecked("188fb6dba3cf5af8a7f6a44d935bb3dd2083a5beb4c020f581739ebc40659c824a4ca8279cf7d852decfbca572e4fa0e"),
            FpElement::from_hex_unchecked("d84d52582fd95bfa7672f7cef9dd4d0b1b4a54d33f244fdb97df71c7d45fd5c5329296b633c9ed23b8475ee47b9d99")
        ).unwrap();

        assert_eq!(vk.ql_1, expected_ql);
        assert_eq!(vk.qr_1, expected_qr);
        assert_eq!(vk.qo_1, expected_qo);
        assert_eq!(vk.qm_1, expected_qm);

        assert_eq!(vk.s1_1, expected_s1);
        assert_eq!(vk.s2_1, expected_s2);
        assert_eq!(vk.s3_1, expected_s3);
    }
}
