use icicle_bn254::curve;
use icicle_core::{error::IcicleError, msm, traits::FieldImpl};
use icicle_cuda_runtime::{memory::HostOrDeviceSlice, stream::CudaStream};

use crate::{
    elliptic_curve::short_weierstrass::{
        curves::bn_254::{
            curve::{BN254Curve, BN254FieldElement},
            field_extension::BN254PrimeField,
        },
        point::ShortWeierstrassProjectivePoint,
    },
    errors::ByteConversionError,
    field::element::FieldElement,
    traits::ByteConversion,
};

impl BN254FieldElement {
    fn to_icicle(&self) -> curve::BaseField {
        curve::BaseField::from_bytes_le(&self.to_bytes_le())
    }

    fn to_icicle_scalar(&self) -> curve::ScalarField {
        curve::ScalarField::from_bytes_le(&self.to_bytes_le())
    }

    fn from_icicle(icicle: &curve::BaseField) -> Result<Self, ByteConversionError> {
        Self::from_bytes_le(&icicle.to_bytes_le())
    }
}

impl ShortWeierstrassProjectivePoint<BN254Curve> {
    fn to_icicle(&self) -> curve::G1Affine {
        let s = self.to_affine();
        curve::G1Affine {
            x: s.x().to_icicle(),
            y: s.y().to_icicle(),
        }
    }

    fn from_icicle(icicle: &curve::G1Projective) -> Result<Self, ByteConversionError> {
        Ok(Self::new([
            FieldElement::<BN254PrimeField>::from_icicle(&icicle.x).unwrap(),
            FieldElement::<BN254PrimeField>::from_icicle(&icicle.y).unwrap(),
            FieldElement::<BN254PrimeField>::from_icicle(&icicle.z).unwrap(),
        ]))
    }
}

pub fn bn254_g1_msm(
    scalars: &[BN254FieldElement],
    points: &[ShortWeierstrassProjectivePoint<BN254Curve>],
    config: Option<msm::MSMConfig>,
) -> Result<ShortWeierstrassProjectivePoint<BN254Curve>, IcicleError> {
    let mut cfg = config.unwrap_or(msm::MSMConfig::default());
    let scalars = HostOrDeviceSlice::Host(
        scalars
            .iter()
            .map(|scalar| scalar.to_icicle_scalar())
            .collect::<Vec<_>>(),
    );
    let points = HostOrDeviceSlice::Host(
        points
            .iter()
            .map(|point| point.to_icicle())
            .collect::<Vec<_>>(),
    );
    let mut msm_results = HostOrDeviceSlice::cuda_malloc(1).unwrap();
    let stream = CudaStream::create().unwrap();
    cfg.ctx.stream = &stream;
    cfg.is_async = true;
    msm::msm(&scalars, &points, &cfg, &mut msm_results).unwrap();
    let mut msm_host_result = vec![curve::G1Projective::zero(); 1];
    stream.synchronize().unwrap();
    msm_results.copy_to_host(&mut msm_host_result[..]).unwrap();
    stream.destroy().unwrap();
    let res =
        ShortWeierstrassProjectivePoint::<BN254Curve>::from_icicle(&msm_host_result[0]).unwrap();
    Ok(res)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        elliptic_curve::{
            short_weierstrass::curves::bn_254::curve::BN254FieldElement, traits::IsEllipticCurve,
        },
        field::element::FieldElement,
        msm::pippenger::msm,
    };

    impl ShortWeierstrassProjectivePoint<BN254Curve> {
        fn from_icicle_affine(
            icicle: &curve::G1Affine,
        ) -> Result<ShortWeierstrassProjectivePoint<BN254Curve>, ByteConversionError> {
            Ok(Self::new([
                FieldElement::<BN254PrimeField>::from_icicle(&icicle.x).unwrap(),
                FieldElement::<BN254PrimeField>::from_icicle(&icicle.y).unwrap(),
                FieldElement::one(),
            ]))
        }
    }

    fn point_times_5() -> ShortWeierstrassProjectivePoint<BN254Curve> {
        let x = BN254FieldElement::from_hex_unchecked(
            "16ab03b69dfb4f870b0143ebf6a71b7b2e4053ca7a4421d09a913b8b834bbfa3",
        );
        let y = BN254FieldElement::from_hex_unchecked(
            "2512347279ba1049ef97d4ec348d838f939d2b7623e88f4826643cf3889599b2",
        );
        BN254Curve::create_point_from_affine(x, y).unwrap()
    }

    #[test]
    fn to_from_icicle() {
        // convert value of 5 to icicle and back again and that icicle 5 matches
        let point = point_times_5();
        let icicle_point = point.to_icicle();
        let res = ShortWeierstrassProjectivePoint::<BN254Curve>::from_icicle_affine(&icicle_point)
            .unwrap();
        assert_eq!(point, res)
    }

    #[test]
    fn to_from_icicle_generator() {
        // Convert generator and see that it matches
        let point = BN254Curve::generator();
        let icicle_point = point.to_icicle();
        let res = ShortWeierstrassProjectivePoint::<BN254Curve>::from_icicle_affine(&icicle_point)
            .unwrap();
        assert_eq!(point, res)
    }

    #[test]
    fn icicle_g1_msm() {
        const LEN: usize = 20;
        let eight: BN254FieldElement = FieldElement::from(8);
        let lambda_scalars = vec![eight; LEN];
        let lambda_points = (0..LEN).map(|_| point_times_5()).collect::<Vec<_>>();
        let expected = msm(
            &lambda_scalars
                .clone()
                .into_iter()
                .map(|x| x.representative())
                .collect::<Vec<_>>(),
            &lambda_points,
        )
        .unwrap();
        let res = bn254_g1_msm(&lambda_scalars, &lambda_points, None).unwrap();
        assert_eq!(res, expected);
    }
}
