use icicle_bls12_377::curve;
use icicle_core::{error::IcicleError, msm, traits::FieldImpl};
use icicle_cuda_runtime::{memory::HostOrDeviceSlice, stream::CudaStream};

use crate::{
    elliptic_curve::short_weierstrass::{
        curves::bls12_377::{
            curve::{BLS12377Curve, BLS12377FieldElement},
            field_extension::BLS12377PrimeField,
        },
        point::ShortWeierstrassProjectivePoint,
    },
    errors::ByteConversionError,
    field::element::FieldElement,
    traits::ByteConversion,
};

impl BLS12377FieldElement {
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

impl ShortWeierstrassProjectivePoint<BLS12377Curve> {
    fn to_icicle(&self) -> curve::G1Affine {
        let s = self.to_affine();
        curve::G1Affine {
            x: s.x().to_icicle(),
            y: s.y().to_icicle(),
        }
    }

    fn from_icicle(icicle: &curve::G1Projective) -> Result<Self, ByteConversionError> {
        Ok(Self::new([
            FieldElement::<BLS12377PrimeField>::from_icicle(&icicle.x).unwrap(),
            FieldElement::<BLS12377PrimeField>::from_icicle(&icicle.y).unwrap(),
            FieldElement::<BLS12377PrimeField>::from_icicle(&icicle.z).unwrap(),
        ]))
    }
}

pub fn bls12_377_g1_msm(
    scalars: &[BLS12377FieldElement],
    points: &[ShortWeierstrassProjectivePoint<BLS12377Curve>],
    config: Option<msm::MSMConfig>,
) -> Result<ShortWeierstrassProjectivePoint<BLS12377Curve>, IcicleError> {
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
        ShortWeierstrassProjectivePoint::<BLS12377Curve>::from_icicle(&msm_host_result[0]).unwrap();
    Ok(res)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        elliptic_curve::{
            short_weierstrass::curves::bls12_377::curve::BLS12377FieldElement,
            traits::IsEllipticCurve,
        },
        field::element::FieldElement,
        msm::pippenger::msm,
    };

    impl ShortWeierstrassProjectivePoint<BLS12377Curve> {
        fn from_icicle_affine(
            icicle: &curve::G1Affine,
        ) -> Result<ShortWeierstrassProjectivePoint<BLS12377Curve>, ByteConversionError> {
            Ok(Self::new([
                FieldElement::<BLS12377PrimeField>::from_icicle(&icicle.x).unwrap(),
                FieldElement::<BLS12377PrimeField>::from_icicle(&icicle.y).unwrap(),
                FieldElement::one(),
            ]))
        }
    }

    fn point_times_5() -> ShortWeierstrassProjectivePoint<BLS12377Curve> {
        let x = BLS12377FieldElement::from_hex_unchecked(
            "3c852d5aab73fbb51e57fbf5a0a8b5d6513ec922b2611b7547bfed74cba0dcdfc3ad2eac2733a4f55d198ec82b9964",
        );
        let y = BLS12377FieldElement::from_hex_unchecked(
            "a71425e68e55299c64d7eada9ae9c3fb87a9626b941d17128b64685fc07d0e635f3c3a512903b4e0a43e464045967b",
        );
        BLS12377Curve::create_point_from_affine(x, y).unwrap()
    }

    #[test]
    fn to_from_icicle() {
        // convert value of 5 to icicle and back again and that icicle 5 matches
        let point = point_times_5();
        let icicle_point = point.to_icicle();
        let res =
            ShortWeierstrassProjectivePoint::<BLS12377Curve>::from_icicle_affine(&icicle_point)
                .unwrap();
        assert_eq!(point, res)
    }

    #[test]
    fn to_from_icicle_generator() {
        // Convert generator and see that it matches
        let point = BLS12377Curve::generator();
        let icicle_point = point.to_icicle();
        let res =
            ShortWeierstrassProjectivePoint::<BLS12377Curve>::from_icicle_affine(&icicle_point)
                .unwrap();
        assert_eq!(point, res)
    }

    #[test]
    fn icicle_g1_msm() {
        const LEN: usize = 20;
        let eight: BLS12377FieldElement = FieldElement::from(8);
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
        let res = bls12_377_g1_msm(&lambda_scalars, &lambda_points, None).unwrap();
        assert_eq!(res, expected);
    }
}
