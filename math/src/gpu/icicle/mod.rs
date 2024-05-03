use crate::{
    cyclic_group::IsGroup,
    elliptic_curve::{
        short_weierstrass::{
            curves::{
                bls12_377::curve::BLS12377Curve,
                bls12_381::{curve::BLS12381Curve, twist::BLS12381TwistCurve, field_extension::BLS12381PrimeField},
                bn_254::{curve::BN254Curve, field_extension::BN254PrimeField}
            },
            point::ShortWeierstrassProjectivePoint,
        },
        traits::IsEllipticCurve,
    },
    errors::ByteConversionError,
    field::{element::FieldElement, traits::{IsField, IsSubFieldOf, IsFFTField}},
    msm::naive::MSMError,
    fft::errors::FFTError,
    traits::ByteConversion,
};
use icicle_bls12_377::curve::CurveCfg as IcicleBLS12377Curve;
use icicle_bls12_381::curve::{
    CurveCfg as IcicleBLS12381Curve,
    ScalarCfg as IcicleBLS12381ScalarCfg
};
use icicle_bn254::curve::{
    CurveCfg as IcicleBN254Curve,
    ScalarCfg as IcicleBN254ScalarCfg
};
use icicle_core::{
    curve::{Affine, Curve, Projective},
    msm,
    ntt::{ntt, NTT, NTTConfig, NTTDir},
    traits::FieldImpl,
};
use icicle_cuda_runtime::{memory::HostOrDeviceSlice, stream::CudaStream};

use std::fmt::Debug;

impl GpuMSMPoint for ShortWeierstrassProjectivePoint<BLS12381Curve> {
    type LambdaCurve = BLS12381Curve;
    type GpuCurve = IcicleBLS12381Curve;

    fn curve_name() -> &'static str {
        "BLS12381"
    }

    fn to_icicle_affine(point: &Self) -> Affine<Self::GpuCurve> {
        let s = Self::to_affine(point);
        Affine::<Self::GpuCurve> {
            x: Self::to_icicle_field(s.x()),
            y: Self::to_icicle_field(s.y()),
        }
    }

    fn from_icicle_projective(
        icicle: &Projective<Self::GpuCurve>,
    ) -> Result<Self, ByteConversionError> {
        Ok(Self::new([
            Self::from_icicle_field(&icicle.x).unwrap(),
            Self::from_icicle_field(&icicle.y).unwrap(),
            Self::from_icicle_field(&icicle.z).unwrap(),
        ]))
    }
}

//NOTE THIS IS A PLACEHOLDER COMPILING ICICLE G2 TOOK TO LONG!
impl GpuMSMPoint for ShortWeierstrassProjectivePoint<BLS12381TwistCurve> {
    type LambdaCurve = BLS12381Curve;
    type GpuCurve = IcicleBLS12381Curve;

    fn curve_name() -> &'static str {
        ""
    }

    fn to_icicle_affine(point: &Self) -> Affine<Self::GpuCurve> {
        let s = Self::to_affine(point);
        Affine::<Self::GpuCurve> {
            x: Self::to_icicle_field(s.x()),
            y: Self::to_icicle_field(s.y()),
        }
    }

    fn from_icicle_projective(
        icicle: &Projective<Self::GpuCurve>,
    ) -> Result<Self, ByteConversionError> {
        Ok(Self::new([
            Self::from_icicle_field(&icicle.x).unwrap(),
            Self::from_icicle_field(&icicle.y).unwrap(),
            Self::from_icicle_field(&icicle.z).unwrap(),
        ]))
    }
}

impl GpuMSMPoint for ShortWeierstrassProjectivePoint<BLS12377Curve> {
    type LambdaCurve = BLS12377Curve;
    type GpuCurve = IcicleBLS12377Curve;
    fn curve_name() -> &'static str {
        "BLS12377"
    }

    fn to_icicle_affine(point: &Self) -> Affine<Self::GpuCurve> {
        let s = Self::to_affine(point);
        Affine::<Self::GpuCurve> {
            x: Self::to_icicle_field(s.x()),
            y: Self::to_icicle_field(s.y()),
        }
    }

    fn from_icicle_projective(
        icicle: &Projective<Self::GpuCurve>,
    ) -> Result<Self, ByteConversionError> {
        Ok(Self::new([
            Self::from_icicle_field(&icicle.x).unwrap(),
            Self::from_icicle_field(&icicle.y).unwrap(),
            Self::from_icicle_field(&icicle.z).unwrap(),
        ]))
    }
}

impl GpuMSMPoint for ShortWeierstrassProjectivePoint<BN254Curve> {
    type LambdaCurve = BN254Curve;
    type GpuCurve = IcicleBN254Curve;
    fn curve_name() -> &'static str {
        "BN254"
    }

    fn to_icicle_affine(point: &Self) -> Affine<Self::GpuCurve> {
        let s = Self::to_affine(point);
        Affine::<Self::GpuCurve> {
            x: Self::to_icicle_field(s.x()),
            y: Self::to_icicle_field(s.y()),
        }
    }

    fn from_icicle_projective(
        icicle: &Projective<Self::GpuCurve>,
    ) -> Result<Self, ByteConversionError> {
        Ok(Self::new([
            Self::from_icicle_field(&icicle.x).unwrap(),
            Self::from_icicle_field(&icicle.y).unwrap(),
            Self::from_icicle_field(&icicle.z).unwrap(),
        ]))
    }
}

pub trait GpuMSMPoint: IsGroup {
    type LambdaCurve: IsEllipticCurve + Clone + Debug;
    type GpuCurve: Curve + msm::MSM<Self::GpuCurve>;
    //type FE: ByteConversion;
    /// Used for searching this field's implementation in other languages, e.g in MSL
    /// for executing parallel operations with the Metal API.
    fn curve_name() -> &'static str {
        ""
    }

    fn to_icicle_affine(point: &Self) -> Affine<Self::GpuCurve>;

    fn from_icicle_projective(
        icicle: &Projective<Self::GpuCurve>,
    ) -> Result<Self, ByteConversionError>;

    fn to_icicle_field<FE: ByteConversion>(element: &FE) -> <Self::GpuCurve as Curve>::BaseField {
        <Self::GpuCurve as Curve>::BaseField::from_bytes_le(&element.to_bytes_le())
    }

    fn to_icicle_scalar<FE: ByteConversion>(
        element: &FE,
    ) -> <Self::GpuCurve as Curve>::ScalarField {
        <Self::GpuCurve as Curve>::ScalarField::from_bytes_le(&element.to_bytes_le())
    }

    fn from_icicle_field<FE: ByteConversion>(
        icicle: &<Self::GpuCurve as Curve>::BaseField,
    ) -> Result<FE, ByteConversionError> {
        FE::from_bytes_le(&icicle.to_bytes_le())
    }
}

pub trait IcicleFFT: IsField 
where
    FieldElement<Self>: ByteConversion,
    <Self::ScalarField as FieldImpl>::Config: NTT<Self::ScalarField>
{
    type ScalarField: FieldImpl;

    fn to_icicle_scalar(element: &FieldElement<Self>) -> Self::ScalarField {
        Self::ScalarField::from_bytes_le(&element.to_bytes_le())
    }

    fn from_icicle_scalar(
        icicle: &Self::ScalarField,
    ) -> Result<FieldElement<Self>, ByteConversionError> {
        FieldElement::<Self>::from_bytes_le(&icicle.to_bytes_le())
    }
}

impl IcicleFFT for BLS12381PrimeField {
    type ScalarField = <IcicleBLS12381Curve as Curve>::ScalarField;
}

impl IcicleFFT for BN254PrimeField {
    type ScalarField = <IcicleBN254Curve as Curve>::ScalarField;
}

pub fn icicle_msm<F: IsField, G: GpuMSMPoint>(
    cs: &[FieldElement<F>],
    points: &[G],
) -> Result<G, MSMError>
where
    FieldElement<F>: ByteConversion,
{
    let mut cfg = msm::MSMConfig::default();
    let scalars = HostOrDeviceSlice::Host(
        cs.iter()
            .map(|scalar| G::to_icicle_scalar(scalar))
            .collect::<Vec<_>>(),
    );
    let points = HostOrDeviceSlice::Host(
        points
            .iter()
            .map(|point| G::to_icicle_affine(point))
            .collect::<Vec<_>>(),
    );
    let mut msm_results = HostOrDeviceSlice::cuda_malloc(1).unwrap();
    let stream = CudaStream::create().unwrap();
    cfg.ctx.stream = &stream;
    cfg.is_async = true;
    msm::msm(&scalars, &points, &cfg, &mut msm_results).unwrap();
    let mut msm_host_result = [Projective::<G::GpuCurve>::zero(); 1];
    stream.synchronize().unwrap();
    msm_results.copy_to_host(&mut msm_host_result[..]).unwrap();
    stream.destroy().unwrap();
    let res = G::from_icicle_projective(&msm_host_result[0]).unwrap();
    Ok(res)
}

pub fn evaluate_fft_icicle<F, E>(
    coeffs: &Vec<FieldElement<E>>,
) -> Result<Vec<FieldElement<E>>, FFTError> 
where
    E: IsSubFieldOf<E>,
    FieldElement<E>: ByteConversion,
    <<E as IcicleFFT>::ScalarField as FieldImpl>::Config: NTT<<E as IcicleFFT>::ScalarField>,
    E: IsField + IcicleFFT + IsFFTField,
{
    let size = coeffs.len();
    let mut cfg = NTTConfig::default();
    let order = coeffs.len() as u64;
    let dir = NTTDir::kForward;
    let scalars = HostOrDeviceSlice::Host(
        coeffs
            .iter()
            .map(|scalar| E::to_icicle_scalar(&scalar))
            .collect::<Vec<_>>(),
    );
    let mut ntt_results = HostOrDeviceSlice::cuda_malloc(size).unwrap();
    let stream = CudaStream::create().unwrap();
    cfg.ctx.stream = &stream;
    cfg.is_async = true;
    let root_of_unity = E::to_icicle_scalar(&E::get_primitive_root_of_unity(order).unwrap());
    <E::ScalarField as FieldImpl>::Config::initialize_domain(root_of_unity, &cfg.ctx).unwrap();
    ntt(&scalars, dir, &cfg, &mut ntt_results).unwrap();
    stream.synchronize().unwrap();
    let mut ntt_host_results = vec![E::ScalarField::zero(); size];
    ntt_results.copy_to_host(&mut ntt_host_results[..]).unwrap();
    stream.destroy().unwrap();
    let res = ntt_host_results
        .as_slice()
        .iter()
        .map(|scalar| E::from_icicle_scalar(&scalar).unwrap())
        .collect::<Vec<_>>();
    Ok(res)
}

pub fn interpolate_fft_icicle<F, E>(
    coeffs: &Vec<FieldElement<E>>,
) -> Result<Vec<FieldElement<E>>, FFTError> 
where
    F: IsSubFieldOf<E>,
    FieldElement<E>: ByteConversion,
    <<E as IcicleFFT>::ScalarField as FieldImpl>::Config: NTT<<E as IcicleFFT>::ScalarField>,
    E: IsField + IcicleFFT + IsFFTField,
{
    let size = coeffs.len();
    let mut cfg = NTTConfig::default();
    let order = coeffs.len() as u64;
    let dir = NTTDir::kInverse;
    let scalars = HostOrDeviceSlice::Host(
        coeffs
            .iter()
            .map(|scalar| E::to_icicle_scalar(scalar))
            .collect::<Vec<_>>(),
    );
    let mut ntt_results = HostOrDeviceSlice::cuda_malloc(size).unwrap();
    let stream = CudaStream::create().unwrap();
    cfg.ctx.stream = &stream;
    cfg.is_async = true;
    let root_of_unity = E::to_icicle_scalar(&E::get_primitive_root_of_unity(order).unwrap());
    <E::ScalarField as FieldImpl>::Config::initialize_domain(root_of_unity, &cfg.ctx).unwrap();
    ntt(&scalars, dir, &cfg, &mut ntt_results).unwrap();
    stream.synchronize().unwrap();
    let mut ntt_host_results = vec![E::ScalarField::zero(); size];
    ntt_results.copy_to_host(&mut ntt_host_results[..]).unwrap();
    stream.destroy().unwrap();
    let res = ntt_host_results
        .as_slice()
        .iter()
        .map(|scalar| E::from_icicle_scalar(&scalar).unwrap())
        .collect::<Vec<_>>();
    Ok(res)
}


#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        elliptic_curve::{
            short_weierstrass::curves::bls12_381::curve::BLS12381FieldElement,
            traits::IsEllipticCurve,
        },
        field::element::FieldElement,
        msm::pippenger::msm,
    };

    impl ShortWeierstrassProjectivePoint<BLS12381Curve> {
        fn from_icicle_affine(
            icicle: &curve::G1Affine,
        ) -> Result<ShortWeierstrassProjectivePoint<BLS12381Curve>, ByteConversionError> {
            Ok(Self::new([
                FieldElement::<BLS12381PrimeField>::from_icicle(&icicle.x).unwrap(),
                FieldElement::<BLS12381PrimeField>::from_icicle(&icicle.y).unwrap(),
                FieldElement::one(),
            ]))
        }
    }

    fn point_times_5() -> ShortWeierstrassProjectivePoint<BLS12381Curve> {
        let x = BLS12381FieldElement::from_hex_unchecked(
            "32bcce7e71eb50384918e0c9809f73bde357027c6bf15092dd849aa0eac274d43af4c68a65fb2cda381734af5eecd5c",
        );
        let y = BLS12381FieldElement::from_hex_unchecked(
            "11e48467b19458aabe7c8a42dc4b67d7390fdf1e150534caadddc7e6f729d8890b68a5ea6885a21b555186452b954d88",
        );
        BLS12381Curve::create_point_from_affine(x, y).unwrap()
    }

    #[test]
    fn to_from_icicle() {
        // convert value of 5 to icicle and back again and that icicle 5 matches
        let point = point_times_5();
        let icicle_point = point.to_icicle();
        let res =
            ShortWeierstrassProjectivePoint::<BLS12381Curve>::from_icicle_affine(&icicle_point)
                .unwrap();
        assert_eq!(point, res)
    }

    #[test]
    fn to_from_icicle_generator() {
        // Convert generator and see that it matches
        let point = BLS12381Curve::generator();
        let icicle_point = point.to_icicle();
        let res =
            ShortWeierstrassProjectivePoint::<BLS12381Curve>::from_icicle_affine(&icicle_point)
                .unwrap();
        assert_eq!(point, res)
    }

    #[test]
    fn icicle_g1_msm() {
        const LEN: usize = 20;
        let eight: BLS12381FieldElement = FieldElement::from(8);
        let lambda_scalars = vec![eight; LEN];
        let lambda_points = (0..LEN).map(|_| point_times_5()).collect::<Vec<_>>();
        let expected = msm(&lambda_scalars, &lambda_points).unwrap();
        let res = bls12_381_g1_msm(&lambda_scalars, &lambda_points, None).unwrap();
        assert_eq!(res, expected);
    }
}
