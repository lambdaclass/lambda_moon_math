use lambdaworks_math::{
    elliptic_curve::{
        short_weierstrass::curves::bn_254::{
            curve::BN254Curve, default_types::FrElement as FE, default_types::FrField as FrF,
            pairing::BN254AtePairing, twist::BN254TwistCurve,
        },
        traits::{IsEllipticCurve, IsPairing},
    },
    field::element::FieldElement,
    unsigned_integer::element::U256,
};
use rand::{Rng, SeedableRng};

pub type Curve = BN254Curve;
pub type TwistedCurve = BN254TwistCurve;

pub type FrElement = FE;
pub type FrField = FrF;

pub type Pairing = BN254AtePairing;

pub type G1Point = <BN254Curve as IsEllipticCurve>::PointRepresentation;
pub type G2Point = <BN254TwistCurve as IsEllipticCurve>::PointRepresentation;
pub type PairingOutput = FieldElement<<Pairing as IsPairing>::OutputField>;

/// g = ORDER_R_MINUS_1_ROOT_UNITY is a primitive unity root of order r - 1.
/// I.e. g^{r-1} mod r = 1 and g^i != 1 for i < r-1.
/// We calculated g in Sage:
///
/// r = 21888242871839275222246405745257275088548364400416034343698204186575808495617
/// F = GF(r)
/// g = F.primitive_element()
/// print("Primitive generator:", g)
pub const ORDER_R_MINUS_1_ROOT_UNITY: FrElement = FrElement::from_hex_unchecked("5");

pub fn sample_fr_elem() -> FrElement {
    let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(9001);
    FrElement::new(U256 {
        limbs: [
            rng.gen::<u64>(),
            rng.gen::<u64>(),
            rng.gen::<u64>(),
            rng.gen::<u64>(),
        ],
    })
}
