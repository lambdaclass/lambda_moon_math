use lambdaworks_math::{
    elliptic_curve::{
        short_weierstrass::curves::bls12_381::{
            curve::BLS12381Curve, default_types::FrElement as FE, default_types::FrField as FrF,
            pairing::BLS12381AtePairing, twist::BLS12381TwistCurve,
        },
        traits::{IsEllipticCurve, IsPairing},
    },
    field::element::FieldElement,
    unsigned_integer::element::U256,
};
use rand::{Rng, SeedableRng};

pub type Curve = BLS12381Curve;
pub type TwistedCurve = BLS12381TwistCurve;

pub type FrElement = FE;
pub type FrField = FrF;

pub type Pairing = BLS12381AtePairing;

pub type G1Point = <BLS12381Curve as IsEllipticCurve>::PointRepresentation;
pub type G2Point = <BLS12381TwistCurve as IsEllipticCurve>::PointRepresentation;
pub type PairingOutput = FieldElement<<Pairing as IsPairing>::OutputField>;

pub const ORDER_R_MINUS_1_ROOT_UNITY: FrElement = FrElement::from_hex_unchecked("7");

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
