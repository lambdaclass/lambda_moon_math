use crate::config::FriMerkleTree;

use super::config::Commitment;
use lambdaworks_math::{
    field::{element::FieldElement, fields::mersenne31::field::Mersenne31Field}, 
    circle::polynomial::{interpolate_cfft, evaluate_cfft}
};

const BLOW_UP_FACTOR: usize = 2;

pub fn prove(trace: Vec<FieldElement<Mersenne31Field>>) -> Commitment {
    
    let lde_domain_size = trace.len() * BLOW_UP_FACTOR;

    // Returns the coef of the interpolating polinomial of the trace on a natural domain.    
    let mut trace_poly = interpolate_cfft(trace);

    // Padding with zeros the coefficients of the polynomial, so we can evaluate it in the lde domain.
    trace_poly.resize(lde_domain_size, FieldElement::zero());
    let lde_eval = evaluate_cfft(trace_poly);

    let tree = FriMerkleTree::<Mersenne31Field>::build(&lde_eval).unwrap();
    let commitment = tree.root;

    commitment
}


#[cfg(test)]
mod tests {
    
    use super::*;

    type FE = FieldElement<Mersenne31Field>;

    #[test]
    fn basic_test() {
        let trace = vec![
            FE::from(1),
            FE::from(2),
            FE::from(3),
            FE::from(4),
            FE::from(5),
            FE::from(6),
            FE::from(7),
            FE::from(8),
        ];

        let commitmet = prove(trace);
        println!("{:?}", commitmet); 
    }
}