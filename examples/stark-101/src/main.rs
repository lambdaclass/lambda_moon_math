use lambdaworks_crypto::{
    fiat_shamir::{default_transcript::DefaultTranscript, is_transcript::IsTranscript},
    merkle_tree::{backends::types::Sha2_256Backend, merkle::MerkleTree},
};
use lambdaworks_math::{field::element::FieldElement, polynomial::Polynomial};
use stark_101::field::U64Stark101PrimeField;

type Fe = FieldElement<U64Stark101PrimeField>;
type Merkle = MerkleTree<Sha2_256Backend<U64Stark101PrimeField>>;

const X: u64 = 3141592;
const N: u64 = 1024;

fn fibonacci_square_trace(a_0: Fe, a_1: Fe, n: usize) -> Vec<Fe> {
    let mut a = vec![a_0, a_1];
    while a.len() < n {
        a.push(a[a.len() - 1].square() + a[a.len() - 2].square());
    }
    a
}

fn get_interpolating_polynomial(a: &Vec<Fe>, G: Vec<Fe>) -> Polynomial<Fe> {
    Polynomial::interpolate(&G.as_slice()[..G.len() - 1], a).expect("err")
}

fn get_group_of_order(generator: Fe, order: u64) -> Vec<Fe> {
    let k = 3221225472 / order;
    let new_generaor = generator.pow(k);
    (0..order)
        .into_iter()
        .map(|i| new_generaor.pow(i as u64))
        .collect()
}

fn main() {
    // We know that 5 is a generator of the multiplicative group of F by brute force.
    let G: Vec<Fe> = get_group_of_order(Fe::from(5), N);
    let a = fibonacci_square_trace(Fe::one(), Fe::from(X), 1023);

    // let f = get_interpolating_polynomial(&a, G);

    // let H: Vec<Fe> = get_group_of_order(Fe::from(5), N * 8);
    // let eval_domain: Vec<Fe> = H.into_iter().map(|h| Fe::from(5) * h).collect();

    // let f_eval = eval_domain
    //     .clone()
    //     .into_iter()
    //     .map(|x| f.evaluate(&x))
    //     .collect::<Vec<_>>();

    // let f_merkle = Merkle::build(&f_eval).unwrap();

    // assert_eq!(
    //     f_merkle.root.to_vec(),
    //     hex::decode("6c266a104eeaceae93c14ad799ce595ec8c2764359d7ad1b4b7c57a4da52be04").unwrap(),
    // );

    let mkt = Merkle::build(&vec![Fe::from(0u64), Fe::from(0u64)]).unwrap();
    // println!("{:?}", Fe::from(30u64));
    println!("{:?}", hex::encode(mkt.root));
}

#[cfg(test)]
mod test_stark_101 {
    use super::*;

    #[test]
    fn test_fibonacci_square_trace() {
        let trace = fibonacci_square_trace(Fe::one(), Fe::from(3141592u64), 1023);
        assert_eq!(trace[1022], Fe::from(2338775057u64));
    }

    #[test]
    fn g_has_order_1024() {
        let g = Fe::from(5).pow(3 * 2u64.pow(20));
        assert_eq!(g.pow(1024u64), Fe::one())
    }

    #[test]
    fn f_interpolates() {
        let g = Fe::from(5).pow(3 * 2u64.pow(20));
        let G: Vec<Fe> = (0..N).into_iter().map(|i| g.pow(i as u64)).collect();
        let a = fibonacci_square_trace(Fe::one(), Fe::from(X), 1023);
        let f = get_interpolating_polynomial(&a, G);
        assert_eq!(f.evaluate(&g), Fe::from(X))
    }
}
