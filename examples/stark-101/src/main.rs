use lambdaworks_math::field::element::FieldElement;
use stark_101::field::U64Stark101PrimeField;

type Fe = FieldElement<U64Stark101PrimeField>;

fn fibonacci_square_trace(a_0: Fe, a_1: Fe, n: usize) -> Vec<Fe> {
    let mut a = vec![a_0, a_1];
    while a.len() < n {
        a.push(a[a.len() - 1].square() + a[a.len() - 1].square());
    }
    a
}

fn main() {
    println!("hello");
}

#[cfg(test)]
mod test_stark_101 {
    use super::*;

    #[test]
    fn test_fibonacci_square_trace() {
        let trace = fibonacci_square_trace(Fe::one(), Fe::from(5u64), 3);
        println!("{:?}", trace);
    }
}
