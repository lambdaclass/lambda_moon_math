use core::fmt::Display;

use crate::{
    field::{
        element::FieldElement,
        traits::{IsField, IsPrimeField},
    },
    polynomial::traits::term::Term,
};
use std::ops;

/// Struct for (coeff: FieldElement<F>, terms: Vec<usize>) representing a multilinear
/// monomial in a sparse format.
// TODO: add check that var labels are 0 indexed
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultilinearMonomial<F: IsField + IsPrimeField>
where
    <F as IsField>::BaseType: Send + Sync,
{
    pub coeff: FieldElement<F>,
    pub vars: Vec<usize>,
}

impl<F: IsField + IsPrimeField> MultilinearMonomial<F>
where
    <F as IsField>::BaseType: Send + Sync,
{
    /// Create a new `Term` from a tuple of the form `(coeff, (variables))`
    pub fn new(term: (FieldElement<F>, Vec<usize>)) -> Self {
        // sort variables in increasing order
        let mut vars = term.1;
        vars.sort();

        Self {
            coeff: term.0,
            vars,
        }
    }
}

impl<F: IsField + IsPrimeField> Term<F> for MultilinearMonomial<F>
where
    <F as IsField>::BaseType: Send + Sync,
{
    /// Returns the total degree of `self`. This is the count of all variables
    fn degree(&self) -> usize {
        self.vars.len()
    }

    /// Returns a list of the powers of each variable in `self` i.e. numbers representing the id of
    /// the specific variable
    fn vars(&self) -> Vec<usize> {
        self.vars.clone()
    }

    fn powers(&self) -> Vec<usize> {
        vec![1; self.vars.len()]
    }

    /// Fetches the max variable by id from the sparse list of id's this is used to ensure the upon
    /// evaluation the correct number of points are supplied
    fn max_var(&self) -> usize {
        // Sparse variables are stored in increasing var_id therefore we grab the last one
        match self.vars.last() {
            Some(&max_var) => max_var,
            None => 0,
        }
    }

    // TODO: test this
    /// Evaluates `self` at the point `p`.
    fn evaluate(&self, p: &[FieldElement<F>]) -> FieldElement<F> {
        if self.vars.is_empty() || p.is_empty() {
            return self.coeff.clone();
        }

        // var_id is index of p
        let eval = self
            .vars
            .iter()
            .fold(FieldElement::<F>::one(), |acc, x| acc * &p[*x]);
        eval * &self.coeff
    }

    /// Assign values to one or more variables in the monomial
    // TODO: can we change this to modify in place to remove the extract allocation?
    //TODO: add valid variable check
    fn partial_evaluate(&self, assignments: &[(usize, FieldElement<F>)]) -> Self {
        //constant term
        if self.vars.is_empty() || assignments.is_empty() {
            //Since term is already constant we return itself
            return self.clone();
        }

        let mut new_coefficient = self.coeff.clone();
        let mut unassigned_variables = self.vars.to_vec();

        for (var_id, assignment) in assignments {
            if unassigned_variables.contains(var_id) {
                new_coefficient = new_coefficient * assignment;
                unassigned_variables.retain(|&id| id != *var_id);
            }
        }

        Self::new((new_coefficient, unassigned_variables))
    }

    fn zero() -> Self {
        Self {
            coeff: FieldElement::<F>::zero(),
            vars: vec![],
        }
    }
}
// Add
impl<F: IsField + IsPrimeField> ops::Add<MultilinearMonomial<F>> for MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn add(self, rhs: MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> ops::Add<&MultilinearMonomial<F>> for &MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn add(self, rhs: &MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> ops::Add<&MultilinearMonomial<F>> for MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn add(self, rhs: &MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> ops::Add<MultilinearMonomial<F>> for &MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn add(self, rhs: MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

// Sub
impl<F: IsField + IsPrimeField> ops::Sub<MultilinearMonomial<F>> for MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn sub(self, rhs: MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> ops::Sub<&MultilinearMonomial<F>> for &MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn sub(self, rhs: &MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> ops::Sub<&MultilinearMonomial<F>> for MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn sub(self, rhs: &MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> ops::Sub<MultilinearMonomial<F>> for &MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn sub(self, rhs: MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

// Mul
impl<F: IsField + IsPrimeField> ops::Mul<MultilinearMonomial<F>> for MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn mul(self, rhs: MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> ops::Mul<&MultilinearMonomial<F>> for &MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn mul(self, rhs: &MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> ops::Mul<&MultilinearMonomial<F>> for MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn mul(self, rhs: &MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> ops::Mul<MultilinearMonomial<F>> for &MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn mul(self, rhs: MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

// Mul: Mono * Field

impl<F: IsField + IsPrimeField> ops::Mul<MultilinearMonomial<F>> for FieldElement<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn mul(self, rhs: MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> ops::Mul<&MultilinearMonomial<F>> for FieldElement<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn mul(self, rhs: &MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> ops::Mul<MultilinearMonomial<F>> for &FieldElement<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn mul(self, rhs: MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> ops::Mul<&MultilinearMonomial<F>> for &FieldElement<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn mul(self, rhs: &MultilinearMonomial<F>) -> Self::Output {
        todo!()
    }
}

// Mul: Field * Mono
impl<F: IsField + IsPrimeField> ops::Mul<FieldElement<F>> for MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn mul(self, rhs: FieldElement<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> ops::Mul<&FieldElement<F>> for MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn mul(self, rhs: &FieldElement<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> ops::Mul<FieldElement<F>> for &MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn mul(self, rhs: FieldElement<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> ops::Mul<&FieldElement<F>> for &MultilinearMonomial<F>
where
    F: IsField,
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = MultilinearMonomial<F>;

    fn mul(self, rhs: &FieldElement<F>) -> Self::Output {
        todo!()
    }
}

impl<F: IsField + IsPrimeField> Display for MultilinearMonomial<F>
where
    <F as IsField>::BaseType: Send + Sync,
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}", self.coeff.representative())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::fields::u64_prime_field::U64PrimeField;

    const ORDER: u64 = 101;
    type F = U64PrimeField<ORDER>;
    type FE = FieldElement<F>;

    #[test]
    fn build_multilinear_monomial() {
        let monomial = MultilinearMonomial::new((FE::new(5), vec![10, 5, 6]));

        // should build and sort the var_id's
        assert_eq!(
            monomial,
            MultilinearMonomial {
                coeff: FE::new(5),
                vars: vec![5, 6, 10]
            }
        );
    }

    #[test]
    fn evaluate_constant_multilinear_monomial() {
        let monomial = MultilinearMonomial::new((FE::new(5), vec![]));

        assert_eq!(monomial.evaluate(&[FE::new(1)]), FE::new(5));
    }

    #[test]
    fn test_partial_evaluation() {
        // 5ab partially evaluate b = 2
        // expected result = 10a
        let five_a_b = MultilinearMonomial::new((FE::new(5), vec![1, 2]));
        let maybe_10_a = five_a_b.partial_evaluate(&[(2, FE::new(2))]);
        assert_eq!(
            maybe_10_a,
            MultilinearMonomial {
                coeff: FE::new(10),
                vars: vec![1]
            }
        );

        // 6abcd evaluate a = 5, c = 3
        // expected = 90bd
        let six_a_b_c_d = MultilinearMonomial::new((FE::new(6), vec![1, 2, 3, 4]));
        let maybe_90_b_d = six_a_b_c_d.partial_evaluate(&[(1, FE::new(5)), (3, FE::new(3))]);
        assert_eq!(
            maybe_90_b_d,
            MultilinearMonomial {
                coeff: FE::new(90),
                vars: vec![2, 4]
            }
        );

        // assign every variable
        // 5ab partially evaluate a= 3, b = 2
        // expected result = 30
        let five_a_b = MultilinearMonomial::new((FE::new(5), vec![1, 2]));
        let maybe_30 = five_a_b.partial_evaluate(&[(1, FE::new(3)), (2, FE::new(2))]);
        assert_eq!(
            maybe_30,
            MultilinearMonomial {
                coeff: FE::new(30),
                vars: vec![]
            }
        );

        // ignore repeated assignments
        // 6abcd evaluate a = 5, c = 3, a = 9
        // expected = 90bd
        // should ignore the second assignment for a, as first already got rid of a
        let six_a_b_c_d = MultilinearMonomial::new((FE::new(6), vec![1, 2, 3, 4]));
        let maybe_90_b_d =
            six_a_b_c_d.partial_evaluate(&[(1, FE::new(5)), (3, FE::new(3)), (1, FE::new(9))]);
        assert_eq!(
            maybe_90_b_d,
            MultilinearMonomial {
                coeff: FE::new(90),
                vars: vec![2, 4]
            }
        );
    }
}