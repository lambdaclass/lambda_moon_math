use crate::{
    field::{element::FieldElement, traits::IsField},
    polynomial::error::MultilinearError,
};
use alloc::{vec, vec::Vec};
use core::ops::{Add, Index, Mul};
#[cfg(feature = "parallel")]
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

/// Represents a multilinear polynomials as a vector of evaluations (FieldElements) in lagrange basis
#[derive(Debug, PartialEq, Clone)]
pub struct DenseMultilinearPolynomial<F: IsField>
where
    <F as IsField>::BaseType: Send + Sync,
{
    evals: Vec<FieldElement<F>>,
    n_vars: usize,
    len: usize,
}

impl<F: IsField> DenseMultilinearPolynomial<F>
where
    <F as IsField>::BaseType: Send + Sync,
{
    /// Build a new multilinear polynomial, from collection of evaluations
    pub fn new(evals: Vec<FieldElement<F>>) -> Self {
        // Pad non-power-2 evaluations to fill out the dense multilinear polynomial
        let mut poly_evals = evals;
        while !(poly_evals.len().is_power_of_two()) {
            poly_evals.push(FieldElement::zero());
        }

        DenseMultilinearPolynomial {
            evals: poly_evals.clone(),
            n_vars: log_2(poly_evals.len()),
            len: poly_evals.len(),
        }
    }

    pub fn num_vars(&self) -> usize {
        self.n_vars
    }

    pub fn evals(&self) -> &Vec<FieldElement<F>> {
        &self.evals
    }

    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.len
    }

    /// Evaluates `self` at the point `p` in O(n) time.
    pub fn evaluate(&self, r: Vec<FieldElement<F>>) -> Result<FieldElement<F>, MultilinearError> {
        // r must have a value for each variable
        if r.len() != self.num_vars() {
            return Err(MultilinearError::IncorrectNumberofEvaluationPoints(
                r.len(),
                self.num_vars(),
            ));
        }

        let mut chis: Vec<FieldElement<F>> =
            vec![FieldElement::one(); (2usize).pow(r.len() as u32)];
        let mut size = 1;
        for j in r {
            size *= 2;
            for i in (0..size).rev().step_by(2) {
                let scalar = &chis[i / 2].clone();
                chis[i] = scalar * &j;
                chis[i - 1] = scalar - &chis[i];
            }
        }
        #[cfg(feature = "parallel")]
        let iter = (0..chis.len()).into_par_iter();

        #[cfg(not(feature = "parallel"))]
        let iter = 0..chis.len();
        Ok(iter.map(|i| &self.evals[i] * &chis[i]).sum())
    }

    pub fn evaluate_with(
        evals: &[FieldElement<F>],
        r: &[FieldElement<F>],
    ) -> Result<FieldElement<F>, MultilinearError> {
        let mut chis: Vec<FieldElement<F>> =
            vec![FieldElement::one(); (2usize).pow(r.len() as u32)];
        if chis.len() != evals.len() {
            return Err(MultilinearError::ChisAndEvalsLengthMismatch(
                chis.len(),
                evals.len(),
            ));
        }
        let mut size = 1;
        for j in r {
            size *= 2;
            for i in (0..size).rev().step_by(2) {
                let scalar = &chis[i / 2].clone();
                chis[i] = scalar * j;
                chis[i - 1] = scalar - &chis[i];
            }
        }
        Ok((0..evals.len()).map(|i| &evals[i] * &chis[i]).sum())
    }

    /// Extends a DenseMultilinearPolynomial by concatenating `other` polynomial of the same length.
    pub fn extend(&mut self, other: &DenseMultilinearPolynomial<F>) {
        debug_assert_eq!(self.evals.len(), self.len);
        let other = other.evals.clone();
        debug_assert_eq!(other.len(), self.len);
        self.evals.extend(other);
        self.n_vars += 1;
        self.len *= 2;
        debug_assert_eq!(self.evals.len(), self.len);
    }

    /// Merges a series of DenseMultilienarPolynomials into one polynomial. Zero-pads the final merged polynomial to the next power_of_two length if necessary.
    pub fn merge(polys: &[DenseMultilinearPolynomial<F>]) -> DenseMultilinearPolynomial<F> {
        // TODO (performance): pre-allocate vector we are resizing two bench to see if it is faster than naively calling resize.
        let mut z: Vec<FieldElement<F>> = Vec::new();
        for poly in polys.iter() {
            z.extend(poly.evals().clone().into_iter());
        }

        // pad the polynomial with zero polynomial at the end
        z.resize(z.len().next_power_of_two(), FieldElement::zero());

        DenseMultilinearPolynomial::new(z)
    }

    pub fn from_u64(evals: &[u64]) -> Self {
        DenseMultilinearPolynomial::new(
            (0..evals.len())
                .map(|i| FieldElement::from(evals[i]))
                .collect::<Vec<FieldElement<F>>>(),
        )
    }

    pub fn scalar_mul(&self, scalar: &FieldElement<F>) -> Self {
        let mut new_poly = self.clone();
        new_poly.evals.iter_mut().for_each(|eval| *eval *= scalar);
        new_poly
    }
}

impl<F: IsField> Index<usize> for DenseMultilinearPolynomial<F>
where
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = FieldElement<F>;

    #[inline(always)]
    fn index(&self, _index: usize) -> &FieldElement<F> {
        &(self.evals[_index])
    }
}

/// Adds another multilinear polynomial to `self`.
/// Assumes the two polynomials have the same number of variables.
impl<F: IsField> Add for DenseMultilinearPolynomial<F>
where
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = Result<Self, &'static str>;

    fn add(self, other: Self) -> Self::Output {
        if self.num_vars() != other.num_vars() {
            return Err("Polynomials must have the same number of variables");
        }

        #[cfg(feature = "parallel")]
        let evals = self.evals.into_par_iter().zip(other.evals.into_par_iter());

        #[cfg(not(feature = "parallel"))]
        let evals = self.evals.iter().zip(other.evals.iter());
        let sum: Vec<FieldElement<F>> = evals.map(|(a, b)| a + b).collect();

        Ok(DenseMultilinearPolynomial::new(sum))
    }
}

impl<F: IsField> Mul<FieldElement<F>> for DenseMultilinearPolynomial<F>
where
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = DenseMultilinearPolynomial<F>;

    fn mul(self, rhs: FieldElement<F>) -> Self::Output {
        Self::scalar_mul(&self, &rhs)
    }
}

impl<F: IsField> Mul<&FieldElement<F>> for DenseMultilinearPolynomial<F>
where
    <F as IsField>::BaseType: Send + Sync,
{
    type Output = DenseMultilinearPolynomial<F>;

    fn mul(self, rhs: &FieldElement<F>) -> Self::Output {
        Self::scalar_mul(&self, rhs)
    }
}

// returns 0 if n is 0
pub fn log_2(n: usize) -> usize {
    if n == 0 {
        return 0;
    }

    if n.is_power_of_two() {
        (1usize.leading_zeros() - n.leading_zeros()) as usize
    } else {
        (0usize.leading_zeros() - n.leading_zeros()) as usize
    }
}

#[cfg(test)]
mod tests {

    use crate::field::fields::u64_prime_field::U64PrimeField;

    use super::*;

    const ORDER: u64 = 101;
    type F = U64PrimeField<ORDER>;
    type FE = FieldElement<F>;

    pub fn evals(r: Vec<FE>) -> Vec<FE> {
        let mut evals: Vec<FE> = vec![FE::one(); (2usize).pow(r.len() as u32)];
        let mut size = 1;
        for j in r {
            size *= 2;
            for i in (0..size).rev().step_by(2) {
                let scalar = evals[i / 2];
                evals[i] = scalar * j;
                evals[i - 1] = scalar - evals[i];
            }
        }
        evals
    }

    pub fn compute_factored_evals(r: Vec<FE>) -> (Vec<FE>, Vec<FE>) {
        let size = r.len();
        let (left_num_vars, _right_num_vars) = (size / 2, size - size / 2);

        let l = evals(r[..left_num_vars].to_vec());
        let r = evals(r[left_num_vars..size].to_vec());

        (l, r)
    }

    fn evaluate_with_lr(z: &[FE], r: &[FE]) -> FE {
        let (l, r) = compute_factored_evals(r.to_vec());

        let size = r.len();
        // ensure size is even
        //TODO
        assert!(size % 2 == 0);
        // n = 2^size
        let n = (2usize).pow(size as u32);
        // compute m = sqrt(n) = 2^{\ell/2}
        let m = (n as f64).sqrt() as usize;

        // compute vector-matrix product between L and Z viewed as a matrix
        let lz = (0..m)
            .map(|i| {
                (0..m).fold(FE::zero(), |mut acc, j| {
                    acc += l[j] * z[j * m + i];
                    acc
                })
            })
            .collect::<Vec<FE>>();

        // compute dot product between LZ and R
        (0..lz.len()).map(|i| lz[i] * r[i]).sum()
    }

    #[test]
    fn evaluation() {
        // Z = [1, 2, 1, 4]
        let z = vec![FE::one(), FE::from(2u64), FE::one(), FE::from(4u64)];

        // r = [4,3]
        let r = vec![FE::from(4u64), FE::from(3u64)];

        let eval_with_lr = evaluate_with_lr(&z, &r);
        let poly = DenseMultilinearPolynomial::new(z);

        let eval = poly.evaluate(r).unwrap();
        assert_eq!(eval, FE::from(28u64));
        assert_eq!(eval_with_lr, eval);
    }

    #[test]
    fn evaluate_with() {
        let two = FE::from(2);

        let z = vec![
            FE::zero(),
            FE::zero(),
            FE::zero(),
            FE::one(),
            FE::one(),
            FE::one(),
            FE::zero(),
            two,
        ];
        let x = vec![FE::one(), FE::one(), FE::one()];

        let y = DenseMultilinearPolynomial::<F>::evaluate_with(z.as_slice(), x.as_slice()).unwrap();
        assert_eq!(y, two);
    }

    #[test]
    fn add() {
        let a = DenseMultilinearPolynomial::new(vec![FE::from(3); 4]);
        let b = DenseMultilinearPolynomial::new(vec![FE::from(7); 4]);

        let c = a.add(b).unwrap();

        assert_eq!(*c.evals(), vec![FE::from(10); 4]);
    }

    #[test]
    fn mul() {
        let a = DenseMultilinearPolynomial::new(vec![FE::from(3); 4]);
        let b = a.mul(&FE::from(2));

        assert_eq!(*b.evals(), vec![FE::from(6); 4]);
    }

    // Take a multilinear polynomial of length 2^2 and merge with a polynomial of 2^1. The resulting polynomial should be padded to len 2^3 = 8 and the last two evals should be FE::zero().
    #[test]
    fn merge() {
        let a = DenseMultilinearPolynomial::new(vec![FE::from(3); 4]);
        let b = DenseMultilinearPolynomial::new(vec![FE::from(3); 2]);

        let c = DenseMultilinearPolynomial::merge(&[a, b]);

        assert_eq!(c.len(), 8);
        assert_eq!(c[c.len() - 1], FE::zero());
        assert_eq!(c[c.len() - 2], FE::zero());
    }

    #[test]
    fn extend() {
        let mut a = DenseMultilinearPolynomial::new(vec![FE::from(3); 4]);
        let b = DenseMultilinearPolynomial::new(vec![FE::from(3); 4]);

        a.extend(&b);

        assert_eq!(a.len(), 8);
        assert_eq!(a.num_vars(), 3);
    }

    #[test]
    #[should_panic]
    fn extend_unequal() {
        let mut a = DenseMultilinearPolynomial::new(vec![FE::from(3); 4]);
        let b = DenseMultilinearPolynomial::new(vec![FE::from(3); 2]);
        a.extend(&b);
    }
}
