use core::fmt::Display;
use std::marker::PhantomData;

use crate::fiat_shamir::is_transcript::IsTranscript;
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::traits::{IsField, IsPrimeField};
use lambdaworks_math::polynomial::{
    dense_multilinear_poly::DenseMultilinearPolynomial, Polynomial,
};
use lambdaworks_math::traits::AsBytes;
use lambdaworks_math::traits::ByteConversion;
#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelIterator, IntoParallelIterator, ParallelIterator};

fn eval_points_quadratic<F: IsField + IsPrimeField, E>(
    poly_a: &DenseMultilinearPolynomial<F>,
    poly_b: &DenseMultilinearPolynomial<F>,
    comb_func: E,
) -> (FieldElement<F>, FieldElement<F>)
where
    <F as IsField>::BaseType: Send + Sync,
    E: Fn(&FieldElement<F>, &FieldElement<F>) -> FieldElement<F> + Sync + Send,
{
    let len = poly_a.len() / 2;
    #[cfg(not(feature = "parallel"))]
    let iter = 0..len;
    #[cfg(feature = "parallel")]
    let iter = (0..len).into_par_iter();
    let res = iter.map(|i| {
        // eval_0: A(low)
        let eval_0 = comb_func(&poly_a[i], &poly_b[i]);

        // eval_2: -A(low) + 2*A(high)
        let poly_a_eval_2 = &poly_a[len + i] + &poly_a[len + i] - &poly_a[i];
        let poly_b_eval_2 = &poly_b[len + i] + &poly_b[len + i] - &poly_b[i];
        let eval_2 = comb_func(&poly_a_eval_2, &poly_b_eval_2);
        (eval_0, eval_2)
    });
    #[cfg(not(feature = "parallel"))]
    let res = res.fold((FieldElement::zero(), FieldElement::zero()), |a, b| {
        (a.0 + b.0, a.1 + b.1)
    });
    #[cfg(feature = "parallel")]
    let res = res.reduce(
        || (FieldElement::zero(), FieldElement::zero()),
        |a, b| (a.0 + b.0, a.1 + b.1),
    );

    res
}

fn eval_points_cubic<F: IsField, E>(
    poly_a: &DenseMultilinearPolynomial<F>,
    poly_b: &DenseMultilinearPolynomial<F>,
    poly_c: &DenseMultilinearPolynomial<F>,
    comb_func: E,
) -> (FieldElement<F>, FieldElement<F>, FieldElement<F>)
where
    <F as IsField>::BaseType: Send + Sync,
    E: Fn(&FieldElement<F>, &FieldElement<F>, &FieldElement<F>) -> FieldElement<F> + Sync,
{
    let len = poly_a.len() / 2;
    #[cfg(not(feature = "parallel"))]
    let iter = 0..len;
    #[cfg(feature = "parallel")]
    let iter = (0..len).into_par_iter();
    let res = iter.map(|i| {
        // eval_0: A(low)
        let eval_0 = comb_func(&poly_a[i], &poly_b[i], &poly_c[i]);

        // eval_2: -A(low) + 2*A(high)
        let poly_a_eval_2 = &poly_a[len + i] + &poly_a[len + i] - &poly_a[i];
        let poly_b_eval_2 = &poly_b[len + i] + &poly_b[len + i] - &poly_b[i];
        let poly_c_eval_2 = &poly_c[len + i] + &poly_c[len + i] - &poly_c[i];
        let eval_2 = comb_func(&poly_a_eval_2, &poly_b_eval_2, &poly_c_eval_2);

        // eval 3: bound_func is -2A(low) + 3A(high); computed incrementally with bound_func applied to eval(2)
        let poly_a_eval_3 = poly_a_eval_2 + &poly_a[len + i] - &poly_a[i];
        let poly_b_eval_3 = poly_b_eval_2 + &poly_b[len + i] - &poly_b[i];
        let poly_c_eval_3 = poly_c_eval_2 + &poly_c[len + i] - &poly_c[i];
        let eval_3 = comb_func(&poly_a_eval_3, &poly_b_eval_3, &poly_c_eval_3);

        (eval_0, eval_2, eval_3)
    });
    #[cfg(not(feature = "parallel"))]
    let res = res.fold(
        (
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
        ),
        |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2),
    );
    #[cfg(feature = "parallel")]
    let res = res.reduce(
        || {
            (
                FieldElement::zero(),
                FieldElement::zero(),
                FieldElement::zero(),
            )
        },
        |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2),
    );

    res
}

#[derive(Debug)]
pub enum SumcheckError {
    InvalidProof,
}

impl Display for SumcheckError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            SumcheckError::InvalidProof => write!(f, "Sumcheck Proof Invalid"),
        }
    }
}

// Proof attesting to sum over the boolean hypercube
#[derive(Debug)]
pub struct SumcheckProof<F: IsPrimeField>
where
    <F as IsField>::BaseType: Send + Sync,
{
    // Sum the proof is attesting to
    pub sum: FieldElement<F>,
    // Univariate polynomial oracles the prover sends to the verifier each round
    pub round_uni_polys: Vec<Polynomial<FieldElement<F>>>,
}

#[derive(Clone, Copy, Debug)]
pub struct Sumcheck<F: IsField + IsPrimeField>
where
    <F as IsField>::BaseType: Send + Sync,
{
    _p: PhantomData<F>,
}

impl<F: IsField + IsPrimeField> Sumcheck<F>
where
    <F as IsField>::BaseType: Send + Sync,
    FieldElement<F>: ByteConversion + AsBytes,
{
    //Used for sum_{(a * b)}
    pub fn prove_quadratic<E>(
        sum: &FieldElement<F>,
        poly_a: &mut DenseMultilinearPolynomial<F>,
        poly_b: &mut DenseMultilinearPolynomial<F>,
        comb_func: E,
        transcript: &mut impl IsTranscript<F>,
    ) -> (SumcheckProof<F>, Vec<FieldElement<F>>)
    where
        E: Fn(&FieldElement<F>, &FieldElement<F>) -> FieldElement<F> + Sync,
    {
        let mut round_uni_polys: Vec<Polynomial<FieldElement<F>>> =
            Vec::with_capacity(poly_a.num_vars());
        let mut challenges = Vec::with_capacity(poly_a.num_vars());
        let mut prev_round_claim = sum.clone();

        for _ in 0..poly_a.num_vars() {
            let round_poly = {
                let (eval_0, eval_2) = eval_points_quadratic(poly_a, poly_b, &comb_func);
                let evals = vec![eval_0.clone(), prev_round_claim - eval_0, eval_2];
                Polynomial::new(&evals)
            };

            // append round's Univariate polynomial to transcript
            transcript.append_bytes(&round_poly.as_bytes());

            // Squeeze Verifier Challenge for next round
            let challenge = &transcript.sample_field_element();
            challenges.push(challenge.clone());

            // compute next claim
            prev_round_claim = round_poly.evaluate(&challenge);

            // add univariate polynomial for this round to the proof
            round_uni_polys.push(round_poly);

            // fix next variable of poly
            poly_a.fix_variable(&challenge);
            poly_b.fix_variable(&challenge);
        }

        (
            SumcheckProof {
                sum: sum.clone(),
                round_uni_polys,
            },
            challenges,
        )
    }

    pub fn prove_quadratic_batched<E>(
        sum: &FieldElement<F>,
        poly_a: &mut Vec<DenseMultilinearPolynomial<F>>,
        poly_b: &mut Vec<DenseMultilinearPolynomial<F>>,
        powers: Option<&[FieldElement<F>]>,
        comb_func: E,
        transcript: &mut impl IsTranscript<F>,
    ) -> SumcheckProof<F>
    where
        E: Fn(&FieldElement<F>, &FieldElement<F>) -> FieldElement<F> + Sync,
    {
        let mut round_uni_polys: Vec<Polynomial<FieldElement<F>>> =
            Vec::with_capacity(poly_a[0].num_vars());
        let mut challenges = Vec::with_capacity(poly_a[0].num_vars());
        let mut prev_round_claim = sum.clone();

        for _ in 0..poly_a[0].num_vars() {
            let mut evals: Vec<(FieldElement<F>, FieldElement<F>)> = Vec::new();

            for (poly_a, poly_b) in poly_a.iter().zip(poly_b.iter()) {
                let (eval_point_0, eval_point_2) =
                    eval_points_quadratic(poly_a, poly_b, &comb_func);
                evals.push((eval_point_0, eval_point_2));
            }

            // TODO: make optional as we want to perform a batched check outside of this
            let evals_combined_0: FieldElement<F>;
            let evals_combined_2;
            if let Some(powers) = powers {
                evals_combined_0 = (0..evals.len()).map(|i| &evals[i].0 * &powers[i]).sum();
                evals_combined_2 = (0..evals.len()).map(|i| &evals[i].1 * &powers[i]).sum();
            } else {
                //TODO: Implement Sum
                evals_combined_0 = (0..evals.len()).map(|i| evals[i].0.clone()).sum();
                evals_combined_2 = (0..evals.len()).map(|i| evals[i].1.clone()).sum();
            }

            let evals = vec![
                evals_combined_0.clone(),
                prev_round_claim - evals_combined_0,
                evals_combined_2,
            ];
            let round_poly = Polynomial::new(&evals);

            // TODO append the prover's message to the transcript
            transcript.append_bytes(&round_poly.as_bytes());

            // Squeeze Verifier Challenge for next round
            let challenge = &transcript.sample_field_element();
            challenges.push(challenge.clone());

            // bound all tables to the verifier's challenege
            for (poly_a, poly_b) in poly_a.iter_mut().zip(poly_b.iter_mut()) {
                poly_a.fix_variable(&challenge);
                poly_b.fix_variable(&challenge);
            }

            prev_round_claim = round_poly.evaluate(&challenge);
            round_uni_polys.push(round_poly);
        }

        SumcheckProof {
            sum: sum.clone(),
            round_uni_polys,
        }
    }

    pub fn prove_cubic<E>(
        sum: &FieldElement<F>,
        poly_a: &mut DenseMultilinearPolynomial<F>,
        poly_b: &mut DenseMultilinearPolynomial<F>,
        poly_c: &mut DenseMultilinearPolynomial<F>,
        comb_func: E,
        transcript: &mut impl IsTranscript<F>,
    ) -> (SumcheckProof<F>, Vec<FieldElement<F>>)
    where
        E: Fn(&FieldElement<F>, &FieldElement<F>, &FieldElement<F>) -> FieldElement<F> + Sync,
    {
        let mut round_uni_polys: Vec<Polynomial<FieldElement<F>>> =
            Vec::with_capacity(poly_a.num_vars());
        let mut challenges = Vec::with_capacity(poly_a.num_vars());
        let mut prev_round_claim = sum.clone();

        for _ in 0..poly_a.num_vars() {
            let round_poly = {
                let (eval_point_0, eval_point_2, eval_point_3) =
                    eval_points_cubic(poly_a, poly_b, poly_c, &comb_func);
                let evals = vec![
                    eval_point_0.clone(),
                    prev_round_claim - eval_point_0,
                    eval_point_2,
                    eval_point_3,
                ];
                Polynomial::new(&evals)
            };

            // TODO append the prover's message to the transcript
            transcript.append_bytes(&round_poly.as_bytes());

            // Squeeze Verifier Challenge for next round
            let challenge = transcript.sample_field_element();
            challenges.push(challenge.clone());

            // bound all tables to the verifier's challenege
            poly_a.fix_variable(&challenge);
            poly_b.fix_variable(&challenge);
            poly_c.fix_variable(&challenge);

            prev_round_claim = round_poly.evaluate(&challenge);
            round_uni_polys.push(round_poly);
        }

        (
            SumcheckProof {
                sum: sum.clone(),
                round_uni_polys,
            },
            challenges,
        )
    }

    pub fn prove_cubic_batched<E>(
        sum: &FieldElement<F>,
        poly_a: &mut Vec<DenseMultilinearPolynomial<F>>,
        poly_b: &mut Vec<DenseMultilinearPolynomial<F>>,
        poly_c: &DenseMultilinearPolynomial<F>,
        powers: Option<&[FieldElement<F>]>,
        comb_func: E,
        transcript: &mut impl IsTranscript<F>,
    ) -> (SumcheckProof<F>, Vec<FieldElement<F>>)
    where
        E: Fn(&FieldElement<F>, &FieldElement<F>, &FieldElement<F>) -> FieldElement<F> + Sync,
    {
        let mut round_uni_polys: Vec<Polynomial<FieldElement<F>>> =
            Vec::with_capacity(poly_a[0].num_vars());
        let mut challenges = Vec::with_capacity(poly_a[0].num_vars());
        let mut prev_round_claim = sum.clone();

        for _ in 0..poly_a[0].num_vars() {
            let mut evals: Vec<(FieldElement<F>, FieldElement<F>, FieldElement<F>)> = Vec::new();

            for (poly_a, poly_b) in poly_a.iter().zip(poly_b.iter()) {
                let (eval_point_0, eval_point_2, eval_point_3) =
                    eval_points_cubic(poly_a, poly_b, poly_c, &comb_func);
                evals.push((eval_point_0, eval_point_2, eval_point_3));
            }

            // TODO: make optional as we want to perform a batched check outside of this
            let evals_combined_0: FieldElement<F>;
            let evals_combined_2: FieldElement<F>;
            let evals_combined_3: FieldElement<F>;
            if let Some(powers) = powers {
                evals_combined_0 = (0..evals.len()).map(|i| &evals[i].0 * &powers[i]).sum();
                evals_combined_2 = (0..evals.len()).map(|i| &evals[i].1 * &powers[i]).sum();
                evals_combined_3 = (0..evals.len()).map(|i| &evals[i].2 * &powers[i]).sum();
            } else {
                evals_combined_0 = (0..evals.len()).map(|i| evals[i].0.clone()).sum();
                evals_combined_2 = (0..evals.len()).map(|i| evals[i].1.clone()).sum();
                evals_combined_3 = (0..evals.len()).map(|i| evals[i].2.clone()).sum();
            }

            let evals = vec![
                evals_combined_0.clone(),
                prev_round_claim - evals_combined_0,
                evals_combined_2,
                evals_combined_3,
            ];
            let round_poly = Polynomial::new(&evals);

            // TODO: Check if order matters
            transcript.append_bytes(&round_poly.as_bytes());

            // Squeeze Verifier Challenge for next round
            let challenge = &transcript.sample_field_element();
            challenges.push(challenge.clone());

            // TODO: rayon::join and gate
            // bound all tables to the verifier's challenege
            for (poly_a, poly_b) in poly_a.iter_mut().zip(poly_b.iter_mut()) {
                poly_a.fix_variable(&challenge);
                poly_b.fix_variable(&challenge);
            }

            prev_round_claim = round_poly.evaluate(&challenge);
            round_uni_polys.push(round_poly);
        }

        (
            SumcheckProof {
                sum: sum.clone(),
                round_uni_polys,
            },
            challenges,
        )
    }

    // Special instance of sumcheck for a cubic polynomial with an additional additive term:
    // this is used in Spartan: (a * ((b * c) - d))
    pub fn prove_cubic_additive_term<E>(
        sum: &FieldElement<F>,
        poly_a: &mut DenseMultilinearPolynomial<F>,
        poly_b: &mut DenseMultilinearPolynomial<F>,
        poly_c: &mut DenseMultilinearPolynomial<F>,
        poly_d: &mut DenseMultilinearPolynomial<F>,
        comb_func: E,
        transcript: &mut impl IsTranscript<F>,
    ) -> (SumcheckProof<F>, Vec<FieldElement<F>>)
    where
        E: Fn(
                &FieldElement<F>,
                &FieldElement<F>,
                &FieldElement<F>,
                &FieldElement<F>,
            ) -> FieldElement<F>
            + Sync,
    {
        let mut round_uni_polys: Vec<Polynomial<FieldElement<F>>> =
            Vec::with_capacity(poly_a.num_vars());
        let mut challenges = Vec::with_capacity(poly_a.num_vars());
        let mut prev_round_claim = sum.clone();

        for _ in 0..poly_a.num_vars() {
            let round_poly = {
                let (eval_point_0, eval_point_2, eval_point_3) = {
                    //TODO: remove this dedup if possible
                    let len = poly_a.len() / 2;
                    #[cfg(not(feature = "parallel"))]
                    let iter = 0..len;
                    #[cfg(feature = "parallel")]
                    let iter = (0..len).into_par_iter();
                    let res = iter.map(|i| {
                        // eval 0: bound_func is A(low)
                        let eval_point_0 =
                            comb_func(&poly_a[i], &poly_b[i], &poly_c[i], &poly_d[i]);

                        // eval 2: bound_func is -A(low) + 2*A(high)
                        let poly_a_point_2 = &poly_a[len + i] + &poly_a[len + i] - &poly_a[i];
                        let poly_b_point_2 = &poly_b[len + i] + &poly_b[len + i] - &poly_b[i];
                        let poly_c_point_2 = &poly_c[len + i] + &poly_c[len + i] - &poly_c[i];
                        let poly_d_point_2 = &poly_d[len + i] + &poly_d[len + i] - &poly_c[i];
                        let eval_point_2 = comb_func(
                            &poly_a_point_2,
                            &poly_b_point_2,
                            &poly_c_point_2,
                            &poly_d_point_2,
                        );

                        // eval 3: bound_func is -2A(low) + 3A(high); computed incrementally with bound_func applied to eval(2)
                        let poly_a_point_3 = poly_a_point_2 + &poly_a[len + i] - &poly_a[i];
                        let poly_b_point_3 = poly_b_point_2 + &poly_b[len + i] - &poly_b[i];
                        let poly_c_point_3 = poly_c_point_2 + &poly_c[len + i] - &poly_c[i];
                        let poly_d_point_3 = poly_d_point_2 + &poly_d[len + i] - &poly_d[i];
                        let eval_point_3 = comb_func(
                            &poly_a_point_3,
                            &poly_b_point_3,
                            &poly_c_point_3,
                            &poly_d_point_3,
                        );
                        (eval_point_0, eval_point_2, eval_point_3)
                    });
                    #[cfg(not(feature = "parallel"))]
                    let res = res.fold(
                        (
                            FieldElement::zero(),
                            FieldElement::zero(),
                            FieldElement::zero(),
                        ),
                        |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2),
                    );
                    #[cfg(feature = "parallel")]
                    let res = res.reduce(
                        || {
                            (
                                FieldElement::zero(),
                                FieldElement::zero(),
                                FieldElement::zero(),
                            )
                        },
                        |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2),
                    );
                    res
                };
                let evals = vec![
                    eval_point_0.clone(),
                    prev_round_claim - eval_point_0,
                    eval_point_2,
                    eval_point_3,
                ];
                Polynomial::new(&evals)
            };

            // TODO: Does it matter that its before the challenge???? -> Should be I believe
            transcript.append_bytes(&round_poly.as_bytes());

            // Squeeze Verifier Challenge for next round
            let challenge = transcript.sample_field_element();
            challenges.push(challenge.clone());

            prev_round_claim = round_poly.evaluate(&challenge);
            round_uni_polys.push(round_poly);

            // TODO: rayon::join and gate
            // bound all tables to the verifier's challenege
            poly_a.fix_variable(&challenge);
            poly_b.fix_variable(&challenge);
            poly_c.fix_variable(&challenge);
            poly_d.fix_variable(&challenge);
        }

        (
            SumcheckProof {
                sum: sum.clone(),
                round_uni_polys,
            },
            challenges,
        )
    }

    // Create a test for this
    pub fn prove_single(
        poly: &mut DenseMultilinearPolynomial<F>,
        sum: &FieldElement<F>,
        transcript: &mut impl IsTranscript<F>,
    ) -> (SumcheckProof<F>, Vec<FieldElement<F>>) {
        let mut round_uni_polys: Vec<Polynomial<FieldElement<F>>> =
            Vec::with_capacity(poly.num_vars());
        let mut challenges = Vec::with_capacity(poly.num_vars());

        let mut prev_round_claim = sum.clone();

        // Number round = num vars
        for _ in 0..poly.num_vars() {
            // Compute evaluation points of the Dense Multilinear Poly
            let round_poly = {
                let mle_half = poly.len() / 2;
                let eval_0: FieldElement<F> = (0..mle_half).map(|i| poly[i].clone()).sum();
                // We evaluate the poly at each round and each random challenge at 0, 1 we can compute both of these evaluations by summing over the boolearn hypercube for 0, 1 at the fixed point
                // An additional optimization is to sum over eval_0 and compute eval_1 = prev_round_claim - eval_0;
                let evals = vec![eval_0.clone(), prev_round_claim - eval_0];
                Polynomial::new(&evals)
            };

            // TODO: Append poly to transcript -> Modify Transcript
            transcript.append_bytes(&round_poly.as_bytes());

            let challenge = &transcript.sample_field_element();
            challenges.push(challenge.clone());

            // grab next claim
            prev_round_claim = round_poly.evaluate(&challenge);

            // add univariate polynomial for this round to the proof
            round_uni_polys.push(round_poly);

            // takes mutable reference and fixes poly at challenge
            // On each round we evaluate over the hypercube to generate the univariate polynomial for this round. Then we fix the challenge for the next variable,
            // reassign and start the next round with the fixed variable. Each round the poly decreases in size
            poly.fix_variable(&challenge);
        }

        (
            SumcheckProof {
                sum: sum.clone(),
                round_uni_polys,
            },
            challenges,
        )
    }

    // Verifies a sumcheck proof returning the claimed evaluation and random points used during sumcheck rounds
    /// Note: Verification does not execute the final check of sumcheck protocol: g_v(r_v) = oracle_g(r),
    /// as the oracle is not passed in. Expected that the caller will implement.
    ///
    pub fn verify(
        proof: SumcheckProof<F>,
        num_vars: usize,
        transcript: &mut impl IsTranscript<F>,
    ) -> Result<(FieldElement<F>, Vec<FieldElement<F>>), SumcheckError> {
        let mut e = proof.sum.clone();
        let mut r: Vec<FieldElement<F>> = Vec::with_capacity(num_vars);

        // verify there is a univariate polynomial for each round
        if proof.round_uni_polys.len() != num_vars {
            return Err(SumcheckError::InvalidProof);
        }

        for poly in proof.round_uni_polys {

            // check if G_k(0) + G_k(1) = e
            if poly.eval_at_one() + poly.eval_at_zero() != e
            {
                println!("Oh No");
                return Err(SumcheckError::InvalidProof);
            }
            transcript.append_bytes(&poly.as_bytes());

            let challenge = &transcript.sample_field_element();
            r.push(challenge.clone());

            e = poly.evaluate(&challenge);
        }
        Ok((proof.sum, r))
    }
}

#[cfg(test)]
mod test {
    use crate::fiat_shamir::default_transcript::DefaultTranscript;
    use crate::subprotocols::sumcheck::Sumcheck;
    use lambdaworks_math::field::element::FieldElement;
    use lambdaworks_math::field::fields::fft_friendly::u64_goldilocks::U64GoldilocksPrimeField;
    use lambdaworks_math::field::traits::IsField;
    use lambdaworks_math::polynomial::dense_multilinear_poly::DenseMultilinearPolynomial;

    type F = U64GoldilocksPrimeField;

    pub fn index_to_field_bitvector<F: IsField>(value: usize, bits: usize) -> Vec<FieldElement<F>> {
        let mut vec: Vec<FieldElement<F>> = Vec::with_capacity(bits);

        for i in (0..bits).rev() {
            if (value >> i) & 1 == 1 {
                vec.push(FieldElement::one());
            } else {
                vec.push(FieldElement::zero());
            }
        }
        vec
    }

    #[test]
    fn prove_cubic() {
        // Create three dense polynomials (all the same)
        let num_vars = 3;
        let num_evals = (2usize).pow(num_vars as u32);
        let mut evals: Vec<FieldElement<F>> = Vec::with_capacity(num_evals);
        for i in 0..num_evals {
            evals.push(FieldElement::from(8 + i as u64));
        }

        let mut a: DenseMultilinearPolynomial<F> = DenseMultilinearPolynomial::new(evals.clone());
        let mut b: DenseMultilinearPolynomial<F> = DenseMultilinearPolynomial::new(evals.clone());
        let mut c: DenseMultilinearPolynomial<F> = DenseMultilinearPolynomial::new(evals.clone());

        let mut claim = FieldElement::<F>::zero();
        for i in 0..num_evals {
            claim += a.evaluate(index_to_field_bitvector(i, num_vars)).unwrap()
                * b.evaluate(index_to_field_bitvector(i, num_vars)).unwrap()
                * c.evaluate(index_to_field_bitvector(i, num_vars)).unwrap();
        }

        let comb_func_prod = |a: &FieldElement<F>,
                              b: &FieldElement<F>,
                              c: &FieldElement<F>|
         -> FieldElement<F> { a * b * c };

        let r = vec![
            FieldElement::from(3),
            FieldElement::from(1),
            FieldElement::from(3),
        ]; // point 0,0,0 within the boolean hypercube

        let mut transcript = DefaultTranscript::new(b"prove_cubic");
        let (proof, challenges) = Sumcheck::<F>::prove_cubic(
            &claim,
            &mut a,
            &mut b,
            &mut c,
            comb_func_prod,
            &mut transcript,
        );

        let mut transcript = DefaultTranscript::new(b"prove cubic");
        let verify_result = Sumcheck::verify(proof, num_vars, &mut transcript);
        assert!(verify_result.is_ok());

        let (verify_evaluation, verify_randomness) = verify_result.unwrap();
        assert_eq!(challenges, verify_randomness);
        assert_eq!(challenges, r);

        let a = a.evaluate(challenges.clone()).unwrap();
        let b = b.evaluate(challenges.clone()).unwrap();
        let c = c.evaluate(challenges.clone()).unwrap();

        let oracle_query = a * b * c;
        assert_eq!(verify_evaluation, oracle_query);
    }

    #[test]
    #[ignore]
    fn prove_cubic_batched() {}

    #[test]
    #[ignore]
    fn prove_cubic_additive() {}

    #[test]
    fn prove_quad() {
        let num_vars = 3;
        let num_evals = (2usize).pow(num_vars as u32);
        let mut evals: Vec<FieldElement<F>> = Vec::with_capacity(num_evals);
        for i in 0..num_evals {
            evals.push(FieldElement::from(8 + i as u64));
        }

        let mut a: DenseMultilinearPolynomial<F> = DenseMultilinearPolynomial::new(evals.clone());
        let mut b: DenseMultilinearPolynomial<F> = DenseMultilinearPolynomial::new(evals.clone());

        let mut claim = FieldElement::<F>::zero();
        for i in 0..num_evals {
            claim += a.evaluate(index_to_field_bitvector(i, num_vars)).unwrap()
                * b.evaluate(index_to_field_bitvector(i, num_vars)).unwrap();
        }

        let comb_func_prod =
            |a: &FieldElement<F>, b: &FieldElement<F>| -> FieldElement<F> { a * b };

        /*
        let r = vec![
            FieldElement::from(3),
            FieldElement::from(1),
            FieldElement::from(3),
        ]; // point 0,0,0 within the boolean hypercube
        */

        let mut transcript = DefaultTranscript::new(b"prove_quad");
        let (proof, challenges) =
            Sumcheck::<F>::prove_quadratic(&claim, &mut a, &mut b, comb_func_prod, &mut transcript);

        let mut transcript = DefaultTranscript::new(b"prove_quad");
        let verify = Sumcheck::verify(proof, num_vars,  &mut transcript).unwrap();

        /*
        let (verify_evaluation, verify_randomness) = verify_result.unwrap();
        assert_eq!(challenges, verify_randomness);
        assert_eq!(challenges, r);

        // Consider this the opening proof to a(r) * b(r)
        let a = a.evaluate(challenges.clone()).unwrap();
        let b = b.evaluate(challenges).unwrap();

        let oracle_query = a * b;
        assert_eq!(verify_evaluation, oracle_query);
        */
    }

    #[test]
    #[ignore]
    fn prove_quad_batched() {}

    #[test]
    #[ignore]
    fn prove_single() {
        let num_vars = 3;
        let num_evals = (2usize).pow(num_vars as u32);
        let mut evals: Vec<FieldElement<F>> = Vec::with_capacity(num_evals);
        for i in 0..num_evals {
            evals.push(FieldElement::from(8 + i as u64));
        }

        let mut a: DenseMultilinearPolynomial<F> = DenseMultilinearPolynomial::new(evals.clone());

        let mut claim = FieldElement::<F>::zero();
        for i in 0..num_evals {
            claim += a.evaluate(index_to_field_bitvector(i, num_vars)).unwrap()
        }

        let r = vec![
            FieldElement::from(3),
            FieldElement::from(1),
            FieldElement::from(3),
        ]; // point 0,0,0 within the boolean hypercube

        let mut transcript = DefaultTranscript::new(b"prove_single");
        let (proof, challenges) = Sumcheck::<F>::prove_single(&mut a, &claim, &mut transcript);

        let mut transcript = DefaultTranscript::new(b"prove_single");
        let verify_result = Sumcheck::verify(proof, a.num_vars(), &mut transcript);
        assert!(verify_result.is_ok());

        let (verify_evaluation, verify_randomness) = verify_result.unwrap();
        assert_eq!(challenges, verify_randomness);
        assert_eq!(challenges, r);

        assert_eq!(verify_evaluation, a.evaluate(challenges).unwrap());
    }
}
