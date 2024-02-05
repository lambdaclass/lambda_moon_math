use std::ops::Div;

use crate::domain::Domain;
use crate::frame::Frame;
use crate::prover::evaluate_polynomial_on_lde_domain;
use itertools::Itertools;
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::traits::{IsFFTField, IsField, IsSubFieldOf};
use lambdaworks_math::polynomial::Polynomial;
use num_integer::Integer;

/// TransitionConstraint represents the behaviour that a transition constraint
/// over the computation that wants to be proven must comply with.
pub trait TransitionConstraint<F, E>: Send + Sync
where
    F: IsSubFieldOf<E> + IsFFTField + Send + Sync,
    E: IsField + Send + Sync,
{
    /// The degree of the constraint interpreting it as a multivariate polynomial.
    fn degree(&self) -> usize;

    /// The index of the constraint.
    /// Each transition constraint should have one index in the range [0, N),
    /// where N is the total number of transition constraints.
    fn constraint_idx(&self) -> usize;

    /// The function representing the evaluation of the constraint over elements
    /// of the trace table.
    ///
    /// Elements of the trace table are found in the `frame` input, and depending on the
    /// constraint, elements of `periodic_values` and `rap_challenges` may be used in
    /// the evaluation.
    /// Once computed, the evaluation should be inserted in the `transition_evaluations`
    /// vector, in the index corresponding to the constraint as given by `constraint_idx()`.
    fn evaluate(
        &self,
        frame: &Frame<F, E>,
        transition_evaluations: &mut [FieldElement<E>],
        periodic_values: &[FieldElement<F>],
        rap_challenges: &[FieldElement<E>],
    );

    /// The periodicity the constraint is applied over the trace.
    ///
    /// Default value is 1, meaning that the constraint is applied to every
    /// step of the trace.
    fn period(&self) -> usize {
        1
    }

    /// The offset with respect to the first trace row, where the constraint
    /// is applied.
    /// For example, if the constraint has periodicity 2 and offset 1, this means
    /// the constraint will be applied over trace rows of index 1, 3, 5, etc.
    ///
    /// Default value is 0, meaning that the constraint is applied from the first
    /// element of the trace on.
    fn offset(&self) -> usize {
        0
    }

    /// For a more fine-grained description of where the constraint should apply,
    /// an exemptions period can be defined.
    /// This specifies the periodicity of the row indexes where the constraint should
    /// NOT apply, within the row indexes where the constraint applies, as specified by
    /// `period()` and `offset()`.
    ///
    /// Default value is None.
    fn exemptions_period(&self) -> Option<usize> {
        None
    }

    /// The offset value for periodic exemptions. Check documentation of `period()`,
    /// `offset()` and `exemptions_period` for a better understanding.
    fn periodic_exemptions_offset(&self) -> Option<usize> {
        None
    }

    /// The number of exemptions at the end of the trace.
    ///
    /// This method's output defines what trace elements should not be considered for
    /// the constraint evaluation at the end of the trace. For example, for a fibonacci
    /// computation that has to use the result 2 following steps, this method is defined
    /// to return the value 2.
    fn end_exemptions(&self) -> usize;

    /// Method for calculating the end exemptions polynomial.
    ///
    /// This polynomial is used to compute zerofiers of the constraint, and the default
    /// implementation should normally not be changed.
    fn end_exemptions_poly(
        &self,
        trace_primitive_root: &FieldElement<F>,
        trace_length: usize,
    ) -> Polynomial<FieldElement<F>> {
        let one_poly = Polynomial::new_monomial(FieldElement::<F>::one(), 0);
        if self.end_exemptions() == 0 {
            return one_poly;
        }
        let period = self.period();
        // FIXME: CHECK IF WE NEED TO CHANGE THE NEW MONOMIAL'S ARGUMENTS TO trace_root^(offset * trace_length / period) INSTEAD OF ONE!!!!
        (1..=self.end_exemptions())
            .map(|exemption| trace_primitive_root.pow(trace_length - exemption * period))
            .fold(one_poly, |acc, offset| {
                acc * (Polynomial::new_monomial(FieldElement::<F>::one(), 1) - offset)
            })
    }

    /// Compute evaluations of the constraints zerofier over a LDE domain.
    fn zerofier_evaluations_on_extended_domain(&self, domain: &Domain<F>) -> Vec<FieldElement<F>> {
        let blowup_factor = domain.blowup_factor;
        let trace_length = domain.trace_roots_of_unity.len();
        let trace_primitive_root = &domain.trace_primitive_root;
        let coset_offset = &domain.coset_offset;
        let lde_root_order = u64::from((blowup_factor * trace_length).trailing_zeros());
        let lde_root = F::get_primitive_root_of_unity(lde_root_order).unwrap();

        let end_exemptions_poly = self.end_exemptions_poly(trace_primitive_root, trace_length);

        // If there is an exemptions period defined for this constraint, the evaluations are calculated directly
        // by computing P_exemptions(x) / Zerofier(x)
        if let Some(exemptions_period) = self.exemptions_period() {
            // FIXME: Rather than making this assertions here, it would be better to handle these
            // errors or make these checks when the AIR is initialized.
            debug_assert!(exemptions_period.is_multiple_of(&self.period()));
            debug_assert!(self.periodic_exemptions_offset().is_some());

            // The elements of the domain have order `trace_length * blowup_factor`, so the zerofier evaluations
            // without the end exemptions, repeat their values after `blowup_factor * exemptions_period` iterations,
            // so we only need to compute those.
            let last_exponent = blowup_factor * exemptions_period;

            let evaluations: Vec<_> = (0..last_exponent)
                .map(|exponent| {
                    let x = lde_root.pow(exponent);
                    let offset_times_x = coset_offset * &x;
                    let offset_exponent = trace_length * self.periodic_exemptions_offset().unwrap()
                        / exemptions_period;

                    let numerator = offset_times_x.pow(trace_length / exemptions_period)
                        - trace_primitive_root.pow(offset_exponent);
                    let denominator = offset_times_x.pow(trace_length / self.period())
                        - trace_primitive_root.pow(self.offset() * trace_length / self.period());

                    numerator.div(denominator)
                })
                .collect();

            let end_exemption_evaluations = evaluate_polynomial_on_lde_domain(
                &end_exemptions_poly,
                blowup_factor,
                domain.interpolation_domain_size,
                coset_offset,
            )
            .unwrap();

            let cycled_evaluations = evaluations
                .iter()
                .cycle()
                .take(end_exemption_evaluations.len());

            std::iter::zip(cycled_evaluations, end_exemption_evaluations)
                .map(|(eval, exemption_eval)| eval * exemption_eval)
                .collect()

        // In this else branch, the zerofiers are computed as the numerator, then inverted
        // using batch inverse and then multiplied by P_exemptions(x). This way we don't do
        // useless divisions.
        } else {
            let last_exponent = blowup_factor * self.period();

            let mut evaluations = (0..last_exponent)
                .map(|exponent| {
                    let x = lde_root.pow(exponent);
                    (coset_offset * &x).pow(trace_length / self.period())
                        - trace_primitive_root.pow(self.offset() * trace_length / self.period())
                })
                .collect_vec();

            FieldElement::inplace_batch_inverse(&mut evaluations).unwrap();

            let end_exemption_evaluations = evaluate_polynomial_on_lde_domain(
                &end_exemptions_poly,
                blowup_factor,
                domain.interpolation_domain_size,
                coset_offset,
            )
            .unwrap();

            let cycled_evaluations = evaluations
                .iter()
                .cycle()
                .take(end_exemption_evaluations.len());

            std::iter::zip(cycled_evaluations, end_exemption_evaluations)
                .map(|(eval, exemption_eval)| eval * exemption_eval)
                .collect()
        }
    }

    /// Returns the evaluation of the zerofier corresponding to this constraint in some point
    /// `z`, which could be in a field extension.
    fn evaluate_zerofier(
        &self,
        z: &FieldElement<E>,
        trace_primitive_root: &FieldElement<F>,
        trace_length: usize,
    ) -> FieldElement<E> {
        let end_exemptions_poly = self.end_exemptions_poly(trace_primitive_root, trace_length);

        if let Some(exemptions_period) = self.exemptions_period() {
            debug_assert!(exemptions_period.is_multiple_of(&self.period()));
            debug_assert!(self.periodic_exemptions_offset().is_some());

            let periodic_exemptions_offset = self.periodic_exemptions_offset().unwrap();
            let offset_exponent = trace_length * periodic_exemptions_offset / exemptions_period;

            let numerator = -trace_primitive_root.pow(offset_exponent)
                + z.pow(trace_length / exemptions_period);
            let denominator = -trace_primitive_root
                .pow(self.offset() * trace_length / self.period())
                + z.pow(trace_length / self.period());

            return numerator.div(denominator) * end_exemptions_poly.evaluate(z);
        }

        (-trace_primitive_root.pow(self.offset() * trace_length / self.period())
            + z.pow(trace_length / self.period()))
        .inv()
        .unwrap()
            * end_exemptions_poly.evaluate(z)
    }
}
