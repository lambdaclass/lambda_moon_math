use lambdaworks_crypto::fiat_shamir::transcript::Transcript;
use lambdaworks_math::errors::DeserializationError;
use lambdaworks_math::field::traits::IsFFTField;
use lambdaworks_math::traits::{AsBytes, Deserializable, IsRandomFieldElementGenerator};
use std::borrow::Borrow;
use std::marker::PhantomData;
use std::mem::size_of;

use crate::setup::{
    new_strong_fiat_shamir_transcript, CommonPreprocessedInput, VerificationKey, Witness,
};
use lambdaworks_crypto::commitments::traits::IsPolynomialCommitmentScheme;
use lambdaworks_math::{
    field::element::FieldElement,
    polynomial::{self, Polynomial},
};
use lambdaworks_math::{field::traits::IsField, traits::ByteConversion};

/// Plonk proof.
/// The challenges are denoted
///     Round 2: β,γ,
///     Round 3: α,
///     Round 4: ζ,
///     Round 5: υ.
/// Here `Z_H` denotes the domain polynomial, `z` is the polynomial
/// that encodes the copy constraints, and `p` is the sum of `z` and
/// the polynomial that encodes the gates constraints.
/// The polynomial `t` is defined as `p / Z_H`.
/// `a`, `b`, and `c` are the wire assignment polynomials.
/// `S_σ1(ζ), S_σ2(ζ) and S_σ3(ζ)` are the copy permutation polynomials.
/// The polynomial `p` can be "linearized" and the result can be written as
/// `linearized_p = p_non_constant + p_constant`, where
/// `p_non_constant` is the sum of all the terms with a "non-constant"
/// polynomial factor, such as `b(ζ)Q_R(X)`, and `p_constant` is the
/// sum of all the rest (such as `PI(ζ)`).
pub struct Proof<F: IsField, CS: IsPolynomialCommitmentScheme<F>> {
    // Round 1.
    /// Commitment to the wire polynomial `a(x)`
    pub a_1: CS::Commitment,
    /// Commitment to the wire polynomial `b(x)`
    pub b_1: CS::Commitment,
    /// Commitment to the wire polynomial `c(x)`
    pub c_1: CS::Commitment,

    // Round 2.
    /// Commitment to the copy constraints polynomial `z(x)`
    pub z_1: CS::Commitment,

    // Round 3.
    /// Commitment to the low part of the quotient polynomial t(X)
    pub t_lo_1: CS::Commitment,
    /// Commitment to the middle part of the quotient polynomial t(X)
    pub t_mid_1: CS::Commitment,
    /// Commitment to the high part of the quotient polynomial t(X)
    pub t_hi_1: CS::Commitment,

    // Round 4.
    /// Value of `a(ζ)`.
    pub a_zeta: FieldElement<F>,
    /// Value of `b(ζ)`.
    pub b_zeta: FieldElement<F>,
    /// Value of `c(ζ)`.
    pub c_zeta: FieldElement<F>,
    /// Value of `S_σ1(ζ)`.
    pub s1_zeta: FieldElement<F>,
    /// Value of `S_σ2(ζ)`.
    pub s2_zeta: FieldElement<F>,
    /// Value of `z(ζω)`.
    pub z_zeta_omega: FieldElement<F>,

    // Round 5
    /// Value of `p_non_constant(ζ)`.
    pub p_non_constant_zeta: FieldElement<F>,
    ///  Value of `t(ζ)`.
    pub t_zeta: FieldElement<F>,
    /// Batch opening proof for all the evaluations at ζ
    pub w_zeta_1: CS::Proof,
    /// Single opening proof for `z(ζω)`.
    pub w_zeta_omega_1: CS::Proof,
}

impl<F, CS> AsBytes for Proof<F, CS>
where
    F: IsField,
    CS: IsPolynomialCommitmentScheme<F>,
    FieldElement<F>: ByteConversion,
    CS::Commitment: AsBytes,
    CS::Proof: AsBytes,
{
    fn as_bytes(&self) -> Vec<u8> {
        let field_elements = [
            &self.a_zeta,
            &self.b_zeta,
            &self.c_zeta,
            &self.s1_zeta,
            &self.s2_zeta,
            &self.z_zeta_omega,
            &self.p_non_constant_zeta,
            &self.t_zeta,
        ];
        let commitments = [
            &self.a_1,
            &self.b_1,
            &self.c_1,
            &self.z_1,
            &self.t_lo_1,
            &self.t_mid_1,
            &self.t_hi_1,
        ];

        let proofs = [&self.w_zeta_1, &self.w_zeta_omega_1];

        let mut serialized_proof: Vec<u8> = Vec::new();

        field_elements.iter().for_each(|element| {
            let serialized_element = element.to_bytes_be();
            serialized_proof.extend_from_slice(&(serialized_element.len() as u32).to_be_bytes());
            serialized_proof.extend_from_slice(&serialized_element);
        });

        commitments.iter().for_each(|commitment| {
            let serialized_commitment = commitment.as_bytes();
            serialized_proof.extend_from_slice(&(serialized_commitment.len() as u32).to_be_bytes());
            serialized_proof.extend_from_slice(&serialized_commitment);
        });

        proofs.iter().for_each(|proof| {
            let serialized_prf = proof.as_bytes();
            serialized_proof.extend_from_slice(&(serialized_prf.len() as u32).to_be_bytes());
            serialized_proof.extend_from_slice(&serialized_prf);
        });

        serialized_proof
    }
}

// TODO: Remove this once FieldElements implement Serializable
fn deserialize_field_element<F>(
    bytes: &[u8],
    offset: usize,
) -> Result<(usize, FieldElement<F>), DeserializationError>
where
    F: IsField,
    FieldElement<F>: ByteConversion,
{
    let mut offset = offset;
    let element_size_bytes: [u8; size_of::<u32>()] = bytes
        .get(offset..offset + size_of::<u32>())
        .ok_or(DeserializationError::InvalidAmountOfBytes)?
        .try_into()
        .map_err(|_| DeserializationError::InvalidAmountOfBytes)?;
    let element_size = u32::from_be_bytes(element_size_bytes) as usize;
    offset += size_of::<u32>();
    let field_element = FieldElement::from_bytes_be(
        bytes
            .get(offset..offset + element_size)
            .ok_or(DeserializationError::InvalidAmountOfBytes)?,
    )?;
    offset += element_size;
    Ok((offset, field_element))
}

fn deserialize_commitment<Commitment>(
    bytes: &[u8],
    offset: usize,
) -> Result<(usize, Commitment), DeserializationError>
where
    Commitment: Deserializable,
{
    let mut offset = offset;
    let element_size_bytes: [u8; size_of::<u32>()] = bytes
        .get(offset..offset + size_of::<u32>())
        .ok_or(DeserializationError::InvalidAmountOfBytes)?
        .try_into()
        .map_err(|_| DeserializationError::InvalidAmountOfBytes)?;
    let element_size = u32::from_be_bytes(element_size_bytes) as usize;
    offset += size_of::<u32>();
    let commitment = Commitment::deserialize(
        bytes
            .get(offset..offset + element_size)
            .ok_or(DeserializationError::InvalidAmountOfBytes)?,
    )?;
    offset += element_size;
    Ok((offset, commitment))
}

impl<F, CS> Deserializable for Proof<F, CS>
where
    F: IsField,
    CS: IsPolynomialCommitmentScheme<F>,
    FieldElement<F>: ByteConversion,
    CS::Commitment: Deserializable,
    CS::Proof: Deserializable,
{
    fn deserialize(bytes: &[u8]) -> Result<Self, DeserializationError>
    where
        Self: Sized,
    {
        let (offset, a_zeta) = deserialize_field_element(bytes, 0)?;
        let (offset, b_zeta) = deserialize_field_element(bytes, offset)?;
        let (offset, c_zeta) = deserialize_field_element(bytes, offset)?;
        let (offset, s1_zeta) = deserialize_field_element(bytes, offset)?;
        let (offset, s2_zeta) = deserialize_field_element(bytes, offset)?;
        let (offset, z_zeta_omega) = deserialize_field_element(bytes, offset)?;
        let (offset, p_non_constant_zeta) = deserialize_field_element(bytes, offset)?;
        let (offset, t_zeta) = deserialize_field_element(bytes, offset)?;

        let (offset, a_1) = deserialize_commitment(bytes, offset)?;
        let (offset, b_1) = deserialize_commitment(bytes, offset)?;
        let (offset, c_1) = deserialize_commitment(bytes, offset)?;
        let (offset, z_1) = deserialize_commitment(bytes, offset)?;
        let (offset, t_lo_1) = deserialize_commitment(bytes, offset)?;
        let (offset, t_mid_1) = deserialize_commitment(bytes, offset)?;
        let (offset, t_hi_1) = deserialize_commitment(bytes, offset)?;
        let (offset, w_zeta_1) = deserialize_commitment(bytes, offset)?;
        let (_, w_zeta_omega_1) = deserialize_commitment(bytes, offset)?;

        Ok(Proof {
            a_1,
            b_1,
            c_1,
            z_1,
            t_lo_1,
            t_mid_1,
            t_hi_1,
            a_zeta,
            b_zeta,
            c_zeta,
            s1_zeta,
            s2_zeta,
            z_zeta_omega,
            p_non_constant_zeta,
            t_zeta,
            w_zeta_1,
            w_zeta_omega_1,
        })
    }
}

pub struct Prover<
    F: IsField,
    CS: IsPolynomialCommitmentScheme<F>,
    R: IsRandomFieldElementGenerator<F>,
> {
    commitment_scheme: CS,
    random_generator: R,
    phantom: PhantomData<F>,
}

struct Round1Result<F: IsField, Hiding> {
    a_1: Hiding,
    b_1: Hiding,
    c_1: Hiding,
    p_a: Polynomial<FieldElement<F>>,
    p_b: Polynomial<FieldElement<F>>,
    p_c: Polynomial<FieldElement<F>>,
}

struct Round2Result<F: IsField, Hiding> {
    z_1: Hiding,
    p_z: Polynomial<FieldElement<F>>,
    beta: FieldElement<F>,
    gamma: FieldElement<F>,
}

struct Round3Result<F: IsField, Hiding> {
    t_lo_1: Hiding,
    t_mid_1: Hiding,
    t_hi_1: Hiding,
    p_t_lo: Polynomial<FieldElement<F>>,
    p_t_mid: Polynomial<FieldElement<F>>,
    p_t_hi: Polynomial<FieldElement<F>>,
    alpha: FieldElement<F>,
}

struct Round4Result<F: IsField> {
    a_zeta: FieldElement<F>,
    b_zeta: FieldElement<F>,
    c_zeta: FieldElement<F>,
    s1_zeta: FieldElement<F>,
    s2_zeta: FieldElement<F>,
    z_zeta_omega: FieldElement<F>,
    zeta: FieldElement<F>,
}

struct Round5Result<F: IsField, Hiding> {
    w_zeta_1: Hiding,
    w_zeta_omega_1: Hiding,
    p_non_constant_zeta: FieldElement<F>,
    t_zeta: FieldElement<F>,
}

impl<F, CS, R> Prover<F, CS, R>
where
    F: IsField + IsFFTField,
    CS: IsPolynomialCommitmentScheme<F>,
    FieldElement<F>: ByteConversion,
    CS::Commitment: AsBytes,
    R: IsRandomFieldElementGenerator<F>,
{
    pub fn new(commitment_scheme: CS, random_generator: R) -> Self {
        Self {
            commitment_scheme,
            random_generator,
            phantom: PhantomData,
        }
    }

    fn blind_polynomial(
        &self,
        target: &Polynomial<FieldElement<F>>,
        blinder: &Polynomial<FieldElement<F>>,
        n: u64,
    ) -> Polynomial<FieldElement<F>>
    where
        F: IsField,
        R: IsRandomFieldElementGenerator<F>,
    {
        let bs: Vec<FieldElement<F>> = (0..n).map(|_| self.random_generator.generate()).collect();
        let random_part = Polynomial::new(&bs);
        target + blinder * random_part
    }

    fn round_1(
        &self,
        witness: &Witness<F>,
        common_preprocessed_input: &CommonPreprocessedInput<F>,
    ) -> Round1Result<F, CS::Commitment>
    where
        CS: IsPolynomialCommitmentScheme<F, Polynomial = Polynomial<FieldElement<F>>>,
    {
        let p_a = Polynomial::interpolate_fft::<F>(&witness.a)
            .expect("xs and ys have equal length and xs are unique");
        let p_b = Polynomial::interpolate_fft::<F>(&witness.b)
            .expect("xs and ys have equal length and xs are unique");
        let p_c = Polynomial::interpolate_fft::<F>(&witness.c)
            .expect("xs and ys have equal length and xs are unique");

        let z_h = Polynomial::new_monomial(FieldElement::one(), common_preprocessed_input.n)
            - FieldElement::<F>::one();
        let p_a = self.blind_polynomial(&p_a, &z_h, 2);
        let p_b = self.blind_polynomial(&p_b, &z_h, 2);
        let p_c = self.blind_polynomial(&p_c, &z_h, 2);

        let a_1 = self.commitment_scheme.commit(&p_a);
        let b_1 = self.commitment_scheme.commit(&p_b);
        let c_1 = self.commitment_scheme.commit(&p_c);

        Round1Result {
            a_1,
            b_1,
            c_1,
            p_a,
            p_b,
            p_c,
        }
    }

    fn round_2(
        &self,
        witness: &Witness<F>,
        common_preprocessed_input: &CommonPreprocessedInput<F>,
        beta: FieldElement<F>,
        gamma: FieldElement<F>,
    ) -> Round2Result<F, CS::Commitment>
    where
        CS: IsPolynomialCommitmentScheme<F, Polynomial = Polynomial<FieldElement<F>>>,
    {
        let cpi = common_preprocessed_input;
        let mut coefficients: Vec<FieldElement<F>> = vec![FieldElement::one()];
        let (s1, s2, s3) = (&cpi.s1_lagrange, &cpi.s2_lagrange, &cpi.s3_lagrange);

        let k2 = &cpi.k1 * &cpi.k1;

        let lp = |w: &FieldElement<F>, eta: &FieldElement<F>| w + &beta * eta + &gamma;

        for i in 0..&cpi.n - 1 {
            let (a_i, b_i, c_i) = (&witness.a[i], &witness.b[i], &witness.c[i]);
            let num = lp(a_i, &cpi.domain[i])
                * lp(b_i, &(&cpi.domain[i] * &cpi.k1))
                * lp(c_i, &(&cpi.domain[i] * &k2));
            let den = lp(a_i, &s1[i]) * lp(b_i, &s2[i]) * lp(c_i, &s3[i]);
            let new_factor = num / den;
            let new_term = coefficients.last().unwrap() * &new_factor;
            coefficients.push(new_term);
        }

        let p_z = Polynomial::interpolate_fft::<F>(&coefficients)
            .expect("xs and ys have equal length and xs are unique");
        let z_h = Polynomial::new_monomial(FieldElement::one(), common_preprocessed_input.n)
            - FieldElement::<F>::one();
        let p_z = self.blind_polynomial(&p_z, &z_h, 3);
        let z_1 = self.commitment_scheme.commit(&p_z);
        Round2Result {
            z_1,
            p_z,
            beta,
            gamma,
        }
    }

    fn round_3(
        &self,
        common_preprocessed_input: &CommonPreprocessedInput<F>,
        public_input: &[FieldElement<F>],
        Round1Result { p_a, p_b, p_c, .. }: &Round1Result<F, CS::Commitment>,
        Round2Result {
            p_z, beta, gamma, ..
        }: &Round2Result<F, CS::Commitment>,
        alpha: FieldElement<F>,
    ) -> Round3Result<F, CS::Commitment>
    where
        CS: IsPolynomialCommitmentScheme<F, Polynomial = Polynomial<FieldElement<F>>>,
    {
        let cpi = common_preprocessed_input;
        let k2 = &cpi.k1 * &cpi.k1;

        let one = Polynomial::new_monomial(FieldElement::one(), 0);
        let p_x = &Polynomial::new_monomial(FieldElement::<F>::one(), 1);
        let zh = Polynomial::new_monomial(FieldElement::<F>::one(), cpi.n) - &one;

        let z_x_omega_coefficients: Vec<FieldElement<F>> = p_z
            .coefficients()
            .iter()
            .enumerate()
            .map(|(i, x)| x * &cpi.domain[i % cpi.n])
            .collect();
        let z_x_omega = Polynomial::new(&z_x_omega_coefficients);
        let mut e1 = vec![FieldElement::<F>::zero(); cpi.domain.len()];
        e1[0] = FieldElement::one();
        let l1 = Polynomial::interpolate_fft::<F>(&e1)
            .expect("xs and ys have equal length and xs are unique");
        let mut p_pi_y = public_input.to_vec();
        p_pi_y.append(&mut vec![FieldElement::zero(); cpi.n - public_input.len()]);
        let p_pi = Polynomial::interpolate_fft::<F>(&p_pi_y)
            .expect("xs and ys have equal length and xs are unique");

        // Compute p
        // To leverage FFT we work with the evaluation form of every polynomial
        // involved
        // TODO: check a factor of 4 is a sensible upper bound
        let degree = 4 * cpi.n;
        let offset = &cpi.k1;
        let p_a_eval = Polynomial::evaluate_offset_fft(p_a, 1, Some(degree), offset).unwrap();
        let p_b_eval = Polynomial::evaluate_offset_fft(p_b, 1, Some(degree), offset).unwrap();
        let p_c_eval = Polynomial::evaluate_offset_fft(p_c, 1, Some(degree), offset).unwrap();
        let ql_eval = Polynomial::evaluate_offset_fft(&cpi.ql, 1, Some(degree), offset).unwrap();
        let qr_eval = Polynomial::evaluate_offset_fft(&cpi.qr, 1, Some(degree), offset).unwrap();
        let qm_eval = Polynomial::evaluate_offset_fft(&cpi.qm, 1, Some(degree), offset).unwrap();
        let qo_eval = Polynomial::evaluate_offset_fft(&cpi.qo, 1, Some(degree), offset).unwrap();
        let qc_eval = Polynomial::evaluate_offset_fft(&cpi.qc, 1, Some(degree), offset).unwrap();
        let p_pi_eval = Polynomial::evaluate_offset_fft(&p_pi, 1, Some(degree), offset).unwrap();
        let p_x_eval = Polynomial::evaluate_offset_fft(p_x, 1, Some(degree), offset).unwrap();
        let p_z_eval = Polynomial::evaluate_offset_fft(p_z, 1, Some(degree), offset).unwrap();
        let p_z_x_omega_eval =
            Polynomial::evaluate_offset_fft(&z_x_omega, 1, Some(degree), offset).unwrap();
        let p_s1_eval = Polynomial::evaluate_offset_fft(&cpi.s1, 1, Some(degree), offset).unwrap();
        let p_s2_eval = Polynomial::evaluate_offset_fft(&cpi.s2, 1, Some(degree), offset).unwrap();
        let p_s3_eval = Polynomial::evaluate_offset_fft(&cpi.s3, 1, Some(degree), offset).unwrap();
        let l1_eval = Polynomial::evaluate_offset_fft(&l1, 1, Some(degree), offset).unwrap();

        let p_constraints_eval: Vec<_> = p_a_eval
            .iter()
            .zip(p_b_eval.iter())
            .zip(p_c_eval.iter())
            .zip(ql_eval.iter())
            .zip(qr_eval.iter())
            .zip(qm_eval.iter())
            .zip(qo_eval.iter())
            .zip(qc_eval.iter())
            .zip(p_pi_eval.iter())
            .map(|((((((((a, b), c), ql), qr), qm), qo), qc), pi)| {
                a * b * qm + a * ql + b * qr + c * qo + qc + pi
            })
            .collect();

        let f_eval: Vec<_> = p_a_eval
            .iter()
            .zip(p_b_eval.iter())
            .zip(p_c_eval.iter())
            .zip(p_x_eval.iter())
            .map(|(((a, b), c), x)| {
                (a + x * beta + gamma)
                    * (b + x * beta * &cpi.k1 + gamma)
                    * (c + x * beta * &k2 + gamma)
            })
            .collect();

        let g_eval: Vec<_> = p_a_eval
            .iter()
            .zip(p_b_eval.iter())
            .zip(p_c_eval.iter())
            .zip(p_s1_eval.iter())
            .zip(p_s2_eval.iter())
            .zip(p_s3_eval.iter())
            .map(|(((((a, b), c), s1), s2), s3)| {
                (a + s1 * beta + gamma) * (b + s2 * beta + gamma) * (c + s3 * beta + gamma)
            })
            .collect();

        let p_permutation_1_eval: Vec<_> = g_eval
            .iter()
            .zip(f_eval.iter())
            .zip(p_z_eval.iter())
            .zip(p_z_x_omega_eval.iter())
            .map(|(((g, f), z), y)| g * y - f * z)
            .collect();

        let p_permutation_2_eval: Vec<_> = p_z_eval
            .iter()
            .zip(l1_eval.iter())
            .map(|(z, l)| (z - FieldElement::<F>::one()) * l)
            .collect();

        let p_eval: Vec<_> = p_permutation_2_eval
            .iter()
            .zip(p_permutation_1_eval.iter())
            .zip(p_constraints_eval.iter())
            .map(|((p2, p1), co)| (p2 * &alpha + p1) * &alpha + co)
            .collect();

        let mut zh_eval = Polynomial::evaluate_offset_fft(&zh, 1, Some(degree), offset).unwrap();
        FieldElement::inplace_batch_inverse(&mut zh_eval).unwrap();
        let c: Vec<_> = p_eval
            .iter()
            .zip(zh_eval.iter())
            .map(|(a, b)| a * b)
            .collect();
        let mut t = Polynomial::interpolate_offset_fft(&c, offset).unwrap();

        polynomial::pad_with_zero_coefficients_to_length(&mut t, 3 * (&cpi.n + 2));
        let p_t_lo = Polynomial::new(&t.coefficients[..&cpi.n + 2]);
        let p_t_mid = Polynomial::new(&t.coefficients[&cpi.n + 2..2 * (&cpi.n + 2)]);
        let p_t_hi = Polynomial::new(&t.coefficients[2 * (&cpi.n + 2)..3 * (&cpi.n + 2)]);

        let b_0 = self.random_generator.generate();
        let b_1 = self.random_generator.generate();

        let p_t_lo = &p_t_lo + &b_0 * Polynomial::new_monomial(FieldElement::one(), cpi.n + 2);
        let p_t_mid =
            &p_t_mid - b_0 + &b_1 * Polynomial::new_monomial(FieldElement::one(), cpi.n + 2);
        let p_t_hi = &p_t_hi - b_1;

        let t_lo_1 = self.commitment_scheme.commit(&p_t_lo);
        let t_mid_1 = self.commitment_scheme.commit(&p_t_mid);
        let t_hi_1 = self.commitment_scheme.commit(&p_t_hi);

        Round3Result {
            t_lo_1,
            t_mid_1,
            t_hi_1,
            p_t_lo,
            p_t_mid,
            p_t_hi,
            alpha,
        }
    }

    fn round_4(
        &self,
        CommonPreprocessedInput { s1, s2, omega, .. }: &CommonPreprocessedInput<F>,
        Round1Result { p_a, p_b, p_c, .. }: &Round1Result<F, CS::Commitment>,
        Round2Result { p_z, .. }: &Round2Result<F, CS::Commitment>,
        zeta: FieldElement<F>,
    ) -> Round4Result<F> {
        let a_zeta = p_a.evaluate(&zeta);
        let b_zeta = p_b.evaluate(&zeta);
        let c_zeta = p_c.evaluate(&zeta);
        let s1_zeta = s1.evaluate(&zeta);
        let s2_zeta = s2.evaluate(&zeta);
        let z_zeta_omega = p_z.evaluate(&(&zeta * omega));
        Round4Result {
            a_zeta,
            b_zeta,
            c_zeta,
            s1_zeta,
            s2_zeta,
            z_zeta_omega,
            zeta,
        }
    }

    fn round_5(
        &self,
        common_preprocessed_input: &CommonPreprocessedInput<F>,
        round_1: &Round1Result<F, CS::Commitment>,
        round_2: &Round2Result<F, CS::Commitment>,
        round_3: &Round3Result<F, CS::Commitment>,
        round_4: &Round4Result<F>,
        transcript: &mut impl Transcript,
    ) -> Round5Result<F, CS::Proof>
    where
        CS: IsPolynomialCommitmentScheme<F, Polynomial = Polynomial<FieldElement<F>>>,
        FieldElement<F>: Borrow<<CS as IsPolynomialCommitmentScheme<F>>::Point>,
    {
        let cpi = common_preprocessed_input;
        let (r1, r2, r3, r4) = (round_1, round_2, round_3, round_4);
        // Precompute variables
        let k2 = &cpi.k1 * &cpi.k1;
        let zeta_raised_n = Polynomial::new_monomial(r4.zeta.pow(cpi.n + 2), 0); // TODO: Paper says n and 2n, but Gnark uses n+2 and 2n+4
        let zeta_raised_2n = Polynomial::new_monomial(r4.zeta.pow(2 * cpi.n + 4), 0);

        let l1_zeta = (&r4.zeta.pow(cpi.n as u64) - FieldElement::<F>::one())
            / (&r4.zeta - FieldElement::<F>::one())
            / FieldElement::<F>::from(cpi.n as u64);

        let mut p_non_constant = &cpi.qm * &r4.a_zeta * &r4.b_zeta
            + &r4.a_zeta * &cpi.ql
            + &r4.b_zeta * &cpi.qr
            + &r4.c_zeta * &cpi.qo
            + &cpi.qc;

        let r_2_1 = (&r4.a_zeta + &r2.beta * &r4.zeta + &r2.gamma)
            * (&r4.b_zeta + &r2.beta * &cpi.k1 * &r4.zeta + &r2.gamma)
            * (&r4.c_zeta + &r2.beta * &k2 * &r4.zeta + &r2.gamma)
            * &r2.p_z;
        let r_2_2 = (&r4.a_zeta + &r2.beta * &r4.s1_zeta + &r2.gamma)
            * (&r4.b_zeta + &r2.beta * &r4.s2_zeta + &r2.gamma)
            * &r2.beta
            * &r4.z_zeta_omega
            * &cpi.s3;
        p_non_constant = p_non_constant + (r_2_2 - r_2_1) * &r3.alpha;

        let r_3 = &r2.p_z * l1_zeta;
        p_non_constant = p_non_constant + (r_3 * &r3.alpha * &r3.alpha);

        let partial_t = &r3.p_t_lo + zeta_raised_n * &r3.p_t_mid + zeta_raised_2n * &r3.p_t_hi;

        // TODO: Refactor to remove clones.
        let polynomials = vec![
            partial_t,
            p_non_constant,
            r1.p_a.clone(),
            r1.p_b.clone(),
            r1.p_c.clone(),
            cpi.s1.clone(),
            cpi.s2.clone(),
        ];
        let ys: Vec<FieldElement<F>> = polynomials.iter().map(|p| p.evaluate(&r4.zeta)).collect();
        let w_zeta_1 = self
            .commitment_scheme
            //TODO (pat_stiles): Eliminate this clone
            .open_batch(r4.zeta.clone(), &ys, &polynomials, Some(transcript));

        let w_zeta_omega_1 = self
            .commitment_scheme
            //TODO (pat_stiles): Eliminate this clone
            .open(
                r4.zeta.clone() * cpi.omega.clone(),
                &r4.z_zeta_omega,
                &r2.p_z,
                Some(transcript),
            );

        Round5Result {
            w_zeta_1,
            w_zeta_omega_1,
            p_non_constant_zeta: ys[1].clone(),
            t_zeta: ys[0].clone(),
        }
    }

    pub fn prove(
        &self,
        witness: &Witness<F>,
        public_input: &[FieldElement<F>],
        common_preprocessed_input: &CommonPreprocessedInput<F>,
        vk: &VerificationKey<CS::Commitment>,
    ) -> Proof<F, CS>
    where
        CS: IsPolynomialCommitmentScheme<F, Polynomial = Polynomial<FieldElement<F>>>,
        FieldElement<F>: Borrow<<CS as IsPolynomialCommitmentScheme<F>>::Point>,
    {
        let mut transcript = new_strong_fiat_shamir_transcript::<F, CS>(vk, public_input);

        // Round 1
        let round_1 = self.round_1(witness, common_preprocessed_input);
        transcript.append(&round_1.a_1.as_bytes());
        transcript.append(&round_1.b_1.as_bytes());
        transcript.append(&round_1.c_1.as_bytes());

        // Round 2
        // TODO: Handle error
        let beta = FieldElement::from_bytes_be(&transcript.challenge()).unwrap();
        let gamma = FieldElement::from_bytes_be(&transcript.challenge()).unwrap();

        let round_2 = self.round_2(witness, common_preprocessed_input, beta, gamma);
        transcript.append(&round_2.z_1.as_bytes());

        // Round 3
        let alpha = FieldElement::from_bytes_be(&transcript.challenge()).unwrap();
        let round_3 = self.round_3(
            common_preprocessed_input,
            public_input,
            &round_1,
            &round_2,
            alpha,
        );
        transcript.append(&round_3.t_lo_1.as_bytes());
        transcript.append(&round_3.t_mid_1.as_bytes());
        transcript.append(&round_3.t_hi_1.as_bytes());

        // Round 4
        let zeta = FieldElement::from_bytes_be(&transcript.challenge()).unwrap();
        let round_4 = self.round_4(common_preprocessed_input, &round_1, &round_2, zeta);

        transcript.append(&round_4.a_zeta.to_bytes_be());
        transcript.append(&round_4.b_zeta.to_bytes_be());
        transcript.append(&round_4.c_zeta.to_bytes_be());
        transcript.append(&round_4.s1_zeta.to_bytes_be());
        transcript.append(&round_4.s2_zeta.to_bytes_be());
        transcript.append(&round_4.z_zeta_omega.to_bytes_be());

        // Round 5
        let round_5 = self.round_5(
            common_preprocessed_input,
            &round_1,
            &round_2,
            &round_3,
            &round_4,
            &mut transcript,
        );

        Proof {
            a_1: round_1.a_1,
            b_1: round_1.b_1,
            c_1: round_1.c_1,
            z_1: round_2.z_1,
            t_lo_1: round_3.t_lo_1,
            t_mid_1: round_3.t_mid_1,
            t_hi_1: round_3.t_hi_1,
            a_zeta: round_4.a_zeta,
            b_zeta: round_4.b_zeta,
            c_zeta: round_4.c_zeta,
            s1_zeta: round_4.s1_zeta,
            s2_zeta: round_4.s2_zeta,
            z_zeta_omega: round_4.z_zeta_omega,
            w_zeta_1: round_5.w_zeta_1,
            w_zeta_omega_1: round_5.w_zeta_omega_1,
            p_non_constant_zeta: round_5.p_non_constant_zeta,
            t_zeta: round_5.t_zeta,
        }
    }
}

#[cfg(test)]
mod tests {
    use lambdaworks_crypto::fiat_shamir::default_transcript::DefaultTranscript;
    use lambdaworks_math::{
        cyclic_group::IsGroup,
        elliptic_curve::{
            short_weierstrass::{
                curves::bls12_381::{curve::BLS12381Curve, default_types::FrElement},
                point::ShortWeierstrassProjectivePoint,
            },
            traits::IsEllipticCurve,
        },
    };

    use crate::{
        test_utils::circuit_1::{test_common_preprocessed_input_1, test_witness_1},
        test_utils::utils::{test_srs, FpElement, TestRandomFieldGenerator, KZG},
    };

    use super::*;

    fn alpha() -> FrElement {
        FrElement::from_hex_unchecked(
            "583cfb0df2ef98f2131d717bc6aadd571c5302597c135cab7c00435817bf6e50",
        )
    }

    fn beta() -> FrElement {
        FrElement::from_hex_unchecked(
            "bdda7414bdf5bf42b77cbb3af4a82f32ec7622dd6c71575bede021e6e4609d4",
        )
    }

    fn gamma() -> FrElement {
        FrElement::from_hex_unchecked(
            "58f6690d9b36e62e4a0aef27612819288df2a3ff5bf01597cf06779503f51583",
        )
    }

    fn zeta() -> FrElement {
        FrElement::from_hex_unchecked(
            "2a4040abb941ee5e2a42602a7a60d282a430a4cf099fa3bb0ba8f4da628ec59a",
        )
    }

    fn upsilon() -> FrElement {
        FrElement::from_hex_unchecked(
            "2d15959489a2a8e44693221ca7cbdcab15253d6bae9fd7fe0664cff02fe4f1cf",
        )
    }

    #[test]
    fn test_round_1() {
        let witness = test_witness_1(FrElement::from(2), FrElement::from(2));
        let common_preprocessed_input = test_common_preprocessed_input_1();
        let srs = test_srs(common_preprocessed_input.n);
        let kzg = KZG::new(srs);
        let random_generator = TestRandomFieldGenerator {};

        let prover = Prover::new(kzg, random_generator);
        let round_1 = prover.round_1(&witness, &common_preprocessed_input);
        let a_1_expected = BLS12381Curve::create_point_from_affine(
            FpElement::from_hex_unchecked("17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb"),
            FpElement::from_hex_unchecked("114d1d6855d545a8aa7d76c8cf2e21f267816aef1db507c96655b9d5caac42364e6f38ba0ecb751bad54dcd6b939c2ca"),
        ).unwrap();
        let b_1_expected = BLS12381Curve::create_point_from_affine(
            FpElement::from_hex_unchecked("44ed7c3ed015c6a39c350cd06d03b48d3e1f5eaf7a256c5b6203886e6e78cd9b76623d163da4dfb0f2491e7cc06408"),
            FpElement::from_hex_unchecked("14c4464d2556fdfdc8e31068ef8d953608e511569a236c825f2ddab4fe04af03aba29e38b9b2b6221243124d235f4c67"),
        ).unwrap();
        let c_1_expected = BLS12381Curve::create_point_from_affine(
            FpElement::from_hex_unchecked("7726dc031bd26122395153ca428d5e6dea0a64c1f9b3b1bb2f2508a5eb6ea0ea0363294fad3160858bc87e46d3422fd"),
            FpElement::from_hex_unchecked("8db0c15bfd77df7fe66284c3b04e6043eaba99ef6a845d4f7255fd0da95f2fb8e474df2e7f8e1a38829f7a9612a9b87"),
        ).unwrap();
        assert_eq!(round_1.a_1, a_1_expected);
        assert_eq!(round_1.b_1, b_1_expected);
        assert_eq!(round_1.c_1, c_1_expected);
    }

    #[test]
    fn test_round_2() {
        let witness = test_witness_1(FrElement::from(2), FrElement::from(2));
        let common_preprocessed_input = test_common_preprocessed_input_1();
        let srs = test_srs(common_preprocessed_input.n);
        let kzg = KZG::new(srs);
        let random_generator = TestRandomFieldGenerator {};
        let prover = Prover::new(kzg, random_generator);

        let result_2 = prover.round_2(&witness, &common_preprocessed_input, beta(), gamma());
        let z_1_expected = BLS12381Curve::create_point_from_affine(
            FpElement::from_hex_unchecked("3e8322968c3496cf1b5786d4d71d158a646ec90c14edf04e758038e1f88dcdfe8443fcecbb75f3074a872a380391742"),
            FpElement::from_hex_unchecked("11eac40d09796ff150004e7b858d83ddd9fe995dced0b3fbd7535d6e361729b25d488799da61fdf1d7b5022684053327"),
        ).unwrap();
        assert_eq!(result_2.z_1, z_1_expected);
    }

    #[test]
    fn test_round_3() {
        let witness = test_witness_1(FrElement::from(2), FrElement::from(2));
        let common_preprocessed_input = test_common_preprocessed_input_1();
        let srs = test_srs(common_preprocessed_input.n);
        let kzg = KZG::new(srs);
        let public_input = vec![FieldElement::from(2_u64), FieldElement::from(4)];
        let random_generator = TestRandomFieldGenerator {};
        let prover = Prover::new(kzg, random_generator);
        let round_1 = prover.round_1(&witness, &common_preprocessed_input);
        let round_2 = prover.round_2(&witness, &common_preprocessed_input, beta(), gamma());
        let round_3 = prover.round_3(
            &common_preprocessed_input,
            &public_input,
            &round_1,
            &round_2,
            alpha(),
        );

        let t_lo_1_expected = BLS12381Curve::create_point_from_affine(
            FpElement::from_hex_unchecked("9f511a769e77e87537b0749d65f467532fbf0f9dc1bcc912c333741be9d0a613f61e5fe595996964646ce30794701e5"),
            FpElement::from_hex_unchecked("89fd6bb571323912210517237d6121144fc01ba2756f47c12c9cc94fc9197313867d68530f152dc8d447f10fcf75a6c"),
        ).unwrap();
        let t_mid_1_expected = BLS12381Curve::create_point_from_affine(
            FpElement::from_hex_unchecked("f96d8a93f3f5be2ab2819891f41c9f883cacea63da423e6ed1701765fcd659fc11e056a48c554f5df3a9c6603d48ca8"),
            FpElement::from_hex_unchecked("14fa74fa049b7276007b739f3b8cfeac09e8cfabd4f858b6b99798c81124c34851960bebda90133cb03c981c08c8b6d3"),
        ).unwrap();
        let t_hi_1_expected = ShortWeierstrassProjectivePoint::<BLS12381Curve>::neutral_element();

        assert_eq!(round_3.t_lo_1, t_lo_1_expected);
        assert_eq!(round_3.t_mid_1, t_mid_1_expected);
        assert_eq!(round_3.t_hi_1, t_hi_1_expected);
    }

    #[test]
    fn test_round_4() {
        let witness = test_witness_1(FrElement::from(2), FrElement::from(2));
        let common_preprocessed_input = test_common_preprocessed_input_1();
        let srs = test_srs(common_preprocessed_input.n);
        let kzg = KZG::new(srs);
        let random_generator = TestRandomFieldGenerator {};
        let prover = Prover::new(kzg, random_generator);

        let round_1 = prover.round_1(&witness, &common_preprocessed_input);
        let round_2 = prover.round_2(&witness, &common_preprocessed_input, beta(), gamma());

        let round_4 = prover.round_4(&common_preprocessed_input, &round_1, &round_2, zeta());
        let expected_a_value = FrElement::from_hex_unchecked(
            "2c090a95b57f1f493b7b747bba34fef7772fd72f97d718ed69549641a823eb2e",
        );
        let expected_b_value = FrElement::from_hex_unchecked(
            "5975959d91369ba4e7a03c6ae94b7fe98e8b61b7bf9af63c8ae0759e17ac0c7e",
        );
        let expected_c_value = FrElement::from_hex_unchecked(
            "6bf31edeb4344b7d2df2cb1bd40b4d13e182d9cb09f89591fa043c1a34b4a93",
        );
        let expected_z_value = FrElement::from_hex_unchecked(
            "38e2ec8e7c3dab29e2b8e9c8ea152914b8fe4612e91f2902c80238efcf21f4ee",
        );
        let expected_s1_value = FrElement::from_hex_unchecked(
            "472f66db4fb6947d9ed9808241fe82324bc08aa2a54be93179db8e564e1137d4",
        );
        let expected_s2_value = FrElement::from_hex_unchecked(
            "5588f1239c24efe0538868d0f716984e69c6980e586864f615e4b0621fdc6f81",
        );

        assert_eq!(round_4.a_zeta, expected_a_value);
        assert_eq!(round_4.b_zeta, expected_b_value);
        assert_eq!(round_4.c_zeta, expected_c_value);
        assert_eq!(round_4.z_zeta_omega, expected_z_value);
        assert_eq!(round_4.s1_zeta, expected_s1_value);
        assert_eq!(round_4.s2_zeta, expected_s2_value);
    }

    #[test]
    #[ignore]
    fn test_round_5() {
        let witness = test_witness_1(FrElement::from(2), FrElement::from(2));
        let common_preprocessed_input = test_common_preprocessed_input_1();
        let srs = test_srs(common_preprocessed_input.n);
        let kzg = KZG::new(srs);
        let public_input = vec![FieldElement::from(2_u64), FieldElement::from(4)];
        let random_generator = TestRandomFieldGenerator {};
        let prover = Prover::new(kzg, random_generator);

        let round_1 = prover.round_1(&witness, &common_preprocessed_input);
        let round_2 = prover.round_2(&witness, &common_preprocessed_input, beta(), gamma());

        let round_3 = prover.round_3(
            &common_preprocessed_input,
            &public_input,
            &round_1,
            &round_2,
            alpha(),
        );

        let round_4 = prover.round_4(&common_preprocessed_input, &round_1, &round_2, zeta());

        let expected_w_zeta_1 = BLS12381Curve::create_point_from_affine(
            FpElement::from_hex_unchecked("fa6250b80a418f0548b132ac264ff9915b2076c0c2548da9316ae19ffa35bbcf905d9f02f9274739608045ef83a4757"),
            FpElement::from_hex_unchecked("17713ade2dbd66e923d4092a5d2da98202959dd65a15e9f7791fab3c0dd08788aa9b4a1cb21d04e0c43bd29225472145"),
        ).unwrap();
        let expected_w_zeta_omega_1 = BLS12381Curve::create_point_from_affine(
            FpElement::from_hex_unchecked("4484f08f8eaccf28bab8ee9539e6e7f4059cb1ce77b9b18e9e452f387163dc0b845f4874bf6445399e650d362799ff5"),
            FpElement::from_hex_unchecked("1254347a0fa2ac856917825a5cff5f9583d39a52edbc2be5bb10fabd0c04d23019bcb963404345743120310fd734a61a"),
        ).unwrap();

        let mut transcript = DefaultTranscript::new();
        let round_5 = prover.round_5(
            &common_preprocessed_input,
            &round_1,
            &round_2,
            &round_3,
            &round_4,
            &mut transcript,
        );
        assert_eq!(round_5.w_zeta_1, expected_w_zeta_1);
        assert_eq!(round_5.w_zeta_omega_1, expected_w_zeta_omega_1);
    }
}
