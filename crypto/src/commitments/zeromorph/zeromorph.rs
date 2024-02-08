#![allow(clippy::too_many_arguments)]
#![allow(clippy::type_complexity)]

use crate::{
    commitments::traits::IsPolynomialCommitmentScheme, fiat_shamir::transcript::Transcript,
};
use core::mem;
use lambdaworks_math::{
    cyclic_group::IsGroup,
    elliptic_curve::traits::IsPairing,
    errors::DeserializationError,
    field::{
        element::FieldElement,
        traits::{IsField, IsPrimeField},
    },
    msm::pippenger::msm,
    polynomial::{dense_multilinear_poly::DenseMultilinearPolynomial, Polynomial},
    traits::{AsBytes, ByteConversion, Deserializable},
    unsigned_integer::element::UnsignedInteger,
};
use std::{borrow::Borrow, iter, marker::PhantomData};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use super::structs::{ZeromorphProof, ZeromorphProverKey, ZeromorphVerifierKey};

fn compute_multilinear_quotients<P: IsPairing>(
    poly: &DenseMultilinearPolynomial<P::BaseField>,
    u_challenge: &[FieldElement<P::BaseField>],
) -> (
    Vec<Polynomial<FieldElement<P::BaseField>>>,
    FieldElement<P::BaseField>,
)
where
    <<P as IsPairing>::BaseField as IsField>::BaseType: Send + Sync,
{
    assert_eq!(poly.num_vars(), u_challenge.len());

    let mut g = poly.evals().to_vec();
    let mut quotients: Vec<_> = u_challenge
        .iter()
        .enumerate()
        .map(|(i, x_i)| {
            let (g_lo, g_hi) = g.split_at_mut(1 << (poly.num_vars() - 1 - i));
            let mut quotient = vec![FieldElement::zero(); g_lo.len()];

            #[cfg(feature = "parallel")]
            let quotient_iter = quotient.par_iter_mut();

            #[cfg(not(feature = "parallel"))]
            let quotient_iter = quotient.iter_mut();

            quotient_iter
                .zip(&*g_lo)
                .zip(&*g_hi)
                .for_each(|((q, g_lo), g_hi)| {
                    *q = g_hi - g_lo;
                });

            #[cfg(feature = "parallel")]
            let g_lo_iter = g_lo.par_iter_mut();

            #[cfg(not(feature = "parallel"))]
            let g_lo_iter = g_lo.iter_mut();
            g_lo_iter.zip(g_hi).for_each(|(g_lo, g_hi)| {
                *g_lo += (&*g_hi - g_lo as &_) * x_i;
            });

            g.truncate(1 << (poly.num_vars() - 1 - i));

            Polynomial::new(&quotient)
        })
        .collect();
    quotients.reverse();
    (quotients, g[0].clone())
}

fn compute_batched_lifted_degree_quotient<P: IsPairing>(
    n: usize,
    quotients: &Vec<Polynomial<FieldElement<P::BaseField>>>,
    y_challenge: &FieldElement<P::BaseField>,
) -> Polynomial<FieldElement<P::BaseField>> {
    // Compute \hat{q} = \sum_k y^k * X^{N - d_k - 1} * q_k
    let mut scalar = FieldElement::<P::BaseField>::one(); // y^k
                                                          // Rather than explicitly computing the shifts of q_k by N - d_k - 1 (i.e. multiplying q_k by X^{N - d_k - 1})
                                                          // then accumulating them, we simply accumulate y^k*q_k into \hat{q} at the index offset N - d_k - 1
    let q_hat =
        quotients
            .iter()
            .enumerate()
            .fold(vec![FieldElement::zero(); n], |mut q_hat, (idx, q)| {
                #[cfg(feature = "parallel")]
                let q_hat_iter = q_hat[n - (1 << idx)..].par_iter_mut();

                #[cfg(not(feature = "parallel"))]
                let q_hat_iter = q_hat[n - (1 << idx)..].iter_mut();
                q_hat_iter.zip(&q.coefficients).for_each(|(q_hat, q)| {
                    *q_hat += &scalar * q;
                });
                scalar *= y_challenge;
                q_hat
            });

    Polynomial::new(&q_hat)
}

fn eval_and_quotient_scalars<P: IsPairing>(
    y_challenge: &FieldElement<P::BaseField>,
    x_challenge: &FieldElement<P::BaseField>,
    z_challenge: &FieldElement<P::BaseField>,
    challenges: &[FieldElement<P::BaseField>],
) -> (
    FieldElement<P::BaseField>,
    (
        Vec<FieldElement<P::BaseField>>,
        Vec<FieldElement<P::BaseField>>,
    ),
) {
    let num_vars = challenges.len();

    // squares of x = [x, x^2, .. x^{2^k}, .. x^{2^num_vars}]
    let squares_of_x: Vec<_> = iter::successors(Some(x_challenge.clone()), |x| Some(x.square()))
        .take(num_vars + 1)
        .collect();

    // offsets of x =
    let offsets_of_x = {
        let mut offsets_of_x = squares_of_x
            .iter()
            .rev()
            .skip(1)
            .scan(FieldElement::one(), |acc, pow_x| {
                *acc *= pow_x;
                Some(acc.clone())
            })
            .collect::<Vec<_>>();
        offsets_of_x.reverse();
        offsets_of_x
    };

    let vs = {
        let v_numer: FieldElement<P::BaseField> =
            squares_of_x[num_vars].clone() - FieldElement::one();
        let mut v_denoms = squares_of_x
            .iter()
            .map(|square_of_x| square_of_x - FieldElement::one())
            .collect::<Vec<_>>();
        //TODO: catch this unwrap()
        FieldElement::inplace_batch_inverse(&mut v_denoms).unwrap();
        v_denoms
            .iter()
            .map(|v_denom| &v_numer * v_denom)
            .collect::<Vec<_>>()
    };

    let q_scalars = iter::successors(Some(FieldElement::one()), |acc| Some(acc * y_challenge))
        .take(num_vars)
        .zip(offsets_of_x)
        .zip(squares_of_x)
        .zip(&vs)
        .zip(&vs[1..])
        .zip(challenges.iter().rev())
        .map(
            |(((((power_of_y, offset_of_x), square_of_x), v_i), v_j), u_i)| {
                (
                    -(power_of_y * offset_of_x),
                    -(z_challenge * (square_of_x * v_j - u_i * v_i)),
                )
            },
        )
        .unzip();

    // -vs[0] * z = -z * (x^(2^num_vars) - 1) / (x - 1) = -z Φ_n(x)
    (-vs[0].clone() * z_challenge, q_scalars)
}

pub struct Zeromorph<P: IsPairing> {
    pk: ZeromorphProverKey<P>,
    vk: ZeromorphVerifierKey<P>,
    _phantom: PhantomData<P>,
}

impl<'a, P: IsPairing> Zeromorph<P>
where
    <P as IsPairing>::G1Point: Deserializable,
    <P as IsPairing>::G2Point: Deserializable,
{
    //TODO: should we derive the pk and ck within this function directly from srs
    pub fn new(pk: ZeromorphProverKey<P>, vk: ZeromorphVerifierKey<P>) -> Self {
        Self {
            pk,
            vk,
            _phantom: PhantomData,
        }
    }

    /// Extracts pk and vk from binary file
    pub fn from_file(
        file_path: &str,
    ) -> Result<Self, crate::errors::ProverVerifyKeysFromFileError> {
        let bytes = std::fs::read(file_path)?;
        Ok(Self::deserialize(&bytes)?)
    }
}

impl<P: IsPairing> AsBytes for Zeromorph<P>
where
    <P as IsPairing>::G1Point: AsBytes,
    <P as IsPairing>::G2Point: AsBytes,
{
    fn as_bytes(&self) -> Vec<u8> {
        let mut serialized_data: Vec<u8> = Vec::new();
        // First 4 bytes encodes protocol version
        let protocol_version: [u8; 4] = [0; 4];

        serialized_data.extend(&protocol_version);

        // Second 8 bytes store the amount of G1 elements to be stored, this is more than can be indexed with a 64-bit architecture, and some millions of terabytes of data if the points were compressed
        let mut g1_powers_len_bytes: Vec<u8> = self.pk.g1_powers.len().to_le_bytes().to_vec();

        // For architectures with less than 64 bits for pointers
        // We add extra zeros at the end
        while g1_powers_len_bytes.len() < 8 {
            g1_powers_len_bytes.push(0)
        }

        serialized_data.extend(&g1_powers_len_bytes);

        // third 8 bytes store the amount of G1 offset elements to be stored, this is more than can be indexed with a 64-bit architecture, and some millions of terabytes of data if the points were compressed
        let mut offset_g1_powers_len_bytes: Vec<u8> =
            self.pk.offset_g1_powers.len().to_le_bytes().to_vec();

        // For architectures with less than 64 bits for pointers
        // We add extra zeros at the end
        while offset_g1_powers_len_bytes.len() < 8 {
            offset_g1_powers_len_bytes.push(0)
        }

        serialized_data.extend(&offset_g1_powers_len_bytes);

        // G1 elements
        for point in &self.pk.g1_powers {
            serialized_data.extend(point.as_bytes());
        }

        // G2 elements
        for point in &self.pk.offset_g1_powers {
            serialized_data.extend(point.as_bytes());
        }

        // NOTE: this could potentially be recycled from the g1_powers but being explicit reduces complexity
        serialized_data.extend(&self.vk.g1.as_bytes());
        serialized_data.extend(&self.vk.g2.as_bytes());
        serialized_data.extend(&self.vk.tau_g2.as_bytes());
        serialized_data.extend(&self.vk.tau_n_max_sub_2_n.as_bytes());

        serialized_data
    }
}

impl<P: IsPairing> Deserializable for Zeromorph<P>
where
    <P as IsPairing>::G1Point: Deserializable,
    <P as IsPairing>::G2Point: Deserializable,
{
    fn deserialize(bytes: &[u8]) -> Result<Self, DeserializationError> {
        const VERSION_OFFSET: usize = 4;
        const G1_LEN_OFFSET: usize = 12;
        const OFFSET_G1_LEN_OFFSET: usize = 20;

        let g1_powers_len_u64 = u64::from_le_bytes(
            // This unwrap can't fail since we are fixing the size of the slice
            bytes[VERSION_OFFSET..G1_LEN_OFFSET].try_into().unwrap(),
        );

        let g1_powers_len = usize::try_from(g1_powers_len_u64)
            .map_err(|_| DeserializationError::PointerSizeError)?;

        let offset_g1_powers_len_u64 = u64::from_le_bytes(
            // This unwrap can't fail since we are fixing the size of the slice
            bytes[G1_LEN_OFFSET..OFFSET_G1_LEN_OFFSET]
                .try_into()
                .unwrap(),
        );

        let offset_g1_powers_len = usize::try_from(offset_g1_powers_len_u64)
            .map_err(|_| DeserializationError::PointerSizeError)?;

        let mut g1_powers: Vec<P::G1Point> = Vec::new();
        let mut offset_g1_powers: Vec<P::G1Point> = Vec::new();

        let size_g1_point = mem::size_of::<P::G1Point>();
        let size_g2_point = mem::size_of::<P::G2Point>();

        for i in 0..g1_powers_len {
            // The second unwrap shouldn't fail since the amount of bytes is fixed
            let point = P::G1Point::deserialize(
                bytes[i * size_g1_point + OFFSET_G1_LEN_OFFSET
                    ..i * size_g1_point + size_g1_point + OFFSET_G1_LEN_OFFSET]
                    .try_into()
                    .unwrap(),
            )?;
            g1_powers.push(point);
        }

        let offset_g1_offset = size_g1_point * g1_powers_len + OFFSET_G1_LEN_OFFSET;
        for i in 0..g1_powers_len {
            // The second unwrap shouldn't fail since the amount of bytes is fixed
            let point = P::G1Point::deserialize(
                bytes[i * size_g1_point + offset_g1_offset
                    ..i * size_g1_point + size_g1_point + offset_g1_offset]
                    .try_into()
                    .unwrap(),
            )?;
            offset_g1_powers.push(point);
        }
        let pk = ZeromorphProverKey {
            g1_powers,
            offset_g1_powers,
        };

        let vk_offset = size_g1_point * g1_powers_len
            + size_g1_point * offset_g1_powers_len
            + OFFSET_G1_LEN_OFFSET;
        let g1 = P::G1Point::deserialize(
            bytes[vk_offset..size_g1_point + vk_offset]
                .try_into()
                .unwrap(),
        )?;
        let g2 = P::G2Point::deserialize(
            bytes[size_g2_point + vk_offset..size_g2_point + vk_offset]
                .try_into()
                .unwrap(),
        )?;
        let tau_g2 = P::G2Point::deserialize(
            bytes[2 * size_g2_point + vk_offset..2 * size_g2_point + vk_offset]
                .try_into()
                .unwrap(),
        )?;
        let tau_n_max_sub_2_n = P::G2Point::deserialize(
            bytes[3 * size_g2_point + vk_offset..3 * size_g2_point + vk_offset]
                .try_into()
                .unwrap(),
        )?;

        let vk = ZeromorphVerifierKey {
            g1,
            g2,
            tau_g2,
            tau_n_max_sub_2_n,
        };

        Ok(Zeromorph::new(pk, vk))
    }
}

impl<const N: usize, P: IsPairing> IsPolynomialCommitmentScheme<P::BaseField> for Zeromorph<P>
where
    <<P as IsPairing>::BaseField as IsField>::BaseType: Send + Sync,
    <P as IsPairing>::G1Point: AsBytes,
    <P as IsPairing>::G2Point: AsBytes,
    <P as IsPairing>::BaseField: IsPrimeField<RepresentativeType = UnsignedInteger<N>>,
    FieldElement<<P as IsPairing>::BaseField>: ByteConversion,
    //TODO: how to eliminate this complication in the trait interface
{
    type Polynomial = DenseMultilinearPolynomial<P::BaseField>;
    type Commitment = P::G1Point;
    type Proof = ZeromorphProof<P>;
    type Point = Vec<FieldElement<P::BaseField>>;

    // TODO: errors lengths are valid
    fn commit(&self, poly: &Self::Polynomial) -> Self::Commitment {
        let scalars: Vec<_> = poly
            .evals()
            .iter()
            .map(|eval| eval.representative())
            .collect();
        msm(&scalars, &self.pk.g1_powers[..poly.len()]).unwrap()
    }

    // TODO: errors lengths are valid
    fn open(
        &self,
        point: impl Borrow<Self::Point>,
        eval: &FieldElement<P::BaseField>,
        poly: &Self::Polynomial,
        transcript: Option<&mut dyn Transcript>,
    ) -> Self::Proof {
        let point = point.borrow();
        //TODO: error or interface or something
        let transcript = transcript.unwrap();
        let num_vars = point.len();
        let n: usize = 1 << num_vars;
        let mut pi_poly = Polynomial::new(&poly.evals());

        // Compute the multilinear quotients q_k = q_k(X_0, ..., X_{k-1})
        let (quotients, remainder) = compute_multilinear_quotients::<P>(&poly, &point);
        debug_assert_eq!(quotients.len(), poly.num_vars());
        debug_assert_eq!(remainder, *eval);

        // Compute and send commitments C_{q_k} = [q_k], k = 0, ..., d-1
        let q_k_commitments: Vec<_> = quotients
            .iter()
            .map(|q| {
                let q_k_commitment = {
                    //TODO: we only need to convert the offset scalars
                    let scalars: Vec<_> = q
                        .coefficients
                        .iter()
                        .map(|eval| eval.representative())
                        .collect();
                    msm(&scalars, &self.pk.g1_powers[..q.coeff_len()]).unwrap()
                };
                transcript.append(&q_k_commitment.as_bytes());
                q_k_commitment
            })
            .collect();

        // Get challenge y
        let y_challenge = FieldElement::from_bytes_be(&transcript.challenge()).unwrap();

        // Compute the batched, lifted-degree quotient \hat{q}
        let q_hat = compute_batched_lifted_degree_quotient::<P>(n, &quotients, &y_challenge);

        // Compute and send the commitment C_q = [\hat{q}]
        // commit at offset
        let offset = 1 << (quotients.len() - 1);

        // We perform an offset commitment here; This could be abstracted but its very small
        let q_hat_com = {
            //TODO: we only need to convert the offset scalars
            let scalars: Vec<_> = q_hat
                .coefficients
                .iter()
                .map(|eval| eval.representative())
                .collect();
            msm(
                &scalars[offset..],
                &self.pk.offset_g1_powers[offset..q_hat.coeff_len()],
            )
            .unwrap()
        };

        transcript.append(&q_hat_com.as_bytes());

        // Get challenges x and z
        let x_challenge = FieldElement::from_bytes_be(&transcript.challenge()).unwrap();
        let z_challenge = FieldElement::from_bytes_be(&transcript.challenge()).unwrap();

        let (eval_scalar, (zeta_degree_check_q_scalars, z_zmpoly_q_scalars)) =
            eval_and_quotient_scalars::<P>(&y_challenge, &x_challenge, &z_challenge, &point);
        // f = z * x * poly.Z + q_hat + (-z * x * Φ_n(x) * e) + x * ∑_k (q_scalars_k * q_k)
        pi_poly = pi_poly * &z_challenge;
        pi_poly = pi_poly + &q_hat;
        pi_poly[0] += eval * eval_scalar;
        quotients
            .into_iter()
            .zip(zeta_degree_check_q_scalars)
            .zip(z_zmpoly_q_scalars)
            .for_each(|((mut q, zeta_degree_check_q_scalar), z_zmpoly_q_scalar)| {
                q = q * &(zeta_degree_check_q_scalar + z_zmpoly_q_scalar);
                pi_poly = &pi_poly + &q;
            });

        debug_assert_eq!(
            pi_poly.evaluate(&x_challenge),
            FieldElement::<P::BaseField>::zero()
        );

        // Compute the KZG opening proof pi_poly; -> TODO should abstract into separate trait
        let pi = {
            pi_poly.ruffini_division_inplace(&x_challenge);
            let scalars: Vec<_> = pi_poly
                .coefficients
                .iter()
                .map(|eval| eval.representative())
                .collect();
            msm(&scalars, &self.pk.g1_powers[..pi_poly.coeff_len()]).unwrap()
        };

        ZeromorphProof {
            pi,
            q_hat_com,
            q_k_com: q_k_commitments,
        }
    }

    // TODO: errors lengths are valid
    fn open_batch(
        &self,
        point: impl Borrow<Self::Point>,
        evals: &[FieldElement<P::BaseField>],
        polys: &[Self::Polynomial],
        transcript: Option<&mut dyn Transcript>,
    ) -> Self::Proof {
        let transcript = transcript.unwrap();
        for (poly, eval) in polys.iter().zip(evals.iter()) {
            // Note by evaluating we confirm the number of challenges is valid
            debug_assert_eq!(poly.evaluate(point.borrow().clone()).unwrap(), *eval);
        }

        // Generate batching challenge \rho and powers 1,...,\rho^{m-1}
        let rho = FieldElement::from_bytes_be(&transcript.challenge()).unwrap();
        // Compute batching of unshifted polynomials f_i:
        let mut scalar = FieldElement::one();
        let (f_batched, batched_evaluation) = (0..polys.len()).fold(
            (
                DenseMultilinearPolynomial::new(vec![
                    FieldElement::zero();
                    1 << polys[0].num_vars()
                ]),
                FieldElement::zero(),
            ),
            |(mut f_batched, mut batched_evaluation), i| {
                f_batched = (f_batched + polys[i].clone() * &scalar).unwrap();
                batched_evaluation += &scalar * &evals[i];
                scalar *= &rho;
                (f_batched, batched_evaluation)
            },
        );
        Self::open(
            &self,
            point,
            &batched_evaluation,
            &f_batched,
            Some(transcript),
        )
    }

    // TODO: errors lengths are valid
    fn verify(
        &self,
        point: impl Borrow<Self::Point>,
        eval: &FieldElement<P::BaseField>,
        p_commitment: &Self::Commitment,
        proof: &Self::Proof,
        transcript: Option<&mut dyn Transcript>,
    ) -> bool {
        let ZeromorphProof {
            pi,
            q_k_com,
            q_hat_com,
        } = proof;

        let transcript = transcript.unwrap();
        let point = point.borrow();

        //Receive q_k commitments
        q_k_com
            .iter()
            .for_each(|c| transcript.append(&c.as_bytes()));

        // Challenge y
        let y_challenge = FieldElement::from_bytes_be(&transcript.challenge()).unwrap();

        // Receive commitment C_{q} -> Since our transcriptcript does not support appending and receiving data we instead store these commitments in a zeromorph proof struct
        transcript.append(&q_hat_com.as_bytes());

        // Challenge x, z
        let x_challenge = FieldElement::from_bytes_be(&transcript.challenge()).unwrap();
        let z_challenge = FieldElement::from_bytes_be(&transcript.challenge()).unwrap();

        let (eval_scalar, (mut q_scalars, zm_poly_q_scalars)) =
            eval_and_quotient_scalars::<P>(&y_challenge, &x_challenge, &z_challenge, &point);

        q_scalars
            .iter_mut()
            .zip(zm_poly_q_scalars)
            .for_each(|(scalar, zm_poly_scalar)| {
                *scalar += zm_poly_scalar;
            });

        let scalars = [
            vec![FieldElement::one(), z_challenge, (eval * eval_scalar)],
            q_scalars,
        ]
        .concat();

        let bases = [
            vec![q_hat_com.clone(), p_commitment.clone(), self.vk.g1.clone()],
            q_k_com.to_vec(),
        ]
        .concat();
        let zeta_z_com = {
            let scalars: Vec<_> = scalars
                .iter()
                .map(|scalar| scalar.representative())
                .collect();
            msm(&scalars, &bases).expect("`points` is sliced by `cs`'s length")
        };

        // e(pi, [tau]_2 - x * [1]_2) == e(C_{\zeta,Z}, [X^(N_max - 2^n - 1)]_2) <==> e(C_{\zeta,Z} - x * pi, [X^{N_max - 2^n - 1}]_2) * e(-pi, [tau_2]) == 1
        let e = P::compute_batch(&[
            (
                &pi,
                &self.vk.tau_g2.operate_with(
                    &self
                        .vk
                        .g2
                        .operate_with_self(x_challenge.representative())
                        .neg(),
                ),
            ),
            (&zeta_z_com, &self.vk.tau_n_max_sub_2_n),
        ])
        .unwrap();
        e == FieldElement::one()
    }

    fn verify_batch(
        &self,
        point: impl Borrow<Self::Point>,
        evals: &[FieldElement<P::BaseField>],
        p_commitments: &[Self::Commitment],
        proof: &Self::Proof,
        transcript: Option<&mut dyn Transcript>,
    ) -> bool {
        debug_assert_eq!(evals.len(), p_commitments.len());
        let transcript = transcript.unwrap();
        // Compute powers of batching challenge rho
        let rho = FieldElement::from_bytes_be(&transcript.challenge()).unwrap();

        // Compute batching of unshifted polynomials f_i:
        let mut scalar = FieldElement::<P::BaseField>::one();
        let (batched_eval, batched_commitment) = evals.iter().zip(p_commitments.iter()).fold(
            (FieldElement::zero(), P::G1Point::neutral_element()),
            |(mut batched_eval, batched_commitment), (eval, commitment)| {
                batched_eval += &scalar * eval;
                batched_commitment
                    .operate_with(&commitment.operate_with_self(scalar.representative()));
                scalar *= &rho;
                (batched_eval, batched_commitment)
            },
        );
        Self::verify(
            &self,
            point,
            &batched_eval,
            &batched_commitment,
            proof,
            Some(transcript),
        )
    }
}

#[cfg(test)]
mod test {

    use core::ops::Neg;

    use crate::{
        commitments::zeromorph::structs::ZeromorphSRS,
        fiat_shamir::default_transcript::DefaultTranscript,
    };

    use super::*;
    use lambdaworks_math::{
        elliptic_curve::short_weierstrass::curves::bls12_381::pairing::BLS12381AtePairing,
        polynomial::dense_multilinear_poly::log_2,
    };
    use rand_chacha::{
        rand_core::{RngCore, SeedableRng},
        ChaCha20Rng,
    };

    // Evaluate Phi_k(x) = \sum_{i=0}^k x^i using the direct inefficent formula
    fn phi<P: IsPairing>(
        challenge: &FieldElement<P::BaseField>,
        subscript: usize,
    ) -> FieldElement<P::BaseField> {
        let len = (1 << subscript) as u64;
        (0..len)
            .into_iter()
            .fold(FieldElement::zero(), |mut acc, i| {
                //Note this is ridiculous DevX
                acc += challenge.pow(i);
                acc
            })
    }

    fn rand_fr<P: IsPairing, R: RngCore>(rng: &mut R) -> FieldElement<P::BaseField>
    where
        FieldElement<<P as IsPairing>::BaseField>: ByteConversion,
    {
        let mut bytes = [0u8; 384];
        rng.fill_bytes(&mut bytes);
        FieldElement::<P::BaseField>::from_bytes_be(&bytes).unwrap()
    }

    #[test]
    fn prove_verify_single() {
        let max_vars = 16;
        let mut rng = &mut ChaCha20Rng::from_seed(*b"zeromorph_poly_commitment_scheme");
        let srs = ZeromorphSRS::setup(1 << (max_vars + 1), rng);

        for num_vars in 3..max_vars {
            // Setup
            let (pk, vk) = {
                let poly_size = 1 << (num_vars + 1);
                srs.trim(poly_size - 1).unwrap()
            };
            let zm = Zeromorph::<BLS12381AtePairing>::new(pk, vk);
            let poly = DenseMultilinearPolynomial::new(
                (0..(1 << num_vars))
                    .map(|_| rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng))
                    .collect::<Vec<_>>(),
            );
            let point = (0..num_vars)
                .map(|_| rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng))
                .collect::<Vec<_>>();
            let eval = poly.evaluate(point.clone()).unwrap();

            // Commit and open
            let commitments = zm.commit(&poly);

            let mut prover_transcriptcript = DefaultTranscript::new();
            let proof = zm.open(&point, &eval, &poly, Some(&mut prover_transcriptcript));

            let mut verifier_transcriptcript = DefaultTranscript::new();
            assert!(zm.verify(
                &point,
                &eval,
                &commitments,
                &proof,
                Some(&mut verifier_transcriptcript),
            ));

            //TODO: check both random oracles are synced
        }
    }

    #[test]
    fn prove_verify_batched() {
        let max_vars = 16;
        let num_polys = 8;
        let mut rng = &mut ChaCha20Rng::from_seed(*b"zeromorph_poly_commitment_scheme");
        let srs = ZeromorphSRS::setup(1 << (max_vars + 1), rng);

        for num_vars in 3..max_vars {
            // Setup
            let (pk, vk) = {
                let poly_size = 1 << (num_vars + 1);
                srs.trim(poly_size - 1).unwrap()
            };
            let zm = Zeromorph::<BLS12381AtePairing>::new(pk, vk);
            let polys: Vec<DenseMultilinearPolynomial<_>> = (0..num_polys)
                .map(|_| {
                    DenseMultilinearPolynomial::new(
                        (0..(1 << num_vars))
                            .map(|_| rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng))
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>();
            let point = (0..num_vars)
                .into_iter()
                .map(|_| rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng))
                .collect::<Vec<_>>();
            let evals = polys
                .clone()
                .into_iter()
                .map(|poly| poly.evaluate(point.clone()).unwrap())
                .collect::<Vec<_>>();

            // Commit and open
            let commitments: Vec<_> = polys.iter().map(|poly| zm.commit(poly)).collect();

            let mut prover_transcriptcript = DefaultTranscript::new();
            let proof = zm.open_batch(&point, &evals, &polys, Some(&mut prover_transcriptcript));

            let mut verifier_transcriptcript = DefaultTranscript::new();
            assert!(zm.verify_batch(
                &point,
                &evals,
                &commitments,
                &proof,
                Some(&mut verifier_transcriptcript),
            ))
            //TODO: check both random oracles are synced
        }
    }

    /// Test for computing qk given multilinear f
    /// Given 𝑓(𝑋₀, …, 𝑋ₙ₋₁), and `(𝑢, 𝑣)` such that \f(\u) = \v, compute `qₖ(𝑋₀, …, 𝑋ₖ₋₁)`
    /// such that the following identity holds:
    ///
    /// `𝑓(𝑋₀, …, 𝑋ₙ₋₁) − 𝑣 = ∑ₖ₌₀ⁿ⁻¹ (𝑋ₖ − 𝑢ₖ) qₖ(𝑋₀, …, 𝑋ₖ₋₁)`
    #[test]
    fn quotient_construction() {
        // Define size params
        let num_vars = 4;
        let n: u64 = 1 << num_vars;

        // Construct a random multilinear polynomial f, and (u,v) such that f(u) = v
        let mut rng = &mut ChaCha20Rng::from_seed(*b"zeromorph_poly_commitment_scheme");
        let multilinear_f = DenseMultilinearPolynomial::new(
            (0..n)
                .map(|_| rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng))
                .collect::<Vec<_>>(),
        );
        let u_challenge = (0..num_vars)
            .into_iter()
            .map(|_| rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng))
            .collect::<Vec<_>>();
        let v_evaluation = multilinear_f.evaluate(u_challenge.clone()).unwrap();

        // Compute multilinear quotients `qₖ(𝑋₀, …, 𝑋ₖ₋₁)`
        let (quotients, constant_term) =
            compute_multilinear_quotients::<BLS12381AtePairing>(&multilinear_f, &u_challenge);

        // Assert the constant term is equal to v_evaluation
        assert_eq!(
            constant_term, v_evaluation,
            "The constant term equal to the evaluation of the polynomial at challenge point."
        );

        //To demonstrate that q_k was properly constructd we show that the identity holds at a random multilinear challenge
        // i.e. 𝑓(𝑧) − 𝑣 − ∑ₖ₌₀ᵈ⁻¹ (𝑧ₖ − 𝑢ₖ)𝑞ₖ(𝑧) = 0
        let z_challenge = (0..num_vars)
            .map(|_| rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng))
            .collect::<Vec<_>>();

        let mut res = multilinear_f.evaluate(z_challenge.clone()).unwrap();
        res = res - v_evaluation;

        for (k, q_k_uni) in quotients.iter().enumerate() {
            let z_partial = z_challenge[&z_challenge.len() - k..].to_vec();
            //This is a weird consequence of how things are done.. the univariate polys are of the multilinear commitment in lagrange basis. Therefore we evaluate as multilinear
            let q_k = DenseMultilinearPolynomial::new(q_k_uni.coefficients.clone());
            let q_k_eval = q_k.evaluate(z_partial).unwrap();

            res = res
                - (&z_challenge[z_challenge.len() - k - 1]
                    - &u_challenge[z_challenge.len() - k - 1])
                    * q_k_eval;
        }
        assert_eq!(res, FieldElement::zero());
    }

    /// Test for construction of batched lifted degree quotient:
    ///  ̂q = ∑ₖ₌₀ⁿ⁻¹ yᵏ Xᵐ⁻ᵈᵏ⁻¹ ̂qₖ, 𝑑ₖ = deg(̂q), 𝑚 = 𝑁
    #[test]
    fn batched_lifted_degree_quotient() {
        const NUM_VARS: usize = 3;
        const N: usize = 1 << NUM_VARS;

        // Define mock qₖ with deg(qₖ) = 2ᵏ⁻¹
        let q_0 = Polynomial::new(&[FieldElement::one()]);
        let q_1 = Polynomial::new(&[FieldElement::from(2u64), FieldElement::from(3u64)]);
        let q_2 = Polynomial::new(&[
            FieldElement::from(4u64),
            FieldElement::from(5u64),
            FieldElement::from(6u64),
            FieldElement::from(7u64),
        ]);
        let quotients = vec![q_0, q_1, q_2];

        let mut rng = &mut ChaCha20Rng::from_seed(*b"zeromorph_poly_commitment_scheme");
        let y_challenge = rand_fr::<BLS12381AtePairing, ChaCha20Rng>(&mut rng);

        //Compute batched quptient  ̂q
        let batched_quotient = compute_batched_lifted_degree_quotient::<BLS12381AtePairing>(
            N,
            &quotients,
            &y_challenge,
        );

        //Explicitly define q_k_lifted = X^{N-2^k} * q_k and compute the expected batched result
        let q_0_lifted = Polynomial::new(&[
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::one(),
        ]);
        let q_1_lifted = Polynomial::new(&[
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::from(2u64),
            FieldElement::from(3u64),
        ]);
        let q_2_lifted = Polynomial::new(&[
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::from(4u64),
            FieldElement::from(5u64),
            FieldElement::from(6u64),
            FieldElement::from(7u64),
        ]);

        //Explicitly compute  ̂q i.e. RLC of lifted polys
        let mut batched_quotient_expected = Polynomial::zero();

        batched_quotient_expected = batched_quotient_expected + &q_0_lifted;
        batched_quotient_expected = batched_quotient_expected + &(q_1_lifted * y_challenge.clone());
        batched_quotient_expected =
            batched_quotient_expected + &(q_2_lifted * (&y_challenge * &y_challenge));
        assert_eq!(batched_quotient, batched_quotient_expected);
    }

    /// evaluated quotient \zeta_x
    ///
    /// 𝜁 = 𝑓 − ∑ₖ₌₀ⁿ⁻¹𝑦ᵏ𝑥ʷˢ⁻ʷ⁺¹𝑓ₖ  = 𝑓 − ∑_{d ∈ {d₀, ..., dₙ₋₁}} X^{d* - d + 1}  − ∑{k∶ dₖ=d} yᵏ fₖ , where d* = lifted degree
    ///
    /// 𝜁 =  ̂q - ∑ₖ₌₀ⁿ⁻¹ yᵏ Xᵐ⁻ᵈᵏ⁻¹ ̂qₖ, m = N
    #[test]
    fn partially_evaluated_quotient_zeta() {
        let num_vars = 3;
        let n: u64 = 1 << num_vars;

        let mut rng = &mut ChaCha20Rng::from_seed(*b"zeromorph_poly_commitment_scheme");
        let x_challenge = rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng);
        let y_challenge = rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng);

        let challenges: Vec<_> = (0..num_vars)
            .map(|_| rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng))
            .collect();
        let z_challenge = rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng);

        let (_, (zeta_x_scalars, _)) = eval_and_quotient_scalars::<BLS12381AtePairing>(
            &y_challenge,
            &x_challenge,
            &z_challenge,
            &challenges,
        );

        // To verify we manually compute zeta using the computed powers and expected
        // 𝜁 =  ̂q - ∑ₖ₌₀ⁿ⁻¹ yᵏ Xᵐ⁻ᵈᵏ⁻¹ ̂qₖ, m = N
        assert_eq!(zeta_x_scalars[0], -x_challenge.pow((n - 1) as u64));

        assert_eq!(
            zeta_x_scalars[1],
            -&y_challenge * x_challenge.pow((n - 1 - 1) as u64)
        );

        assert_eq!(
            zeta_x_scalars[2],
            -&y_challenge * y_challenge * x_challenge.pow((n - 3 - 1) as u64)
        );
    }

    /// Test efficiently computing 𝛷ₖ(x) = ∑ᵢ₌₀ᵏ⁻¹xⁱ
    /// 𝛷ₖ(𝑥) = ∑ᵢ₌₀ᵏ⁻¹𝑥ⁱ = (𝑥²^ᵏ − 1) / (𝑥 − 1)
    #[test]
    fn phi_n_x_evaluation() {
        const N: u64 = 8u64;
        let log_n = log_2(N as usize);

        // 𝛷ₖ(𝑥)
        let mut rng = &mut ChaCha20Rng::from_seed(*b"zeromorph_poly_commitment_scheme");
        let x_challenge = rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng);

        let efficient = (x_challenge.pow((1 << log_n) as u64)
            - FieldElement::<<BLS12381AtePairing as IsPairing>::BaseField>::one())
            / (&x_challenge - FieldElement::one());
        let expected: FieldElement<_> = phi::<BLS12381AtePairing>(&x_challenge, log_n);
        assert_eq!(efficient, expected);
    }

    /// Test efficiently computing 𝛷ₖ(x) = ∑ᵢ₌₀ᵏ⁻¹xⁱ
    /// 𝛷ₙ₋ₖ₋₁(𝑥²^ᵏ⁺¹) = (𝑥²^ⁿ − 1) / (𝑥²^ᵏ⁺¹ − 1)
    #[test]
    fn phi_n_k_1_x_evaluation() {
        const N: u64 = 8u64;
        let log_n = log_2(N as usize);

        // 𝛷ₖ(𝑥)
        let mut rng = &mut ChaCha20Rng::from_seed(*b"zeromorph_poly_commitment_scheme");
        let x_challenge = rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng);
        let k = 2;

        //𝑥²^ᵏ⁺¹
        let x_pow = x_challenge.pow((1 << (k + 1)) as u64);

        //(𝑥²^ⁿ − 1) / (𝑥²^ᵏ⁺¹ − 1)
        let efficient = (x_challenge.pow((1 << log_n) as u64)
            - FieldElement::<<BLS12381AtePairing as IsPairing>::BaseField>::one())
            / (x_pow - FieldElement::one());
        let expected: FieldElement<_> = phi::<BLS12381AtePairing>(&x_challenge, log_n - k - 1);
        assert_eq!(efficient, expected);
    }

    /// Test construction of 𝑍ₓ
    /// 𝑍ₓ =  ̂𝑓 − 𝑣 ∑ₖ₌₀ⁿ⁻¹(𝑥²^ᵏ𝛷ₙ₋ₖ₋₁(𝑥ᵏ⁺¹)− 𝑢ₖ𝛷ₙ₋ₖ(𝑥²^ᵏ)) ̂qₖ
    #[test]
    fn partially_evaluated_quotient_z_x() {
        let num_vars = 3;

        // Construct a random multilinear polynomial f, and (u,v) such that f(u) = v.
        let mut rng = &mut ChaCha20Rng::from_seed(*b"zeromorph_poly_commitment_scheme");
        let challenges: Vec<_> = (0..num_vars)
            .into_iter()
            .map(|_| rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng))
            .collect();

        let u_rev = {
            let mut res = challenges.clone();
            res.reverse();
            res
        };

        let x_challenge = rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng);
        let y_challenge = rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng);
        let z_challenge = rand_fr::<BLS12381AtePairing, &mut ChaCha20Rng>(&mut rng);

        // Construct Z_x scalars
        let (_, (_, z_x_scalars)) = eval_and_quotient_scalars::<BLS12381AtePairing>(
            &y_challenge,
            &x_challenge,
            &z_challenge,
            &challenges,
        );

        for k in 0..num_vars {
            let x_pow_2k = x_challenge.pow((1 << k) as u64); // x^{2^k}
            let x_pow_2kp1 = x_challenge.pow((1 << (k + 1)) as u64); // x^{2^{k+1}}
                                                                     // x^{2^k} * \Phi_{n-k-1}(x^{2^{k+1}}) - u_k *  \Phi_{n-k}(x^{2^k})
            let mut scalar = &x_pow_2k * &phi::<BLS12381AtePairing>(&x_pow_2kp1, num_vars - k - 1)
                - &u_rev[k] * &phi::<BLS12381AtePairing>(&x_pow_2k, num_vars - k);
            scalar *= &z_challenge;
            //TODO: this could be a trouble spot
            scalar = scalar.neg();
            assert_eq!(z_x_scalars[k], scalar);
        }
    }
}
