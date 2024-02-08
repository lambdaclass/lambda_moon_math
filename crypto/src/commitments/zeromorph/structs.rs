use core::mem;

use lambdaworks_math::{
    cyclic_group::IsGroup,
    elliptic_curve::traits::IsPairing,
    errors::DeserializationError,
    field::{element::FieldElement, traits::IsPrimeField},
    traits::{AsBytes, ByteConversion, Deserializable},
};

//TODO: have own random gen type
use rand::RngCore;

//TODO: gate with alloc
#[derive(Debug, Clone, Default)]
pub struct ZeromorphSRS<P: IsPairing> {
    pub g1_powers: Vec<P::G1Point>,
    pub g2_powers: Vec<P::G2Point>,
}

impl<P: IsPairing> ZeromorphSRS<P> {
    /// Performs
    // TODO(pat_stiles): accept a Rng or Rng wrapper.
    // TODO(pat_stiles): Quality of life improvements in regards to random sampling
    // TODO: errors lengths are valid
    // TODO: Result
    pub fn setup<R: RngCore>(max_degree: usize, rng: &mut R) -> ZeromorphSRS<P>
    where
        <P as IsPairing>::BaseField: IsPrimeField,
        FieldElement<<P as IsPairing>::BaseField>: ByteConversion,
    {
        let mut bytes = [0u8; 384];
        rng.fill_bytes(&mut bytes);
        let tau = FieldElement::<P::BaseField>::from_bytes_be(&bytes).unwrap();
        rng.fill_bytes(&mut bytes);
        let g1_scalar = FieldElement::<P::BaseField>::from_bytes_be(&bytes).unwrap();
        let g1 = P::G1Point::neutral_element().operate_with_self(g1_scalar.representative());
        rng.fill_bytes(&mut bytes);
        let g2_scalar = FieldElement::<P::BaseField>::from_bytes_be(&bytes).unwrap();
        let g2 = P::G2Point::neutral_element().operate_with_self(g2_scalar.representative());

        //TODO: Add window table for msm this is slow af
        let tau_powers = (0..=max_degree)
            .scan(tau.clone(), |state, _| {
                let val = state.clone();
                *state *= &tau;
                Some(val)
            })
            .map(|scalar| scalar.representative())
            .collect::<Vec<_>>();

        let g1s: Vec<_> = (0..=max_degree).map(|e| g1.operate_with_self(e)).collect();
        let g2s: Vec<_> = (0..=max_degree).map(|e| g2.operate_with_self(e)).collect();

        //TODO: gate with rayon
        let g1_powers = g1s
            .iter()
            .zip(tau_powers.iter())
            .map(|(g1, tau)| g1.operate_with_self(*tau))
            .collect();
        let g2_powers = g2s
            .iter()
            .zip(tau_powers.iter())
            .map(|(g2, tau)| g2.operate_with_self(*tau))
            .collect();

        ZeromorphSRS {
            g1_powers,
            g2_powers,
        }
    }

    pub fn new(g1_powers: &[P::G1Point], g2_powers: &[P::G2Point]) -> Self {
        Self {
            g1_powers: g1_powers.to_vec(),
            g2_powers: g2_powers.to_vec(),
        }
    }

    // TODO: errors lengths are valid
    pub fn trim(
        &self,
        max_degree: usize,
    ) -> Result<(ZeromorphProverKey<P>, ZeromorphVerifierKey<P>), ZeromorphError> {
        if max_degree > self.g1_powers.len() {
            return Err(ZeromorphError::KeyLengthError(
                max_degree,
                self.g1_powers.len(),
            ));
        }
        let offset = self.g1_powers.len() - max_degree;
        let offset_g1_powers = self.g1_powers[offset..].to_vec();
        Ok((
            ZeromorphProverKey {
                g1_powers: self.g1_powers.clone(),
                offset_g1_powers: offset_g1_powers,
            },
            ZeromorphVerifierKey {
                g1: self.g1_powers[0].clone(),
                g2: self.g2_powers[0].clone(),
                tau_g2: self.g2_powers[1].clone(),
                tau_n_max_sub_2_n: self.g2_powers[offset].clone(),
            },
        ))
    }
}

#[cfg(feature = "std")]
impl<P: IsPairing> ZeromorphSRS<P>
where
    <P as IsPairing>::G1Point: Deserializable,
    <P as IsPairing>::G2Point: Deserializable,
{
    pub fn from_file(file_path: &str) -> Result<Self, crate::errors::SrsFromFileError> {
        let bytes = std::fs::read(file_path)?;
        Ok(Self::deserialize(&bytes)?)
    }
}

impl<P: IsPairing> AsBytes for ZeromorphSRS<P>
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
        let mut g1_powers_len_bytes: Vec<u8> = self.g1_powers.len().to_le_bytes().to_vec();

        // For architectures with less than 64 bits for pointers
        // We add extra zeros at the end
        while g1_powers_len_bytes.len() < 8 {
            g1_powers_len_bytes.push(0)
        }

        serialized_data.extend(&g1_powers_len_bytes);

        // third 8 bytes store the amount of G2 elements to be stored, this is more than can be indexed with a 64-bit architecture, and some millions of terabytes of data if the points were compressed
        let mut g2_powers_len_bytes: Vec<u8> = self.g2_powers.len().to_le_bytes().to_vec();

        // For architectures with less than 64 bits for pointers
        // We add extra zeros at the end
        while g2_powers_len_bytes.len() < 8 {
            g2_powers_len_bytes.push(0)
        }

        serialized_data.extend(&g2_powers_len_bytes);

        // G1 elements
        for point in &self.g1_powers {
            serialized_data.extend(point.as_bytes());
        }

        // G2 elements
        for point in &self.g2_powers {
            serialized_data.extend(point.as_bytes());
        }

        serialized_data
    }
}

impl<P: IsPairing> Deserializable for ZeromorphSRS<P>
where
    <P as IsPairing>::G1Point: Deserializable,
    <P as IsPairing>::G2Point: Deserializable,
{
    fn deserialize(bytes: &[u8]) -> Result<Self, DeserializationError> {
        const VERSION_OFFSET: usize = 4;
        const G1_LEN_OFFSET: usize = 12;
        const G2_LEN_OFFSET: usize = 20;

        let g1_powers_len_u64 = u64::from_le_bytes(
            // This unwrap can't fail since we are fixing the size of the slice
            bytes[VERSION_OFFSET..G1_LEN_OFFSET].try_into().unwrap(),
        );

        let g1_powers_len = usize::try_from(g1_powers_len_u64)
            .map_err(|_| DeserializationError::PointerSizeError)?;

        let g2_powers_len_u64 = u64::from_le_bytes(
            // This unwrap can't fail since we are fixing the size of the slice
            bytes[G1_LEN_OFFSET..G2_LEN_OFFSET].try_into().unwrap(),
        );

        let g2_powers_len = usize::try_from(g2_powers_len_u64)
            .map_err(|_| DeserializationError::PointerSizeError)?;

        let mut g1_powers: Vec<P::G1Point> = Vec::new();
        let mut g2_powers: Vec<P::G2Point> = Vec::new();

        let size_g1_point = mem::size_of::<P::G1Point>();
        let size_g2_point = mem::size_of::<P::G2Point>();

        for i in 0..g1_powers_len {
            // The second unwrap shouldn't fail since the amount of bytes is fixed
            let point = P::G1Point::deserialize(
                bytes[i * size_g1_point + G2_LEN_OFFSET
                    ..i * size_g1_point + size_g1_point + G2_LEN_OFFSET]
                    .try_into()
                    .unwrap(),
            )?;
            g1_powers.push(point);
        }

        let g2s_offset = size_g1_point * g1_powers_len + G2_LEN_OFFSET;
        for i in 0..g2_powers_len {
            // The second unwrap shouldn't fail since the amount of bytes is fixed
            let point = P::G2Point::deserialize(
                bytes[i * size_g2_point + g2s_offset..i * size_g2_point + g2s_offset]
                    .try_into()
                    .unwrap(),
            )?;
            g2_powers.push(point);
        }

        let srs = ZeromorphSRS::new(&g1_powers, &g2_powers);
        Ok(srs)
    }
}

#[derive(Clone, Debug)]
pub struct ZeromorphProof<P: IsPairing> {
    pub pi: P::G1Point,
    pub q_hat_com: P::G1Point,
    pub q_k_com: Vec<P::G1Point>,
}

#[derive(Clone, Debug)]
pub struct ZeromorphProverKey<P: IsPairing> {
    pub g1_powers: Vec<P::G1Point>,
    pub offset_g1_powers: Vec<P::G1Point>,
}

#[derive(Copy, Clone, Debug)]
pub struct ZeromorphVerifierKey<P: IsPairing> {
    pub g1: P::G1Point,
    pub g2: P::G2Point,
    pub tau_g2: P::G2Point,
    pub tau_n_max_sub_2_n: P::G2Point,
}

#[derive(Debug)]
pub enum ZeromorphError {
    //#[error("Length Error: SRS Length: {0}, Key Length: {0}")]
    KeyLengthError(usize, usize),
}
