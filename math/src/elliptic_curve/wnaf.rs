use crate::{
    cyclic_group::IsGroup,
    elliptic_curve::short_weierstrass::{
        point::ShortWeierstrassProjectivePoint, traits::IsShortWeierstrass,
    },
    field::{element::FieldElement, traits::IsPrimeField},
    traits::ByteConversion,
};
use alloc::{vec, vec::Vec};
use core::marker::PhantomData;

#[cfg(feature = "parallel")]
use rayon::iter::{
    IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator,
};

extern crate std; // To be able to use f64::ln()

pub struct WnafTable<EC, ScalarField>
where
    EC: IsShortWeierstrass<PointRepresentation = ShortWeierstrassProjectivePoint<EC>>,
    EC::PointRepresentation: Send + Sync,
    ScalarField: IsPrimeField + Sync,
    FieldElement<ScalarField>: ByteConversion + Send + Sync,
{
    table: Vec<Vec<ShortWeierstrassProjectivePoint<EC>>>,
    window_size: usize,
    phantom: PhantomData<ScalarField>,
}

impl<EC, ScalarField> WnafTable<EC, ScalarField>
where
    EC: IsShortWeierstrass<PointRepresentation = ShortWeierstrassProjectivePoint<EC>>,
    EC::PointRepresentation: Send + Sync,
    ScalarField: IsPrimeField + Sync,
    FieldElement<ScalarField>: ByteConversion + Send + Sync,
{
    pub fn new(base: &ShortWeierstrassProjectivePoint<EC>, max_num_of_scalars: usize) -> Self {
        let scalar_field_bit_size = ScalarField::field_bit_size();
        let window = Self::get_mul_window_size(max_num_of_scalars);
        let in_window = 1 << window;
        let outerc = (scalar_field_bit_size + window - 1) / window;
        let last_in_window = 1 << (scalar_field_bit_size - (outerc - 1) * window);

        let mut g_outer = base.clone();
        let mut g_outers = Vec::with_capacity(outerc);
        for _ in 0..outerc {
            g_outers.push(g_outer.clone());
            for _ in 0..window {
                g_outer = g_outer.double();
            }
        }

        let mut table =
            vec![vec![ShortWeierstrassProjectivePoint::<EC>::neutral_element(); in_window]; outerc];

        #[cfg(feature = "parallel")]
        let iter = table.par_iter_mut();
        #[cfg(not(feature = "parallel"))]
        let iter = table.iter_mut();

        iter.enumerate().take(outerc).zip(g_outers).for_each(
            |((outer, multiples_of_g), g_outer)| {
                let curr_in_window = if outer == outerc - 1 {
                    last_in_window
                } else {
                    in_window
                };

                let mut g_inner = ShortWeierstrassProjectivePoint::<EC>::neutral_element();
                for inner in multiples_of_g.iter_mut().take(curr_in_window) {
                    *inner = g_inner.clone();
                    g_inner = g_inner.operate_with(&g_outer);
                }
            },
        );

        Self {
            table,
            window_size: window,
            phantom: PhantomData,
        }
    }

    pub fn multi_scalar_mul(
        &self,
        v: &[FieldElement<ScalarField>],
    ) -> Vec<ShortWeierstrassProjectivePoint<EC>> {
        #[cfg(feature = "parallel")]
        let iter = v.par_iter();
        #[cfg(not(feature = "parallel"))]
        let iter = v.iter();

        iter.map(|e| self.windowed_mul(e.clone())).collect()
    }

    fn windowed_mul(
        &self,
        scalar: FieldElement<ScalarField>,
    ) -> ShortWeierstrassProjectivePoint<EC> {
        let mut res = self.table[0][0].clone();

        let modulus_size = ScalarField::field_bit_size();
        let outerc = (modulus_size + self.window_size - 1) / self.window_size;
        let scalar_bits_le: Vec<bool> = scalar
            .to_bytes_le()
            .iter()
            .flat_map(|byte| (0..8).map(|i| (byte >> i) & 1 == 1).collect::<Vec<_>>())
            .collect();

        for outer in 0..outerc {
            let mut inner = 0usize;
            for i in 0..self.window_size {
                if outer * self.window_size + i < modulus_size
                    && scalar_bits_le[outer * self.window_size + i]
                {
                    inner |= 1 << i;
                }
            }
            res = res.operate_with(&self.table[outer][inner]);
        }

        res.to_affine()
    }

    fn get_mul_window_size(max_num_of_scalars: usize) -> usize {
        if max_num_of_scalars < 32 {
            3
        } else {
            f64::ln(max_num_of_scalars as f64).ceil() as usize
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        elliptic_curve::{
            short_weierstrass::curves::bls12_381::{
                curve::BLS12381Curve,
                default_types::{FrElement, FrField},
            },
            traits::IsEllipticCurve,
        },
        unsigned_integer::element::U256,
    };
    use rand::*;
    use std::time::Instant;

    #[test]
    fn wnaf_works() {
        let point_count = 100;
        let g1 = BLS12381Curve::generator();

        let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(9001);
        let mut scalars = Vec::new();
        for _i in 0..point_count {
            scalars.push(FrElement::new(U256::from(rng.gen::<u128>())));
        }

        let start1 = Instant::now();
        let naive_result: Vec<_> = scalars
            .iter()
            .map(|scalar| {
                g1.operate_with_self(scalar.clone().representative())
                    .to_affine()
            })
            .collect();
        let duration1 = start1.elapsed();
        println!(
            "Time taken for naive ksk with {} scalars: {:?}",
            point_count, duration1
        );

        let start2 = Instant::now();
        let wnaf_result =
            WnafTable::<BLS12381Curve, FrField>::new(&g1, point_count).multi_scalar_mul(&scalars);
        let duration2 = start2.elapsed();
        println!(
            "Time taken for wnaf msm including table generation with {} scalars: {:?}",
            point_count, duration2
        );

        for i in 0..point_count {
            assert_eq!(naive_result[i], wnaf_result[i]);
        }
    }
}
