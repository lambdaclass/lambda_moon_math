use crate::{
    cyclic_group::IsGroup,
    field::{element::FieldElement, traits::IsField},
    unsigned_integer::element::UnsignedInteger,
};

use super::naive::MSMError;
#[cfg(feature = "icicle")]
use crate::gpu::icicle::{icicle_msm, GpuMSMPoint};
use crate::traits::ByteConversion;
use alloc::vec;

/// This function computes the multiscalar multiplication (MSM).
///
/// Assume a group G of order r is given.
/// Let `points = [g_1, ..., g_n]` be a tuple of group points in G and
/// let `cs = [k_1, ..., k_n]` be a tuple of scalars in the Galois field GF(r).
///
/// Then, with additive notation, `msm(cs, points)` computes k_1 * g_1 + .... + k_n * g_n.
///
/// If `points` and `cs` are empty, then `msm` returns the zero element of the group.
///
/// Panics if `cs` and `points` have different lengths.
pub fn msm<const NUM_LIMBS: usize, F: IsField<BaseType = UnsignedInteger<NUM_LIMBS>>, G>(
    cs: &[FieldElement<F>],
    points: &[G],
) -> Result<G, MSMError>
where
    G: IsGroup + GpuMSMPoint,
    FieldElement<F>: ByteConversion,
{
    if cs.len() != points.len() {
        return Err(MSMError::LengthMismatch(cs.len(), points.len()));
    }

    #[cfg(feature = "icicle")]
    {
        if !G::curve_name().is_empty() {
            icicle_msm(cs, points)
        } else {
            println!(
                "Icicle msm failed for field {}. Program will fallback to CPU.",
                core::any::type_name::<F>()
            );
            let window_size = optimum_window_size(cs.len());
            let cs = cs.iter().map(|cs| *cs.value()).collect::<Vec<_>>();
            Ok(msm_with(&cs, points, window_size))
        }
    }

    #[cfg(not(feature = "icicle"))]
    {
        let window_size = optimum_window_size(cs.len());
        let cs = cs.iter().map(|cs| *cs.value()).collect::<Vec<_>>();
        Ok(msm_with(&cs, points, window_size))
    }
}

fn optimum_window_size(data_length: usize) -> usize {
    const SCALE_FACTORS: (usize, usize) = (4, 5);

    // We approximate the optimum window size with: f(n) = k * log2(n), where k is a scaling factor
    let len_isqrt = data_length.checked_ilog2().unwrap_or(0);
    (len_isqrt as usize * SCALE_FACTORS.0) / SCALE_FACTORS.1
}

pub fn msm_with<const NUM_LIMBS: usize, G>(
    cs: &[UnsignedInteger<NUM_LIMBS>],
    points: &[G],
    window_size: usize,
) -> G
where
    G: IsGroup,
{
    // When input is small enough, windows of length 2 seem faster than 1.
    const MIN_WINDOW_SIZE: usize = 2;
    const MAX_WINDOW_SIZE: usize = 32;

    let window_size = window_size.clamp(MIN_WINDOW_SIZE, MAX_WINDOW_SIZE);

    // The number of windows of size `s` is ceil(lambda/s).
    let num_windows = (64 * NUM_LIMBS - 1) / window_size + 1;

    // We define `buckets` outside of the loop so we only have to allocate once, and reuse it.
    //
    // This line forces a heap allocation which might be undesired. We can define this buckets
    // variable in the Pippenger struct to only allocate once, but use a bit of extra memory.
    // If we accept a const window_size, we could make it an array instaed of a vector
    // avoiding the heap allocation. We should be aware if that might be too agressive for
    // the stack and cause a potential stack overflow.
    let n_buckets = (1 << window_size) - 1;
    let mut buckets = vec![G::neutral_element(); n_buckets];

    (0..num_windows)
        .rev()
        .map(|window_idx| {
            // Put in the right bucket the corresponding ps[i] for the current window.
            cs.iter().zip(points).for_each(|(k, p)| {
                // We truncate the number to the least significative limb.
                // This is ok because window_size < usize::BITS.
                let window_unmasked = (k >> (window_idx * window_size)).limbs[NUM_LIMBS - 1];
                let m_ij = window_unmasked & n_buckets as u64;
                if m_ij != 0 {
                    let idx = (m_ij - 1) as usize;
                    buckets[idx] = buckets[idx].operate_with(p);
                }
            });

            // Do the reduction step for the buckets.
            buckets
                .iter_mut()
                // This first part iterates buckets in descending order, generating an iterator with the sum of
                // each bucket and all that came before as its items; i.e: (b_n, b_n + b_n-1, ..., b_n + ... + b_0)
                .rev()
                .scan(G::neutral_element(), |m, b| {
                    *m = m.operate_with(b); // Reduction step.
                    *b = G::neutral_element(); // Cleanup bucket slot to reuse in the next window.
                    Some(m.clone())
                })
                // This next part sums all elements of the iterator: (b_n) + (b_n + b_n-1) + ...
                // This results in: (n + 1) * b_n + n * b_n-1 + ... + b_0
                .reduce(|g, m| g.operate_with(&m))
                .unwrap_or_else(G::neutral_element)
        })
        // NOTE: this operation is non-associative and strictly sequential
        .reduce(|t, g| t.operate_with_self(1_u64 << window_size).operate_with(&g))
        .unwrap_or_else(G::neutral_element)
}

#[cfg(feature = "parallel")]
// It has the following differences with the sequential one:
//  1. It uses one vec per thread to store buckets.
//  2. It reduces all window results via a different method.
pub fn parallel_msm_with<const NUM_LIMBS: usize, G>(
    cs: &[UnsignedInteger<NUM_LIMBS>],
    points: &[G],
    window_size: usize,
) -> G
where
    G: IsGroup + Send + Sync,
{
    use rayon::prelude::*;

    assert!(window_size < usize::BITS as usize); // Program would go OOM anyways

    // The number of windows of size `s` is ceil(lambda/s).
    let num_windows = (64 * NUM_LIMBS - 1) / window_size + 1;
    let n_buckets = (1 << window_size) - 1;

    // TODO: limit the number of threads, and reuse vecs
    (0..num_windows)
        .into_par_iter()
        .map(|window_idx| {
            let mut buckets = vec![G::neutral_element(); n_buckets];
            // Put in the right bucket the corresponding ps[i] for the current window.
            let shift = window_idx * window_size;
            cs.iter().zip(points).for_each(|(k, p)| {
                // We truncate the number to the least significative limb.
                // This is ok because window_size < usize::BITS.
                let window_unmasked = (k >> shift).limbs[NUM_LIMBS - 1];
                let m_ij = window_unmasked & n_buckets as u64;
                if m_ij != 0 {
                    let idx = (m_ij - 1) as usize;
                    buckets[idx] = buckets[idx].operate_with(p);
                }
            });

            let mut m = G::neutral_element();

            // Do the reduction step for the buckets.
            let window_item = buckets
                // NOTE: changing this into a parallel iter drops performance, because of the
                //  need to use multiplication in the `map` step
                .into_iter()
                .rev()
                .map(|b| {
                    m = m.operate_with(&b); // Reduction step.
                    m.clone()
                })
                .reduce(|g, m| g.operate_with(&m))
                .unwrap_or_else(G::neutral_element);

            window_item.operate_with_self(UnsignedInteger::<NUM_LIMBS>::from_u64(1) << shift)
        })
        .reduce(G::neutral_element, |a, b| a.operate_with(&b))
}

#[cfg(test)]
mod tests {
    use crate::cyclic_group::IsGroup;
    use crate::msm::{naive, pippenger};
    use crate::{
        elliptic_curve::{
            short_weierstrass::curves::bls12_381::curve::BLS12381Curve, traits::IsEllipticCurve,
        },
        unsigned_integer::element::UnsignedInteger,
    };
    use alloc::{format, vec::Vec};
    use proptest::{collection, prelude::*, prop_assert_eq, prop_compose, proptest};

    const _CASES: u32 = 20;
    const _MAX_WSIZE: usize = 8;
    const _MAX_LEN: usize = 30;

    prop_compose! {
        fn unsigned_integer()(limbs: [u64; 6]) -> UnsignedInteger<6> {
            UnsignedInteger::from_limbs(limbs)
        }
    }

    prop_compose! {
        fn unsigned_integer_vec()(vec in collection::vec(unsigned_integer(), 0.._MAX_LEN)) -> Vec<UnsignedInteger<6>> {
            vec
        }
    }

    prop_compose! {
        fn point()(power: u128) -> <BLS12381Curve as IsEllipticCurve>::PointRepresentation {
            BLS12381Curve::generator().operate_with_self(power)
        }
    }

    prop_compose! {
        fn points_vec()(vec in collection::vec(point(), 0.._MAX_LEN)) -> Vec<<BLS12381Curve as IsEllipticCurve>::PointRepresentation> {
            vec
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig {
            cases: _CASES, .. ProptestConfig::default()
          })]
        // Property-based test that ensures `pippenger::msm` gives same result as `naive::msm`.
        #[test]
        fn test_pippenger_matches_naive_msm(window_size in 1.._MAX_WSIZE, cs in unsigned_integer_vec(), points in points_vec()) {
            let min_len = cs.len().min(points.len());
            let cs = cs[..min_len].to_vec();
            let points = points[..min_len].to_vec();

            let pippenger = pippenger::msm_with(&cs, &points, window_size);
            let naive = naive::msm(&cs, &points).unwrap();

            prop_assert_eq!(naive, pippenger);
        }

        // Property-based test that ensures `pippenger::msm_with` gives same result as `pippenger::parallel_msm_with`.
        #[test]
        #[cfg(feature = "parallel")]
        fn test_parallel_pippenger_matches_sequential(window_size in 1.._MAX_WSIZE, cs in unsigned_integer_vec(), points in points_vec()) {
            let min_len = cs.len().min(points.len());
            let cs = cs[..min_len].to_vec();
            let points = points[..min_len].to_vec();

            let sequential = pippenger::msm_with(&cs, &points, window_size);
            let parallel = pippenger::parallel_msm_with(&cs, &points, window_size);

            prop_assert_eq!(parallel, sequential);
        }
    }
}
