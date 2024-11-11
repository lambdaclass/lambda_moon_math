extern crate alloc;
use crate::circle::point::CirclePoint;
use crate::field::fields::mersenne31::field::Mersenne31Field;
use alloc::vec::Vec;

/// Given g_n, a generator of the subgroup of size n of the circle, i.e. <g_n>,
/// and given a shift, that is a another point of the circle,
/// we define the coset shift + <g_n> which is the set of all the points in
/// <g_n> plus the shift.
/// For example, if <g_4> = {p1, p2, p3, p4}, then g_8 + <g_4> = {g_8 + p1, g_8 + p2, g_8 + p3, g_8 + p4}.

#[derive(Debug, Clone)]
pub struct Coset {
    // Coset: shift + <g_n> where n = 2^{log_2_size}.
    // Example: g_16 + <g_8>, n = 8, log_2_size = 3, shift = g_16.
    pub log_2_size: u32, //TODO: Change log_2_size to u8 because log_2_size < 31.a
    pub shift: CirclePoint<Mersenne31Field>,
}

impl Coset {
    pub fn new(log_2_size: u32, shift: CirclePoint<Mersenne31Field>) -> Self {
        Coset { log_2_size, shift }
    }

    /// Returns the coset g_2n + <g_n>
    pub fn new_standard(log_2_size: u32) -> Self {
        // shift is a generator of the subgroup of order 2n = 2^{log_2_size + 1}.
        let shift = CirclePoint::get_generator_of_subgroup(log_2_size + 1);
        Coset { log_2_size, shift }
    }

    /// Returns g_n, the generator of the subgroup of order n = 2^log_2_size.
    pub fn get_generator(&self) -> CirclePoint<Mersenne31Field> {
        CirclePoint::GENERATOR.repeated_double(31 - self.log_2_size)
    }

    /// Given a standard coset g_2n + <g_n>, returns the subcoset with half size g_2n + <g_{n/2}>
    pub fn half_coset(coset: Self) -> Self {
        Coset {
            log_2_size: coset.log_2_size - 1,
            shift: coset.shift,
        }
    }  

    /// Given a coset shift + G returns the coset -shift + G.
    /// Note that (g_2n + <g_{n/2}>) U (-g_2n + <g_{n/2}>) = g_2n + <g_n>.
    pub fn conjugate(coset: Self) -> Self {
        Coset {
            log_2_size: coset.log_2_size,
            shift: coset.shift.conjugate(),
        }
    }

    /// Returns the vector of shift + g for every g in <g_n>.
    /// where g = i * g_n for i = 0, ..., n-1.
    #[cfg(feature = "alloc")]
    pub fn get_coset_points(coset: &Self) -> Vec<CirclePoint<Mersenne31Field>> {
        // g_n the generator of the subgroup of order n.
        let generator_n = CirclePoint::get_generator_of_subgroup(coset.log_2_size);
        let size: u8 = 1 << coset.log_2_size;
        core::iter::successors(Some(coset.shift.clone()), move |prev| {
            Some(prev + &generator_n)
        })
        .take(size.into())
        .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coset_points_vector_has_right_size() {
        let coset = Coset::new_standard(3);
        let points = Coset::get_coset_points(&coset);
        assert_eq!(1 << coset.log_2_size, points.len())
    }

    #[test]
    fn antipode_of_coset_point_is_in_coset() {
        let coset = Coset::new_standard(3);
        let points = Coset::get_coset_points(&coset);
        let point = points[2].clone();
        let anitpode_point = points[6].clone();
        assert_eq!(anitpode_point, point.antipode())
    }

    #[test]
    fn coset_generator_has_right_order() {
        let coset = Coset::new(2, CirclePoint::GENERATOR * 3);
        let generator_n = coset.get_generator();
        assert_eq!(generator_n.repeated_double(2), CirclePoint::zero());
    }
}
