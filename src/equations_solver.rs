// -------------------------------------------------------------------------------------------------
//
//  This library implements various numerical algorithms which can be of use for all
//  kinds of scientific and engineering applications.
//
//  Copyright (c) 2025 by Dr. Panos Asproulis (p.asproulis@icloud.com).
//  All Rights Reserved.
//
// -------------------------------------------------------------------------------------------------

//!  This module implements various algorithms capable of solving equations.

//
// Import of external crates
//
use num_complex::Complex;
use num_traits::Float;

///
/// Solve a cubic equation of the form: a*h^3 + b*h^2 + c*h + d = 0
/// using Cardano's method. The cubic equation has three roots.
///
/// - Arguments:
///   - 'coefficients': the Bernoulli coefficients a, b, c, d.
///
/// - Returns:
///   - The roots of the cubic equation.
///
pub fn solve_cubic_equation<T: Float>(coefficients: &(T, T, T, T)) -> Vec<Complex<T>> {
    let (a, b, c, d) = coefficients;
    //
    // Normalize the coefficients
    //
    let b = *b / *a;
    let c = *c / *a;
    let d = *d / *a;
    //
    // Constants
    //
    let zero = T::zero();
    let one = T::one();
    let two = one + one;
    let half = one / two;
    let three = two + one;
    let one_third = one / three;
    let four = two + two;
    let nine = three * three;
    let twenty_seven = nine * three;
    //
    // Calculate the discriminant
    //
    let delta0 = b * b - three * c;
    let delta1 = two * b * b * b - nine * b * c + twenty_seven * d;
    //
    // Calculate ci as in Cardano's formula
    //
    let sqrt_term = Complex::new(delta1 * delta1 - four * delta0 * delta0 * delta0, zero).sqrt();
    let ci = ((Complex::new(delta1, zero) + sqrt_term) / two).powf(one_third);
    //
    // Cube roots of unity: ω and ω²
    //
    let omega = Complex::new(-half, three.sqrt() / two);
    //
    // Compute the roots
    //
    let mut roots: Vec<Complex<T>> = Vec::with_capacity(3);
    for k in 0..3 {
        let omega_k = omega.powi(k);
        let co = omega_k * ci;
        let bc = Complex::new(b, zero);
        let dc = Complex::new(delta0, zero);
        roots.push(-(bc + co + dc / co) / three);
    }

    roots
}
