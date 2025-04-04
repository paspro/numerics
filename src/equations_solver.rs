// -------------------------------------------------------------------------------------------------
//
//  This library implements various numerical algorithms which can be of use for all
//  kinds of scientific and engineering applications.
//
//  Copyright (c) 2025 by Dr. Panos Asproulis (p.asproulis@icloud.com).
//  All Rights Reserved.
//
// -------------------------------------------------------------------------------------------------

//! This module implements various algorithms capable of solving equations.

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
pub fn solve_cubic_equation<T: Float + std::fmt::Debug>(
    coefficients: &(T, T, T, T),
) -> Vec<Complex<T>> {
    let (a, b, c, d) = coefficients;
    //
    // Normalize the coefficients.
    //
    let b = *b / *a;
    let c = *c / *a;
    let d = *d / *a;
    //
    // Constants.
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
    // Calculate the discriminant.
    //
    let delta0 = b * b - three * c;
    let delta1 = two * b * b * b - nine * b * c + twenty_seven * d;
    //
    // Apply Cardano's formula.
    //
    let sqrt_term: Complex<T>;
    if delta0.abs() <= T::epsilon() {
        if delta1.abs() <= T::epsilon() {
            //
            // All roots are equal.
            //
            return vec![Complex::new(-b / three, zero)];
        } else {
            //
            // We need to respect the sigh of delta1.
            //
            sqrt_term = Complex::new(delta1, zero);
        }
    } else {
        //
        // General case.
        //
        sqrt_term = Complex::new(delta1 * delta1 - four * delta0 * delta0 * delta0, zero).sqrt();
    }
    let ci = ((Complex::new(delta1, zero) + sqrt_term) / two).powf(one_third);
    println!(
        "ci: {:?}, delta1: {:?}, sqrt_term: {:?}",
        ci, delta1, sqrt_term
    );
    //
    // Cube roots of unity: ω and ω².
    //
    let omega = Complex::new(-half, three.sqrt() / two);
    //
    // Compute the roots.
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_cubic_real_roots_spaced() {
        //
        // x³ - 6x² + 11x - 6 = 0 has roots 1, 2, 3.
        //
        let coefficients = (1.0, -6.0, 11.0, -6.0);
        let roots = solve_cubic_equation(&coefficients);
        //
        // Sort roots by real part.
        //
        let mut sorted_roots = roots.clone();
        sorted_roots.sort_by(|a, b| a.re.partial_cmp(&b.re).unwrap());
        assert_relative_eq!(sorted_roots[0].re, 1.0, epsilon = 1e-10);
        assert_relative_eq!(sorted_roots[1].re, 2.0, epsilon = 1e-10);
        assert_relative_eq!(sorted_roots[2].re, 3.0, epsilon = 1e-10);
        //
        // Check imaginary parts.
        //
        for root in &sorted_roots {
            assert_relative_eq!(root.im, 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_cubic_complex_conjugate_roots() {
        //
        // x³ - x² + x - 1 = 0 has one real root and two complex conjugate roots.
        //
        let coefficients = (1.0, -1.0, 1.0, -1.0);
        let roots = solve_cubic_equation(&coefficients);
        //
        // Find the real root and the complex roots.
        //
        let mut real_roots = Vec::new();
        let mut complex_roots = Vec::new();

        for root in &roots {
            if root.im.abs() < 1e-10 {
                real_roots.push(root);
            } else {
                complex_roots.push(root);
            }
        }

        assert_eq!(real_roots.len(), 1);
        assert_eq!(complex_roots.len(), 2);
        //
        // Verify real root is approximately 1.
        //
        assert_relative_eq!(real_roots[0].re, 1.0, epsilon = 1e-10);
        //
        // Verify complex roots are conjugates of each other.
        //
        assert_relative_eq!(complex_roots[0].re, complex_roots[1].re, epsilon = 1e-10);
        assert_relative_eq!(complex_roots[0].im, -complex_roots[1].im, epsilon = 1e-10);
    }

    #[test]
    fn test_cubic_with_zero_coefficient() {
        //
        // x³ + 0x² + 0x - 8 = 0 has roots 2, -1+√3i, -1-√3i.
        //
        let coefficients = (1.0, 0.0, 0.0, -8.0);
        let roots = solve_cubic_equation(&coefficients);
        //
        // Find the real root.
        //
        let real_roots: Vec<&Complex<f64>> = roots.iter().filter(|r| r.im.abs() < 1e-10).collect();
        let complex_roots: Vec<&Complex<f64>> =
            roots.iter().filter(|r| r.im.abs() >= 1e-10).collect();

        assert_eq!(real_roots.len(), 1);
        assert_eq!(complex_roots.len(), 2);
        //
        // Real root should be 2.
        //
        assert_relative_eq!(real_roots[0].re, 2.0, epsilon = 1e-10);
        //
        // Complex roots should be -1 ± √3i.
        //
        assert_relative_eq!(complex_roots[0].re, -1.0, epsilon = 1e-10);
        assert_relative_eq!(complex_roots[1].re, -1.0, epsilon = 1e-10);
        //
        // |im| should be √3.
        //
        assert_relative_eq!(
            complex_roots[0].im.abs(),
            1.732050807568877,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_cubic_verify_by_substitution() {
        //
        // Generic test with arbitrary coefficients.
        //
        let coefficients = (2.5, -4.3, 7.1, -9.8);
        let roots = solve_cubic_equation(&coefficients);

        for root in &roots {
            let (a, b, c, d) = coefficients;
            //
            // Compute P(root) = a*root³ + b*root² + c*root + d.
            //
            let result = Complex::new(a, 0.0) * root.powi(3)
                + Complex::new(b, 0.0) * root.powi(2)
                + Complex::new(c, 0.0) * root
                + Complex::new(d, 0.0);
            //
            // Check if P(root) ≈ 0.
            //
            assert_relative_eq!(result.norm(), 0.0, epsilon = 1e-8);
        }
    }
}
