// -------------------------------------------------------------------------------------------------
//
//  This library implements various numerical algorithms which can be of use for all
//  kinds of scientific and engineering applications.
//
//  Copyright (c) 2025 by Dr. Panos Asproulis (p.asproulis@icloud.com).
//  All Rights Reserved.
//
// -------------------------------------------------------------------------------------------------

//!  This module implements various root finding algorithms.

use num_traits::Float;

///
/// This function reverses a generic equation f(x) by computing the value
/// of "x0" when the value of the equation f(x0) is known with the use of
/// an accelerated Secant root-finder algorithm operating on the function
/// g(x) = f(x) - f0, where f0 = f(x0) -> g(x0) = 0.
///
/// - Arguments:
///   - `f`: The function to invert.
///   - `f0`: The value of f0 = f(x0).
///   - `guess1`: A guess for the value of x0.
///   - `guess2`: A second guess for the value of x0.
///   - `tolerance`: The tolerance to utilise for the numerical method.
///   - `max_iterations`: The maximum number of allowed iterations.
///   - `n_order`: The interpolation order used by the method.
///
/// - Returns:
///   - The value of "x0" when f(x0) is known.
///
#[allow(clippy::too_many_arguments)]
pub fn function_inverter<T: Float>(
    f: &dyn Fn(T) -> T,
    f0: T,
    guess1: T,
    guess2: T,
    tolerance: T,
    max_iterations: i32,
    n_order: i32,
) -> T {
    if n_order < 0 {
        eprintln!("\nError in Function: function_inverter");
        eprintln!("The order of interpolation must be a positive integer.");
        eprintln!();
        return T::nan();
    }
    ///
    /// The solution context for the function inverter.
    ///
    struct SolutionContext<'a, T: Float> {
        ///
        /// The function to invert.
        ///
        f: &'a dyn Fn(T) -> T,
        ///
        /// The value of f0 = f(x0).
        ///
        f0: T,
        ///
        /// A guess for the value of x0.
        ///
        guess1: T,
        ///
        /// A second guess for the value of x0.
        ///
        guess2: T,
        ///
        /// The interpolation order used by the method.
        ///
        n_order: i32,
    }
    //
    // Implement the solution context.
    //
    impl<'a, T: Float> SolutionContext<'a, T> {
        ///
        /// Create a new solution context.
        ///
        /// - Arguments:
        ///   - `f`: The function to invert.
        ///   - `f0`: The value of f0 = f(x0).
        ///   - `guess1`: A guess for the value of x0.
        ///   - `guess2`: A second guess for the value of x0.
        ///   - `n_order`: The interpolation order used by the method.
        ///
        /// - Returns:
        ///   - The solution context.
        ///
        fn new(f: &'a dyn Fn(T) -> T, f0: T, guess1: T, guess2: T, n_order: i32) -> Self {
            Self {
                f,
                f0,
                guess1,
                guess2,
                n_order,
            }
        }
        ///
        /// Compute the maximum value for the interpolation.
        ///
        /// - Arguments:
        ///   - `p`: The value of p.
        ///
        /// - Returns:
        ///   - The maximum value for the interpolation.
        ///
        fn i_max(&self, p: i32) -> i32 {
            if p == -1 || p == 0 {
                0
            } else if p <= self.n_order {
                p - 1
            } else {
                self.n_order
            }
        }
        ///
        /// Compute the solution for the function inverter.
        ///
        /// - Arguments:
        ///  - `p`: The value of p.
        ///  - `i`: The value of i.
        ///
        /// - Returns:
        ///   - The solution for the function inverter.
        ///
        fn solution(&mut self, p: i32, i: i32) -> T {
            if p == -1 {
                self.guess1
            } else if p == 0 {
                self.guess2
            } else if i == 0 {
                //
                // Secant step.
                //
                let x1 = self.solution(p - 1, self.i_max(p - 1));
                let x2 = self.solution(p - 2, self.i_max(p - 2));
                let f1 = (self.f)(x1) - self.f0;
                let f2 = (self.f)(x2) - self.f0;
                x1 - f1 * (x1 - x2) / (f1 - f2)
            } else {
                //
                // i-th order approximation.
                //
                let x1 = self.solution(p - 1, i - 1);
                let x2 = self.solution(p - 1, self.i_max(p - 1));
                let x3 = self.solution(p, i - 1);
                let x4 = self.solution(p - i - 2, self.i_max(p - i - 2));
                x3 + (x2 - x3) * (x1 - x3) / (x1 + x2 - x3 - x4)
            }
        }
    }

    let mut context = SolutionContext::new(f, f0, guess1, guess2, n_order);
    let mut s = guess2;

    for iter in (n_order + 1)..=max_iterations {
        //
        // Compute the new value.
        //
        s = context.solution(iter, n_order);
        let error = (f(s) - f0).abs();
        //
        // Check for convergence.
        //
        if error <= tolerance {
            break;
        }
    }

    s
}

///
/// This function reverses a generic equation f(x,y) by computing the value
/// of "x0" when the value of the equation f(x0,y0) and y0 are known using
/// an accelerated Secant root-finder algorithm operating on the function
/// g(x,y0) = f(x,y0) - f0, where f0 = f(x0,y0) -> g(x0,y0) = 0.
///
/// - Arguments:
///   - `f` - The function to invert.
///   - `f0` - The value of f0 = f(x0,y0).
///   - `y0` - The value y0.
///   - `guess1` - A guess for the value of x0.
///   - `guess2` - A second guess for the value of x0.
///   - `tolerance` - The tolerance to utilise for the numerical method.
///   - `max_iterations` - The maximum number of allowed iterations.
///   - `n_order` - The interpolation order used by the method.
///
/// - Returns:
///   - The value of "x0" when f(x0,y0) is known.
///
#[allow(clippy::too_many_arguments)]
pub fn function_inverter_x<T: Float>(
    f: &dyn Fn(T, T) -> T,
    f0: T,
    y0: T,
    guess1: T,
    guess2: T,
    tolerance: T,
    max_iterations: i32,
    n_order: i32,
) -> T {
    if n_order < 0 {
        eprintln!("\nError in Function: function_inverter_x");
        eprintln!("The order of interpolation must be a positive integer.");
        eprintln!();
    }
    ///
    /// The solution context for the function inverter.
    ///
    struct SolutionContext<'a, T: Float> {
        ///
        /// The function to invert.
        ///
        f: &'a dyn Fn(T, T) -> T,
        ///
        /// The value of f0 = f(x0,y0).
        ///
        f0: T,
        ///
        /// The value y0.
        ///
        y0: T,
        ///
        /// A guess for the value of x0.
        ///     
        guess1: T,
        ///
        /// A second guess for the value of x0.
        ///     
        guess2: T,
        ///
        /// The interpolation order used by the method.
        ///
        n_order: i32,
    }
    //
    // Implement the solution context.
    //
    impl<'a, T: Float> SolutionContext<'a, T> {
        ///
        /// Create a new solution context.
        ///
        /// - Arguments:
        ///   - `f`: The function to invert.
        ///   - `f0`: The value of f0 = f(x0,y0).
        ///   - `y0`: The value y0.
        ///   - `guess1`: A guess for the value of x0.
        ///   - `guess2`: A second guess for the value of x0.
        ///   - `n_order`: The interpolation order used by the method.
        ///
        /// - Returns:
        ///   - The solution context.
        ///
        fn new(f: &'a dyn Fn(T, T) -> T, f0: T, y0: T, guess1: T, guess2: T, n_order: i32) -> Self {
            Self {
                f,
                f0,
                y0,
                guess1,
                guess2,
                n_order,
            }
        }
        ///
        /// Compute the maximum value for the interpolation.
        ///
        /// - Arguments:
        ///   - `p`: The value of p.
        ///
        /// - Returns:
        ///   - The maximum value for the interpolation.
        ///
        fn imax(&self, p: i32) -> i32 {
            if p == -1 || p == 0 {
                0
            } else if p <= self.n_order {
                p - 1
            } else {
                self.n_order
            }
        }
        ///
        /// Compute the solution for the function inverter.
        ///
        /// - Arguments:
        ///   - `p`: The value of p.
        ///   - `i`: The value of i.
        ///
        /// - Returns:
        ///   - The solution for the function inverter.
        ///
        fn solution(&mut self, p: i32, i: i32) -> T {
            if p == -1 {
                self.guess1
            } else if p == 0 {
                self.guess2
            } else if i == 0 {
                //
                // Secant step.
                //
                let x1 = self.solution(p - 1, self.imax(p - 1));
                let x2 = self.solution(p - 2, self.imax(p - 2));
                let f1 = (self.f)(x1, self.y0) - self.f0;
                let f2 = (self.f)(x2, self.y0) - self.f0;
                x1 - f1 * (x1 - x2) / (f1 - f2)
            } else {
                //
                // i-th order approximation.
                //
                let x1 = self.solution(p - 1, i - 1);
                let x2 = self.solution(p - 1, self.imax(p - 1));
                let x3 = self.solution(p, i - 1);
                let x4 = self.solution(p - i - 2, self.imax(p - i - 2));
                x3 + (x2 - x3) * (x1 - x3) / (x1 + x2 - x3 - x4)
            }
        }
    }

    let mut context = SolutionContext::new(f, f0, y0, guess1, guess2, n_order);
    let mut s = guess2;

    for iter in (n_order + 1)..=max_iterations {
        //
        // Compute the new value.
        //
        s = context.solution(iter, n_order);
        let error = (f(s, y0) - f0).abs();
        //
        // Check for convergence.
        //
        if error <= tolerance {
            break;
        }
    }

    s
}

///
/// This function reverses a generic equation f(x,y) by computing the value
/// of "y0" when the value of the equation f(x0,y0) and x0 are known using
/// the numerical rootfinder algorithm of Newton-Raphson operating on the
/// function g(x0,y) = f(x0,y) - f0, where f0 = f(x0,y0) -> g(x0,y0) = 0.
///
/// - Arguments:
///   - `f`: The function to invert.
///   - `f0`: The value of f0 = f(x0,y0).
///   - `x0`: The value x0.
///   - `guess1`: A guess for the value of y0.
///   - `guess2`: A second guess for the value of y0.
///   - `tolerance`: The tolerance to utilise for the numerical method.
///   - `max_iterations`: The maximum number of allowed iterations.
///   - `n_order`: The interpolation order used by the method
///
/// - Returns:
///   - The value of "y0" when f(x0,y0) is known.
///
#[allow(clippy::too_many_arguments)]
pub fn function_inverter_y<T: Float>(
    f: &dyn Fn(T, T) -> T,
    f0: T,
    x0: T,
    guess1: T,
    guess2: T,
    tolerance: T,
    max_iterations: i32,
    n_order: i32,
) -> T {
    if n_order < 0 {
        eprintln!("\nError in Function: function_inverter_y");
        eprintln!("The order of interpolation must be a positive integer.");
        eprintln!();
    }
    ///
    /// The solution context for the function inverter.
    ///
    struct SolutionContext<'a, T: Float> {
        ///
        /// The function to invert.
        ///
        f: &'a dyn Fn(T, T) -> T,
        ///
        /// The value of f0 = f(x0,y0).
        ///
        f0: T,
        ///
        /// The value x0.
        ///
        x0: T,
        ///
        /// A guess for the value of y0.
        ///     
        guess1: T,
        ///
        /// A second guess for the value of y0.
        ///
        guess2: T,
        ///
        /// The interpolation order used by the method.
        ///
        n_order: i32,
    }
    //
    // Implement the solution context.
    //
    impl<'a, T: Float> SolutionContext<'a, T> {
        ///
        /// Create a new solution context.
        ///
        /// - Arguments:
        ///   - `f`: The function to invert.
        ///   - `f0`: The value of f0 = f(x0,y0).
        ///   - `x0`: The value x0.
        ///   - `guess1`: A guess for the value of y0.
        ///   - `guess2`: A second guess for the value of y0.
        ///   - `n_order`: The interpolation order used by the method.
        ///
        fn new(f: &'a dyn Fn(T, T) -> T, f0: T, x0: T, guess1: T, guess2: T, n_order: i32) -> Self {
            Self {
                f,
                f0,
                x0,
                guess1,
                guess2,
                n_order,
            }
        }
        ///
        /// Compute the maximum value for the interpolation.
        ///
        /// - Arguments:
        ///   - `p`: The value of p.
        ///
        /// - Returns:
        ///   - The maximum value for the interpolation.
        ///
        fn imax(&self, p: i32) -> i32 {
            if p == -1 || p == 0 {
                0
            } else if p <= self.n_order {
                p - 1
            } else {
                self.n_order
            }
        }
        ///
        /// Compute the solution for the function inverter.
        ///
        /// - Arguments:
        ///   - `p`: The value of p.
        ///   - `i`: The value of i.
        ///
        /// - Returns:
        ///   - The solution for the function inverter.
        ///
        fn solution(&mut self, p: i32, i: i32) -> T {
            if p == -1 {
                self.guess1
            } else if p == 0 {
                self.guess2
            } else if i == 0 {
                //
                // Secant step.
                //
                let y1 = self.solution(p - 1, self.imax(p - 1));
                let y2 = self.solution(p - 2, self.imax(p - 2));
                let f1 = (self.f)(self.x0, y1) - self.f0;
                let f2 = (self.f)(self.x0, y2) - self.f0;
                y1 - f1 * (y1 - y2) / (f1 - f2)
            } else {
                //
                // i-th order approximation.
                //
                let y1 = self.solution(p - 1, i - 1);
                let y2 = self.solution(p - 1, self.imax(p - 1));
                let y3 = self.solution(p, i - 1);
                let y4 = self.solution(p - i - 2, self.imax(p - i - 2));
                y3 + (y2 - y3) * (y1 - y3) / (y1 + y2 - y3 - y4)
            }
        }
    }

    let mut context = SolutionContext::new(f, f0, x0, guess1, guess2, n_order);
    let mut s = guess2;

    for iter in (n_order + 1)..=max_iterations {
        //
        // Compute the new value.
        //
        s = context.solution(iter, n_order);
        let error = (f(x0, s) - f0).abs();
        //
        // Check for convergence.
        //
        if error <= tolerance {
            break;
        }
    }

    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn test_function_inverter_linear() {
        //
        // Test inverting a linear function f(x) = 2x + 3.
        //
        let f = |x: f64| -> f64 { 2.0 * x + 3.0 };
        let f0 = 7.0; // f(2) = 7
        let guess1 = 1.0;
        let guess2 = 3.0;
        let tolerance = 1e-5;
        let max_iterations = 100;
        let n_order = 1;

        let result = function_inverter(&f, f0, guess1, guess2, tolerance, max_iterations, n_order);
        //
        // Expected: x = 2 (where f(x) = 7).
        //
        assert_approx_eq!(result, 2.0, 1e-9);
    }

    #[test]
    fn test_function_inverter_quadratic() {
        //
        // Test inverting a quadratic function f(x) = x² - 4.
        //
        let f = |x: f64| -> f64 { x.powi(2) - 4.0 };
        let f0 = 0.0; // f(2) = 0 or f(-2) = 0
        let guess1 = 1.0;
        let guess2 = 3.0;
        let tolerance = 1e-10;
        let max_iterations = 100;
        let n_order = 2;

        let result = function_inverter(&f, f0, guess1, guess2, tolerance, max_iterations, n_order);
        //
        // Expected: x = 2 (because of our guesses being positive).
        //
        assert_approx_eq!(result, 2.0, 1e-9);
        //
        // Test with negative guesses to find the other root.
        //
        let result = function_inverter(&f, f0, -1.0, -3.0, tolerance, max_iterations, n_order);
        //
        // Expected: x = -2 (because of our guesses being negative).
        //
        assert_approx_eq!(result, -2.0, 1e-9);
    }

    #[test]
    fn test_function_inverter_exponential() {
        //
        // Test inverting an exponential function f(x) = e^x.
        //
        let f = |x: f64| -> f64 { x.exp() };
        let f0 = 1.0; // f(0) = 1
        let guess1 = -1.0;
        let guess2 = 1.0;
        let tolerance = 1e-10;
        let max_iterations = 100;
        let n_order = 2;

        let result = function_inverter(&f, f0, guess1, guess2, tolerance, max_iterations, n_order);
        //
        // Expected: x = 0.
        //
        assert_approx_eq!(result, 0.0, 1e-9);
    }

    #[test]
    fn test_function_inverter_x() {
        //
        // Test inverting a function f(x,y) = x² + y² for x.
        //
        let f = |x: f64, y: f64| -> f64 { x.powi(2) + y.powi(2) };
        let f0 = 5.0; // f(2,1) = 5
        let y0 = 1.0;
        let guess1 = 1.0;
        let guess2 = 3.0;
        let tolerance = 1e-10;
        let max_iterations = 100;
        let n_order = 2;

        let result = function_inverter_x(
            &f,
            f0,
            y0,
            guess1,
            guess2,
            tolerance,
            max_iterations,
            n_order,
        );
        //
        // Expected: x = 2.
        //
        assert_approx_eq!(result, 2.0, 1e-9);
        //
        // Test with negative guesses to find the other root.
        //
        let result =
            function_inverter_x(&f, f0, y0, -1.0, -3.0, tolerance, max_iterations, n_order);
        //
        // Expected: x = -2.
        //
        assert_approx_eq!(result, -2.0, 1e-9);
    }

    #[test]
    fn test_function_inverter_y() {
        //
        // Test inverting a function f(x,y) = x² + y² for y.
        //
        let f = |x: f64, y: f64| -> f64 { x.powi(2) + y.powi(2) };
        let f0 = 5.0; // f(2,1) = 5
        let x0 = 2.0;
        let guess1 = 0.5;
        let guess2 = 1.5;
        let tolerance = 1e-10;
        let max_iterations = 100;
        let n_order = 2;

        let result = function_inverter_y(
            &f,
            f0,
            x0,
            guess1,
            guess2,
            tolerance,
            max_iterations,
            n_order,
        );
        //
        // Expected: y = 1.
        //
        assert_approx_eq!(result, 1.0, 1e-9);
        //
        // Test with negative guesses to find the other root.
        //
        let result =
            function_inverter_y(&f, f0, x0, -0.5, -1.5, tolerance, max_iterations, n_order);
        //
        // Expected: y = -1.
        //
        assert_approx_eq!(result, -1.0, 1e-9);
    }

    #[test]
    fn test_function_inverter_polynomial() {
        //
        // Test inverting a polynomial function f(x) = x³ - 2x² + x - 3.
        //
        let f = |x: f64| -> f64 { x.powi(3) - 2.0 * x.powi(2) + x - 3.0 };
        let f0 = 0.0; // Looking for a root
        let guess1 = 2.0;
        let guess2 = 4.0;
        let tolerance = 1e-10;
        let max_iterations = 100;
        let n_order = 3;

        let result = function_inverter(&f, f0, guess1, guess2, tolerance, max_iterations, n_order);
        //
        // Verify that the result is actually a root.
        //
        let function_value = f(result);
        assert_approx_eq!(function_value, 0.0, 1e-7);
    }

    #[test]
    fn test_function_inverter_with_float32() {
        //
        // Test with f32 type.
        //
        let f = |x: f32| -> f32 { 2.0 * x + 3.0 };
        let f0 = 7.0;
        let guess1 = 1.0;
        let guess2 = 3.0;
        let tolerance = 1e-5; // Lower precision for f32
        let max_iterations = 100;
        let n_order = 1;

        let result = function_inverter(&f, f0, guess1, guess2, tolerance, max_iterations, n_order);

        assert_approx_eq!(result, 2.0, 1e-4);
    }

    #[test]
    fn test_function_inverter_high_order() {
        //
        // Test with a higher order of interpolation.
        //
        let f = |x: f64| -> f64 { x.powi(4) - 16.0 };
        let f0 = 0.0;
        let guess1 = 1.0;
        let guess2 = 3.0;
        let tolerance = 1e-10;
        let max_iterations = 100;
        let n_order = 4; // Higher order

        let result = function_inverter(&f, f0, guess1, guess2, tolerance, max_iterations, n_order);

        assert_approx_eq!(result, 2.0, 1e-9);
    }

    #[test]
    fn test_max_iterations_limit() {
        //
        // Test that the function returns a reasonable approximation
        // even when max iterations is small.
        //
        let f = |x: f64| -> f64 { x.sqrt() };
        let f0 = 3.0; // f(9) = 3
        let guess1 = 4.0;
        let guess2 = 12.0;
        let tolerance = 1e-10;
        let max_iterations = 5; // Very few iterations
        let n_order = 2;

        let result = function_inverter(&f, f0, guess1, guess2, tolerance, max_iterations, n_order);
        //
        // Result might not be exact but should be close.
        //
        let function_value = f(result);
        assert!((function_value - f0).abs() < 0.1);
    }

    #[test]
    fn test_invalid_order() {
        //
        // Test with invalid order (negative).
        //
        let f = |x: f64| -> f64 { x };
        let f0 = 1.0;
        let guess1 = 0.0;
        let guess2 = 2.0;
        let tolerance = 1e-10;
        let max_iterations = 100;
        let n_order = -1; // Invalid

        let result = function_inverter(&f, f0, guess1, guess2, tolerance, max_iterations, n_order);
        //
        // Should return NaN for invalid order.
        //
        assert!(result.is_nan());
    }
}
