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
                let x1 = self.solution(p - 1, self.imax(p - 1));
                let x2 = self.solution(p - 2, self.imax(p - 2));
                let f1 = (self.f)(x1) - self.f0;
                let f2 = (self.f)(x2) - self.f0;
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

    let mut context = SolutionContext::new(f, f0, guess1, guess2, n_order);
    let mut s_old = guess2;

    for iter in (n_order + 1)..=max_iterations {
        //
        // Compute the new value.
        //
        let s = context.solution(iter, n_order);
        let error = (s_old - s).abs() / s.abs();
        s_old = s;
        //
        // Check for convergence.
        //
        if error <= tolerance {
            break;
        }
    }

    s_old
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
    let mut s_old = guess2;

    for iter in (n_order + 1)..=max_iterations {
        //
        // Compute the new value.
        //
        let s = context.solution(iter, n_order);
        let error = (s_old - s).abs() / s.abs();
        s_old = s;
        //
        // Check for convergence.
        //
        if error <= tolerance {
            break;
        }
    }

    s_old
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
    let mut s_old = guess2;

    for iter in (n_order + 1)..=max_iterations {
        //
        // Compute the new value.
        //
        let s = context.solution(iter, n_order);
        let error = (s_old - s).abs() / s.abs();
        s_old = s;
        //
        // Check for convergence.
        //
        if error <= tolerance {
            break;
        }
    }

    s_old
}
