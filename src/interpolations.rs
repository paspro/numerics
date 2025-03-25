// -------------------------------------------------------------------------------------------------
//
//  This library implements various numerical algorithms which can be of use for all
//  kinds of scientific and engineering applications.
//
//  Copyright (c) 2025 by Dr. Panos Asproulis (p.asproulis@icloud.com).
//  All Rights Reserved.
//
// -------------------------------------------------------------------------------------------------

//!  This module implements various interpolation algorithms.

use num_traits::Float;

//
/// This function performs a piecewise cubic interpolation in order to
/// compute the value of a property. The interpolation function has the
/// form:
///
/// s(x) = SUM(Cmj*(x-xj)^m), where the sum is over m = 0, 3
///        and Cmj are the interpolation coefficients computed based on
///        the supplied data in the form of vectors x, u(x), du(x)/dx.
///
/// - Arguments:
///   - `x`: The vector with the x-values (size=2).
///   - `u`: The vector with the property values (size=2).
///   - `uprime`: The vector with the gradient of the property.
///   - `xi`: The x-value to be used for computing the interpolated
///           property value which lies between x0 and x1.
///
/// - Returns:
///   - The interpolated value at location xi.
///
pub fn interpolate_cubic<T: Float>(x: &[T; 2], u: &[T; 2], u_prime: &[T; 2], xi: T) -> T {
    //
    // Calculate key differences.
    //
    let dx = x[1] - x[0]; // x₁ - x₀
    let t = (xi - x[0]) / dx; // Normalized coordinate (0 to 1)
    //
    // Constants for calculation.
    //
    let one = T::one();
    let two = T::from(2.0).unwrap();
    let three = T::from(3.0).unwrap();
    //
    // Hermite basis functions for cubic interpolation.
    //
    let h00 = two * t.powi(3) - three * t.powi(2) + one; // 2t³ - 3t² + 1
    let h10 = t.powi(3) - two * t.powi(2) + t; // t³ - 2t² + t
    let h01 = -two * t.powi(3) + three * t.powi(2); // -2t³ + 3t²
    let h11 = t.powi(3) - t.powi(2); // t³ - t²
    //
    // Compute interpolated value using Hermite basis functions.
    //
    h00 * u[0] + h10 * dx * u_prime[0] + h01 * u[1] + h11 * dx * u_prime[1]
}

///
/// This function performs a cubic interpolation in order to compute
/// the value of a property.
///
/// - Arguments:
///   - `x`: The vector with the x-values (size=4).
///   - `u`: The vector with the property values (size=4).
///   - `xi`: The x-value to be used for computing the interpolated
///           property value which lies between x1 and x2.
///
/// - Returns:
///   - The interpolated value at location xi.
///
pub fn interpolate_cubic_numerical<T: Float>(x: &[T; 4], u: &[T; 4], xi: T) -> T {
    //
    // Compute the interpolation coefficients.
    //
    let dx = [xi - x[0], xi - x[1], xi - x[2], xi - x[3]];

    let div = [
        (x[0] - x[1]) * (x[0] - x[2]) * (x[0] - x[3]),
        (x[1] - x[0]) * (x[1] - x[2]) * (x[1] - x[3]),
        (x[2] - x[0]) * (x[2] - x[1]) * (x[2] - x[3]),
        (x[3] - x[0]) * (x[3] - x[1]) * (x[3] - x[2]),
    ];

    let coeff = [
        dx[1] * dx[2] * dx[3] / div[0],
        dx[0] * dx[2] * dx[3] / div[1],
        dx[0] * dx[1] * dx[3] / div[2],
        dx[0] * dx[1] * dx[2] / div[3],
    ];
    //
    // Perform the interpolation.
    //
    u[0] * coeff[0] + u[1] * coeff[1] + u[2] * coeff[2] + u[3] * coeff[3]
}

///
/// This function performs a piecewise bilinear interpolation in order to
/// compute the value of a property that is a function of 2.0 parameters.
/// The interpolation function has the form:
///
/// s(x,y) = SUMn(SUMm(Cmnjk*(x-xj)^m*(y-yk)^n), where SUMm is over m = 0, 1,
///         SUMn is over n = 0, 1 and Cmnjk are the interpolation coefficients
///         computed based on the supplied data in the form of the table x, y,
///         u(x,y).
///
/// - Arguments:
///   - `x`: The vector with the x-values (dim=2).
///   - `y`: The vector with the y-values (dim=2).
///   - `u`: The array with the property values (dim=2x2).
///   - `xi`: The x-value to be used for computing the interpolated property value.
///   - `yi`: The y-value to be used for computing the interpolated property value.
///
/// - Returns:
///   - The interpolated value at location xi, yi.
///
pub fn interpolate_bilinear<T: Float>(x: &[T; 2], y: &[T; 2], u: &[[T; 2]; 2], xi: T, yi: T) -> T {
    //
    // Compute the interpolation coefficients.
    //
    let dx21 = x[1] - x[0];
    let dy21 = y[1] - y[0];
    let dx2i = x[1] - xi;
    let dy2i = y[1] - yi;
    let dxi1 = xi - x[0];
    let dyi1 = yi - y[0];
    //
    // Perform the interpolation.
    //
    let u1d = (u[0][0] * dx2i + u[1][0] * dxi1) * dy2i;
    let u2d = (u[0][1] * dx2i + u[1][1] * dxi1) * dyi1;

    (u1d + u2d) / (dx21 * dy21)
}

///
/// This function performs a bicubic interpolation in order to
/// compute the value of a property that is a function of 2.0 parameters.
///
/// - Arguments:
///   - `x`: The vector with the x-values (dim=4).
///   - `y`: The vector with the y-values (dim=4).
///   - `u`: The array with the property values (dim=4x4).
///   - `xi`: The x-value to be used for computing the interpolated property value.
///   - `yi`: The y-value to be used for computing the interpolated property value.
///
/// - Returns:
///   - The interpolated value at location xi, yi.
///
pub fn interpolate_bicubic<T: Float>(x: &[T; 4], y: &[T; 4], u: &[[T; 4]; 4], xi: T, yi: T) -> T {
    //
    // Perform a cubic interpolation for each column.
    //
    let mut uv = [T::zero(); 4];

    for i in 0..4 {
        let column = [u[i][0], u[i][1], u[i][2], u[i][3]];
        uv[i] = interpolate_cubic_numerical(y, &column, yi);
    }
    //
    // Perform a horizontal interpolation.
    //
    interpolate_cubic_numerical(x, &uv, xi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn test_interpolate_cubic() {
        //
        // Test case 1: Interpolate a linear function at midpoint.
        //
        let x = [0.0, 1.0];
        let u = [1.0, 3.0];
        let u_prime = [2.0, 2.0]; // Constant gradient
        let xi = 0.5;

        let result = interpolate_cubic(&x, &u, &u_prime, xi);
        assert_approx_eq!(result, 2.0, 1e-6);
        //
        // Test case 2: Interpolate a quadratic function.
        //
        let x = [1.0, 3.0];
        let u = [1.0, 9.0]; // f(x) = x²
        let u_prime = [2.0, 6.0]; // f'(x) = 2x
        let xi = 2.0;

        let result = interpolate_cubic(&x, &u, &u_prime, xi);
        assert_approx_eq!(result, 4.0, 1e-6);
        //
        // Test case 3: Interpolate with non-uniform spacing.
        //
        let x = [1.0, 4.0];
        let u = [2.0, 8.0];
        let u_prime = [1.5, 3.0];
        let xi = 2.5;

        let result = interpolate_cubic(&x, &u, &u_prime, xi);
        println!("result: {}", result);
        assert_approx_eq!(result, 4.4375, 1e-6);
    }

    #[test]
    fn test_interpolate_cubic_numerical() {
        //
        // Test case 1: Interpolate a cubic function.
        // f(x) = x³ at points x = -1, 0, 1, 2.
        //
        let x = [-1.0, 0.0, 1.0, 2.0];
        let u = [-1.0, 0.0, 1.0, 8.0]; // f(x) = x³
        let xi = 0.5;

        let result = interpolate_cubic_numerical(&x, &u, xi);
        assert_approx_eq!(result, 0.125, 1e-6); // 0.5³ = 0.125
        //
        // Test case 2: Interpolate a quadratic function.
        // f(x) = x² at points x = -2, -1, 1, 2.
        //
        let x = [-2.0, -1.0, 1.0, 2.0];
        let u = [4.0, 1.0, 1.0, 4.0]; // f(x) = x²
        let xi = 0.0;

        let result = interpolate_cubic_numerical(&x, &u, xi);
        assert_approx_eq!(result, 0.0, 1e-6); // 0² = 0
        //
        // Test case 3: Interpolation at one of the data points.
        //
        let x = [1.0, 2.0, 3.0, 4.0];
        let u = [10.0, 20.0, 30.0, 40.0];
        let xi = 3.0; // Exact data point

        let result = interpolate_cubic_numerical(&x, &u, xi);
        assert_approx_eq!(result, 30.0, 1e-6);
    }

    #[test]
    fn test_interpolate_bilinear() {
        //
        // Test case 1: Simple bilinear interpolation.
        //
        let x = [0.0, 1.0];
        let y = [0.0, 1.0];
        let u = [
            [1.0, 2.0], // u(0,0) = 1, u(0,1) = 2
            [3.0, 4.0], // u(1,0) = 3, u(1,1) = 4
        ];
        //
        // Interpolate at center point (0.5, 0.5).
        //
        let result = interpolate_bilinear(&x, &y, &u, 0.5, 0.5);
        assert_approx_eq!(result, 2.5, 1e-6); // Should be average of all four corners
        //
        // Test case 2: Interpolate along edges.
        //
        let result_edge_x = interpolate_bilinear(&x, &y, &u, 0.5, 0.0);
        assert_approx_eq!(result_edge_x, 2.0, 1e-6); // Halfway between 1 and 3

        let result_edge_y = interpolate_bilinear(&x, &y, &u, 0.0, 0.5);
        assert_approx_eq!(result_edge_y, 1.5, 1e-6); // Halfway between 1 and 2
        //
        // Test case 3: Non-uniform grid.
        //
        let x = [1.0, 3.0];
        let y = [2.0, 5.0];
        let u = [[10.0, 20.0], [30.0, 40.0]];
        //
        // Interpolate at point (2.0, 3.5).
        //
        let result = interpolate_bilinear(&x, &y, &u, 2.0, 3.5);
        assert_approx_eq!(result, 25.0, 1e-6);
    }

    #[test]
    fn test_interpolate_bicubic() {
        //
        // Test case 1: Interpolate a simple function f(x,y) = x + y.
        //
        let x = [-1.0, 0.0, 1.0, 2.0];
        let y = [-1.0, 0.0, 1.0, 2.0];
        let u = [
            [-2.0, -1.0, 0.0, 1.0], // f(-1,y) = -1 + y
            [-1.0, 0.0, 1.0, 2.0],  // f(0,y) = 0 + y
            [0.0, 1.0, 2.0, 3.0],   // f(1,y) = 1 + y
            [1.0, 2.0, 3.0, 4.0],   // f(2,y) = 2 + y
        ];
        //
        // Interpolate at point (0.5, 0.5).
        //
        let result = interpolate_bicubic(&x, &y, &u, 0.5, 0.5);
        assert_approx_eq!(result, 1.0, 1e-6); // 0.5 + 0.5 = 1.0
        //
        // Test case 2: Interpolate at exact data point.
        //
        let result_exact = interpolate_bicubic(&x, &y, &u, 1.0, 1.0);
        assert_approx_eq!(result_exact, 2.0, 1e-6); // 1 + 1 = 2
        //
        // Test case 3: Interpolate a more complex function (quadratic).
        // f(x,y) = x² + y²
        //
        let u_quadratic = [
            [2.0, 1.0, 2.0, 5.0], // f(-1,y) = 1 + y²
            [1.0, 0.0, 1.0, 4.0], // f(0,y) = 0 + y²
            [2.0, 1.0, 2.0, 5.0], // f(1,y) = 1 + y²
            [5.0, 4.0, 5.0, 8.0], // f(2,y) = 4 + y²
        ];
        //
        // Interpolate at point (0.5, 0.5).
        //
        let result_quad = interpolate_bicubic(&x, &y, &u_quadratic, 0.5, 0.5);
        assert_approx_eq!(result_quad, 0.5, 1e-2); // 0.5² + 0.5² = 0.5 (with some tolerance)
    }

    #[test]
    fn test_float_types() {
        //
        // Test with f32.
        //
        let x_f32 = [0.0f32, 1.0f32];
        let u_f32 = [1.0f32, 3.0f32];
        let u_prime_f32 = [2.0f32, 2.0f32];
        let xi_f32 = 0.5f32;

        let result_f32 = interpolate_cubic(&x_f32, &u_f32, &u_prime_f32, xi_f32);
        assert_approx_eq!(result_f32, 2.0f32, 1e-6);
        //
        // Test with f64.
        //
        let x_f64 = [0.0f64, 1.0f64];
        let u_f64 = [1.0f64, 3.0f64];
        let u_prime_f64 = [2.0f64, 2.0f64];
        let xi_f64 = 0.5f64;

        let result_f64 = interpolate_cubic(&x_f64, &u_f64, &u_prime_f64, xi_f64);
        assert_approx_eq!(result_f64, 2.0f64, 1e-10);
    }
}
