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
pub fn interpolate_cubic(x: &[f64; 2], u: &[f64; 2], uprime: &[f64; 2], xi: f64) -> f64 {
    //
    // Compute the interpolation coefficients.
    //
    let dx = x[1] - x[0];
    let du = u[1] - u[0];
    let dx2 = dx * dx;
    let dx3 = dx2 * dx;

    let coeff = [
        u[0],
        uprime[0],
        3.0 * du / dx2 - (uprime[1] + 2.0 * uprime[0]) / dx,
        -2.0 * du / dx3 + (uprime[1] + uprime[0]) / dx2,
    ];
    //
    // Perform the interpolation.
    //
    let x_diff = xi - x[0];
    coeff[0] + coeff[1] * x_diff + coeff[2] * x_diff.powi(2) + coeff[3] * x_diff.powi(3)
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
pub fn interpolate_cubic_numerical(x: &[f64; 4], u: &[f64; 4], xi: f64) -> f64 {
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
pub fn interpolate_bilinear(
    x: &[f64; 2],
    y: &[f64; 2],
    u: &[[f64; 2]; 2],
    xi: f64,
    yi: f64,
) -> f64 {
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
pub fn interpolate_bicubic(x: &[f64; 4], y: &[f64; 4], u: &[[f64; 4]; 4], xi: f64, yi: f64) -> f64 {
    //
    // Perform a cubic interpolation for each column.
    //
    let mut uv = [0.0; 4];

    for i in 0..4 {
        let column = [u[i][0], u[i][1], u[i][2], u[i][3]];
        uv[i] = interpolate_cubic_numerical(y, &column, yi);
    }
    //
    // Perform a horizontal interpolation.
    //
    interpolate_cubic_numerical(x, &uv, xi)
}
