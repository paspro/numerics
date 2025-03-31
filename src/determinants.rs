// -------------------------------------------------------------------------------------------------
//
//  This library implements various numerical algorithms which can be of use for all
//  kinds of scientific and engineering applications.
//
//  Copyright (c) 2025 by Dr. Panos Asproulis (p.asproulis@icloud.com).
//  All Rights Reserved.
//
// -------------------------------------------------------------------------------------------------

//!
//! This submodule contains functions which compute determinants of matrices.
//!

use num_traits::Float;

///
/// Compute the determinant of a 4x4 matrix.
///
/// - Arguments:
///   - `a`: The 4x4 array to use as an `Array2<f64>` object.
///
/// - Returns:
///   - The value of the determinant.
///
#[inline]
pub fn matrix_determinant_4x4<T: Float>(a: &[[T; 4]]) -> T {
    a[0][0]
        * (a[1][1] * (a[2][2] * a[3][3] - a[2][3] * a[3][2])
            - a[1][2] * (a[2][1] * a[3][3] - a[2][3] * a[3][1])
            + a[1][3] * (a[2][1] * a[3][2] - a[2][2] * a[3][1]))
        - a[0][1]
            * (a[1][0] * (a[2][2] * a[3][3] - a[2][3] * a[3][2])
                - a[1][2] * (a[2][0] * a[3][3] - a[2][3] * a[3][0])
                + a[1][3] * (a[2][0] * a[3][2] - a[2][2] * a[3][0]))
        + a[0][2]
            * (a[1][0] * (a[2][1] * a[3][3] - a[2][3] * a[3][1])
                - a[1][1] * (a[2][0] * a[3][3] - a[2][3] * a[3][0])
                + a[1][3] * (a[2][0] * a[3][1] - a[2][1] * a[3][0]))
        - a[0][3]
            * (a[1][0] * (a[2][1] * a[3][2] - a[2][2] * a[3][1])
                - a[1][1] * (a[2][0] * a[3][2] - a[2][2] * a[3][0])
                + a[1][2] * (a[2][0] * a[3][1] - a[2][1] * a[3][0]))
}

// -------------------------------------------------------------------------------------------------
//
// Unit Tests.
//
// -------------------------------------------------------------------------------------------------

#[cfg(test)]
pub mod tests {
    use super::*;

    #[test]
    fn test_matrix_determinant_4x4_identity() {
        let identity = vec![
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        assert_eq!(matrix_determinant_4x4(&identity), 1.0);
    }

    #[test]
    fn test_matrix_determinant_4x4_zero() {
        let zero = vec![
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ];
        assert!(f64::abs(matrix_determinant_4x4(&zero)) <= f64::EPSILON);
    }

    #[test]
    fn test_matrix_determinant_4x4_diagonal() {
        let diagonal = vec![
            [2.0, 0.0, 0.0, 0.0],
            [0.0, 3.0, 0.0, 0.0],
            [0.0, 0.0, 4.0, 0.0],
            [0.0, 0.0, 0.0, 5.0],
        ];
        assert_eq!(matrix_determinant_4x4(&diagonal), 120.0); // 2 * 3 * 4 * 5
    }

    #[test]
    fn test_matrix_determinant_4x4_triangular() {
        let upper_triangular = vec![
            [1.0, 2.0, 3.0, 4.0],
            [0.0, 2.0, 3.0, 4.0],
            [0.0, 0.0, 3.0, 4.0],
            [0.0, 0.0, 0.0, 4.0],
        ];
        assert_eq!(matrix_determinant_4x4(&upper_triangular), 24.0); // 1 * 2 * 3 * 4
    }

    #[test]
    fn test_matrix_determinant_4x4_singular() {
        //
        // Matrix with linearly dependent rows
        //
        let singular = vec![
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 6.0, 8.0],   // Multiple of first row
            [3.0, 6.0, 9.0, 12.0],  // Multiple of first row
            [4.0, 8.0, 12.0, 16.0], // Multiple of first row
        ];
        assert_eq!(matrix_determinant_4x4(&singular), 0.0);
    }

    #[test]
    fn test_matrix_determinant_4x4_skew_symmetric() {
        let skew_symmetric = vec![
            [0.0, -1.0, 2.0, -3.0],
            [1.0, 0.0, -4.0, 5.0],
            [-2.0, 4.0, 0.0, -6.0],
            [3.0, -5.0, 6.0, 0.0],
        ];
        //
        // Skew-symmetric matrices of even order have determinant ≥ 0
        //
        assert!(matrix_determinant_4x4(&skew_symmetric) >= 0.0);
    }

    #[test]
    fn test_matrix_determinant_4x4_float_precision() {
        let matrix = vec![
            [1e-15, 0.0, 0.0, 0.0],
            [0.0, 1e-15, 0.0, 0.0],
            [0.0, 0.0, 1e-15, 0.0],
            [0.0, 0.0, 0.0, 1e-15],
        ];
        let det = matrix_determinant_4x4(&matrix);
        assert!((det - 1e-60).abs() < 1e-75);
    }

    #[test]
    fn test_matrix_determinant_4x4_large_values() {
        let matrix = vec![
            [1e5, 0.0, 0.0, 0.0],
            [0.0, 1e5, 0.0, 0.0],
            [0.0, 0.0, 1e5, 0.0],
            [0.0, 0.0, 0.0, 1e5],
        ];
        assert_eq!(matrix_determinant_4x4(&matrix), 1e20);
    }
}
