// -------------------------------------------------------------------------------------------------
//
//  This library implements various numerical algorithms which can be of use for all
//  kinds of scientific and engineering applications.
//
//  Copyright (c) 2025 by Dr. Panos Asproulis (p.asproulis@icloud.com).
//  All Rights Reserved.
//
// -------------------------------------------------------------------------------------------------

//! This module implements various search algorithms.

///
/// This function implements a binary search algorithm in order to locate
/// a value in a vector. If the exact value does not exist it will return
/// the closest smaller value available in the table. Note that the
/// supplied vector must be sorted in ascending order.
///
/// - Arguments:
///   - `vec`: The vector to use for searching.
///   - `value`: The value to search for.
///   - `imin`: The minimum vector index to use for searching.
///   - `imax`: The maximum vector index to use for searching.
///
/// - Returns:
///   - The index with the location of the result.
///
pub fn binary_search_vector(vec: &[f64], value: f64, imin: i32, imax: i32) -> i32 {
    let mut high = imax;
    let mut low = imin;
    let mut med = imin + (imax - imin) / 2;

    while high != (low + 1) {
        let vmed = vec[med as usize];

        if vmed <= value {
            low = med;
        } else {
            high = med;
        }

        med = low + (high - low) / 2;
    }

    med
}

///
/// This function implements a binary search algorithm in order to locate
/// a value in a 2.0 dimensional array. If the exact value does not exist
/// it will return the closest smaller value available in the table.
///
/// - Arguments:
///   - `mat`: The matrix (2D) to use for searching.
///   - `value`: The value to search for.
///   - `imin`: The minimum array i-index to use for searching.
///   - `imax`: The maximum array i-index to use for searching.
///   - `jmin`: The minimum array j-index to use for searching.
///   - `jmax`: The maximum array j-index to use for searching.
///
/// - Returns:
///   - A tuple (i, j) with the indices for the location of the value.
///
pub fn binary_search_array(
    mat: &[Vec<f64>],
    value: f64,
    imin: i32,
    imax: i32,
    jmin: i32,
    jmax: i32,
) -> (i32, i32) {
    let mut j_index = binary_search_vector(&mat[imin as usize], value, jmin, jmax);
    let mut value_found = mat[imin as usize][j_index as usize];
    let mut error_previous = (value - value_found) / value_found.abs();

    let mut i = imin;
    let mut j = j_index;

    for i_index in (imin + 1)..=imax {
        j_index = binary_search_vector(&mat[i_index as usize], value, jmin, jmax);
        value_found = mat[i as usize][j as usize];
        let error_new = (value - value_found) / value_found.abs();

        if error_new >= 0.0 && error_new < error_previous {
            i = i_index;
            j = j_index;
        }

        error_previous = error_new;
    }

    (i, j)
}
