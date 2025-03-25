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

use num_traits::real::Real;

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
pub fn binary_search_vector<T: PartialOrd>(vc: &[T], value: &T, imin: usize, imax: usize) -> usize {
    let mut high = imax;
    let mut low = imin;
    let mut med = imin + (imax - imin) / 2;

    while high != (low + 1) {
        let vmed = &vc[med];

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
/// a value in a 2-dimensional array. If the exact value does not exist
/// it will return the closest smaller value available in the table.
/// Note that the supplied matrix must be sorted in ascending order.
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
pub fn binary_search_array<T: PartialOrd + Real>(
    mat: &[Vec<T>],
    value: &T,
    imin: usize,
    imax: usize,
    jmin: usize,
    jmax: usize,
) -> (usize, usize) {
    let mut j_index = binary_search_vector(&mat[imin], value, jmin, jmax);
    let mut value_found = &mat[imin][j_index];
    let mut error_previous = *value - *value_found;

    let mut i = imin;
    let mut j = j_index;

    for i_index in (imin + 1)..=imax {
        j_index = binary_search_vector(&mat[i_index], value, jmin, jmax);
        value_found = &mat[i_index][j_index];
        let error_new = *value - *value_found;

        if error_new >= T::zero() && error_new < error_previous {
            i = i_index;
            j = j_index;
        }

        error_previous = error_new;
    }

    (i, j)
}
