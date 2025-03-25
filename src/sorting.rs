// -------------------------------------------------------------------------------------------------
//
//  This library implements various numerical algorithms which can be of use for all
//  kinds of scientific and engineering applications.
//
//  Copyright (c) 2025 by Dr. Panos Asproulis (p.asproulis@icloud.com).
//  All Rights Reserved.
//
// -------------------------------------------------------------------------------------------------

//! Implementation of various sorting algorithms.

///
/// This function implements the QuickSort algorithm for sorting the
/// elements of a vector to ascending order.
///
/// - Arguments:
///   - `vec`: The vector to sort.
///   - `indx`: A vector of indices that represents the modification
///             of the original vector. It can be used in order to
///             apply similar modifications to other associated vectors.
///
pub fn quicksort(vec: &mut [f64], indx: &mut Option<&mut [i32]>) {
    if vec.len() > 1 {
        let iq = partition(vec, indx);
        //
        // Sort the partitions.
        //
        if let Some(idx) = indx {
            let (left, right) = idx.split_at_mut(iq);
            quicksort(&mut vec[..iq], &mut Some(left));
            quicksort(&mut vec[iq..], &mut Some(right));
        } else {
            quicksort(&mut vec[..iq], &mut None);
            quicksort(&mut vec[iq..], &mut None);
        }
    }
}

///
/// This function partitions a vector for sorting purposes.
///
/// - Arguments:
///   - `vec`: The vector to partition.
///   - `indx`: An optional vector of indices to be partitioned in the same manner.
///
/// - Returns:
///   - The partition marker.
///
fn partition(vec: &mut [f64], indx: &mut Option<&mut [i32]>) -> usize {
    let x = vec[0];
    let mut i = 0;
    let mut j = vec.len();

    loop {
        j -= 1;
        while vec[j] > x {
            j -= 1;
        }

        i += 1;
        while i < j && vec[i] < x {
            i += 1;
        }

        if i < j {
            vec.swap(i, j);
            if let Some(idx) = indx {
                idx.swap(i, j);
            }
        } else {
            return if i == j { i + 1 } else { i };
        }
    }
}
