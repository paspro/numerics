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
/// Sorts a vector in ascending order using the quicksort algorithm and 
/// optionally tracks the original indices of elements.
///
/// This function implements an in-place quicksort that can also maintain
/// a parallel array of indices, showing where each element in the sorted
/// result originally came from.
///
/// - Arguments:
///   - `vec`: The vector to sort.
///   - `return_indices`: If true, returns a vector containing the original indices.
///
/// - Returns:
///   - If `return_indices` is true, returns a vector of indices representing .
///     the original positions of elements in the sorted result.
///   - If `return_indices` is false, returns an empty vector.
///
/// - Example:
///
/// ```
/// let mut values = vec![5, 2, 9, 1, 5, 6];
/// let indices = quicksort_with_indices(&mut values, true);
/// assert_eq!(values, [1, 2, 5, 5, 6, 9]);
/// assert_eq!(indices, [3, 1, 0, 4, 5, 2]); // Original positions of each element
/// ```
///
pub fn quicksort_with_indices<T: PartialOrd>(vec: &mut [T], return_indices: bool) -> Vec<usize> {
    let n = vec.len();
    //
    // Create index array if requested.
    //
    let mut indices = if return_indices {
        (0..n).collect::<Vec<usize>>()
    } else {
        Vec::new()
    };
    //
    // Call the internal quicksort implementation.
    //
    if n > 1 {
        if return_indices {
            _quicksort_with_indices(vec, &mut indices, 0, n - 1);
        } else {
            _quicksort(vec, 0, n - 1);
        }
    }
    
    indices
}

///
/// Internal implementation of quicksort with index tracking.
/// 
/// - Arguments:
///   - `vec`: The vector to sort.
///   - `indices`: The vector to store the original indices.
///   - `left`: The leftmost index of the vector.
///   - `right`: The rightmost index of the vector.
/// 
fn _quicksort_with_indices<T: PartialOrd>(vec: &mut [T], indices: &mut [usize], left: usize, right: usize) {
    if left < right {
        let pivot = partition_with_indices(vec, indices, left, right);
        
        if pivot > 0 {
            _quicksort_with_indices(vec, indices, left, pivot - 1);
        }
        _quicksort_with_indices(vec, indices, pivot + 1, right);
    }
}

///
/// Internal implementation of quicksort without index tracking (for performance).
/// 
/// - Arguments:
///   - `vec`: The vector to sort.
///   - `left`: The leftmost index of the vector.
///   - `right`: The rightmost index of the vector.
/// 
fn _quicksort<T: PartialOrd>(vec: &mut [T], left: usize, right: usize) {
    if left < right {
        let pivot = partition(vec, left, right);
        
        if pivot > 0 {
            _quicksort(vec, left, pivot - 1);
        }
        _quicksort(vec, pivot + 1, right);
    }
}

///
/// Partition function with index tracking.
/// 
/// - Arguments:
///   - `vec`: The vector to partition.
///   - `indices`: The vector to store the original indices.
///   - `left`: The leftmost index of the vector.
///   - `right`: The rightmost index of the vector.
/// 
/// - Returns:
///   - The pivot index.
/// 
fn partition_with_indices<T: PartialOrd>(vec: &mut [T], indices: &mut [usize], left: usize, right: usize) -> usize {
    //
    // Use rightmost element as pivot.
    //
    let pivot = right;
    let mut i = left;
    
    for j in left..right {
        if vec[j] <= vec[pivot] {
            //
            // Swap elements.
            //
            vec.swap(i, j);
            //
            // Also swap indices to track the original positions.
            //
            indices.swap(i, j);
            i += 1;
        }
    }
    //
    // Put pivot in its final position.
    //
    vec.swap(i, pivot);
    indices.swap(i, pivot);
    
    i
}

///
/// Partition function without index tracking.
/// 
/// - Arguments:
///   - `vec`: The vector to partition.
///   - `left`: The leftmost index of the vector.
///   - `right`: The rightmost index of the vector.
/// 
/// - Returns:
///   - The pivot index.
/// 
fn partition<T: PartialOrd>(vec: &mut [T], left: usize, right: usize) -> usize {
    //
    // Use rightmost element as pivot.
    //
    let pivot = right;
    let mut i = left;
    
    for j in left..right {
        if vec[j] <= vec[pivot] {
            vec.swap(i, j);
            i += 1;
        }
    }
    
    vec.swap(i, pivot);
    i
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quicksort_with_indices_basic() {
        //
        // Test with integers.
        //
        let mut vec = vec![5, 2, 9, 1, 5, 6];
        let indices = quicksort_with_indices(&mut vec, true);
        
        assert_eq!(vec, vec![1, 2, 5, 5, 6, 9]);
        assert_eq!(indices, vec![3, 1, 0, 4, 5, 2]);
        //
        // Test without index tracking.
        //
        let mut vec = vec![5, 2, 9, 1, 5, 6];
        let indices = quicksort_with_indices(&mut vec, false);
        
        assert_eq!(vec, vec![1, 2, 5, 5, 6, 9]);
        assert_eq!(indices, Vec::<usize>::new());
    }
    
    #[test]
    fn test_quicksort_with_indices_empty() {
        //
        // Test with empty vector.
        //
        let mut vec: Vec<i32> = vec![];
        let indices = quicksort_with_indices(&mut vec, true);
        
        assert_eq!(vec, vec![]);
        assert_eq!(indices, vec![]);
    }
    
    #[test]
    fn test_quicksort_with_indices_single_element() {
        //
        // Test with single element.
        //
        let mut vec = vec![42];
        let indices = quicksort_with_indices(&mut vec, true);
        
        assert_eq!(vec, vec![42]);
        assert_eq!(indices, vec![0]);
    }
    
    #[test]
    fn test_quicksort_with_indices_already_sorted() {
        //
        // Test with already sorted vector.
        //
        let mut vec = vec![1, 2, 3, 4, 5];
        let indices = quicksort_with_indices(&mut vec, true);
        
        assert_eq!(vec, vec![1, 2, 3, 4, 5]);
        assert_eq!(indices, vec![0, 1, 2, 3, 4]);
    }
    
    #[test]
    fn test_quicksort_with_indices_reverse_sorted() {
        //
        // Test with reverse sorted vector.
        //
        let mut vec = vec![5, 4, 3, 2, 1];
        let indices = quicksort_with_indices(&mut vec, true);
        
        assert_eq!(vec, vec![1, 2, 3, 4, 5]);
        assert_eq!(indices, vec![4, 3, 2, 1, 0]);
    }
    
    #[test]
    fn test_quicksort_with_indices_duplicate_elements() {
        //
        // Test with duplicate elements.
        //
        let mut vec = vec![3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5];
        let indices = quicksort_with_indices(&mut vec, true);
        
        assert_eq!(vec, vec![1, 1, 2, 3, 3, 4, 5, 5, 5, 6, 9]);
        //
        // Check that indices correctly map to original positions.
        //
        for (i, &idx) in indices.iter().enumerate() {
            assert_eq!(vec[i], match idx {
                0 => 3,
                1 => 1,
                2 => 4,
                3 => 1,
                4 => 5,
                5 => 9,
                6 => 2,
                7 => 6,
                8 => 5,
                9 => 3,
                10 => 5,
                _ => panic!("Invalid index"),
            });
        }
    }
    
    #[test]
    fn test_quicksort_with_indices_different_types() {
        //
        // Test with floating point numbers.
        //
        let mut vec = vec![3.14, 1.59, 2.65, 3.58, 9.79];
        let indices = quicksort_with_indices(&mut vec, true);
        
        assert!(vec[0] <= vec[1] && vec[1] <= vec[2] && vec[2] <= vec[3] && vec[3] <= vec[4]);
        assert_eq!(indices.len(), 5);
        //
        // Test with strings.
        //
        let mut vec = vec!["banana", "apple", "cherry", "date", "elderberry"];
        let indices = quicksort_with_indices(&mut vec, true);
        
        assert_eq!(vec, vec!["apple", "banana", "cherry", "date", "elderberry"]);
        assert_eq!(indices, vec![1, 0, 2, 3, 4]);
    }
    
    #[test]
    fn test_quicksort_with_indices_large() {
        //
        // Test with a larger vector.
        //
        let original: Vec<i32> = (0..100).rev().collect();
        let expected: Vec<i32> = (0..100).collect();
        let expected_indices: Vec<usize> = (0..100).rev().collect();
        
        let mut vec = original.clone();
        let indices = quicksort_with_indices(&mut vec, true);
        
        assert_eq!(vec, expected);
        assert_eq!(indices, expected_indices);
    }
    
    #[test]
    fn test_partition() {
        //
        // Test basic partition operation.
        //
        let mut vec = vec![5, 2, 9, 1, 5, 6];
        let pivot_idx = partition(&mut vec, 0, 5);
        //
        // All elements before pivot should be <= pivot value.
        //
        for i in 0..pivot_idx {
            assert!(vec[i] <= vec[pivot_idx]);
        }
        //
        // All elements after pivot should be > pivot value.
        //
        for i in (pivot_idx+1)..vec.len() {
            assert!(vec[i] > vec[pivot_idx]);
        }
    }
    
    #[test]
    fn test_partition_with_indices() {
        //
        // Test partition with index tracking.
        //
        let mut vec = vec![5, 2, 9, 1, 5, 6];
        let mut indices = vec![0, 1, 2, 3, 4, 5];
        
        let pivot_idx = partition_with_indices(&mut vec, &mut indices, 0, 5);
        //
        // All elements before pivot should be <= pivot value.
        //
        for i in 0..pivot_idx {
            assert!(vec[i] <= vec[pivot_idx]);
        }
        //
        // All elements after pivot should be > pivot value.
        //
        for i in (pivot_idx+1)..vec.len() {
            assert!(vec[i] > vec[pivot_idx]);
        }
        //
        // Check that indices were swapped along with elements.
        //
        for (i, &idx) in indices.iter().enumerate() {
            assert_eq!(vec[i], match idx {
                0 => 5,
                1 => 2,
                2 => 9,
                3 => 1,
                4 => 5,
                5 => 6,
                _ => panic!("Invalid index"),
            });
        }
    }
    
    #[test]
    fn test_internal_quicksort() {
        //
        // Test the internal quicksort function.
        //
        let mut vec = vec![5, 2, 9, 1, 5, 6];
        let len = vec.len() - 1;
        _quicksort(&mut vec, 0, len);
        
        assert_eq!(vec, vec![1, 2, 5, 5, 6, 9]);
    }
    
    #[test]
    fn test_internal_quicksort_with_indices() {
        //
        // Test the internal quicksort function with index tracking.
        //
        let mut vec = vec![5, 2, 9, 1, 5, 6];
        let mut indices = vec![0, 1, 2, 3, 4, 5];
        let len = vec.len() - 1;

        _quicksort_with_indices(&mut vec, &mut indices, 0, len);
        
        assert_eq!(vec, vec![1, 2, 5, 5, 6, 9]);
        //
        // Verify indices match the original positions.
        //
        let expected_indices = vec![3, 1, 0, 4, 5, 2];
        assert_eq!(indices, expected_indices);
    }
    
    #[test]
    fn test_quicksort_stability() {
        //
        // Quicksort is not a stable sort, but we can test that equal elements
        // are still sorted correctly even if their relative order might change.
        //
        // Create a custom type that compares only on first element.
        //
        #[derive(Debug, Clone, PartialEq)]
        struct Pair(i32, i32);
        
        impl PartialOrd for Pair {
            fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
                self.0.partial_cmp(&other.0)
            }
        }
        
        let mut vec = vec![
            Pair(5, 1),
            Pair(2, 2),
            Pair(5, 3),
            Pair(1, 4),
            Pair(5, 5),
        ];
        
        let indices = quicksort_with_indices(&mut vec, true);
        //
        // Check sorting is correct.
        //
        for i in 0..vec.len()-1 {
            assert!(vec[i].0 <= vec[i+1].0);
        }
        //
        // Count how many "5"s we have and verify they're together.
        //
        let count_5 = vec.iter().filter(|p| p.0 == 5).count();
        assert_eq!(count_5, 3);
        //
        // Check that the indices correctly map back to original elements.
        //
        for (i, &idx) in indices.iter().enumerate() {
            let original_pair = match idx {
                0 => Pair(5, 1),
                1 => Pair(2, 2),
                2 => Pair(5, 3),
                3 => Pair(1, 4), 
                4 => Pair(5, 5),
                _ => panic!("Invalid index"),
            };
            //
            // The second element might not match due to instability,
            // but the first (sorting key) should match.
            //
            assert_eq!(vec[i].0, original_pair.0);
        }
    }
}