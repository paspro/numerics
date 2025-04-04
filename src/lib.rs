// -------------------------------------------------------------------------------------------------
//
//  This library implements various numerical algorithms which can be of use for all
//  kinds of scientific and engineering applications.
//
//  Copyright (c) 2025 by Dr. Panos Asproulis (p.asproulis@icloud.com).
//  All Rights Reserved.
//
// -------------------------------------------------------------------------------------------------

//!  This library implements various numerical algorithms which can be of use for all
//!  kinds of scientific and engineering applications.
//!
//!  Copyright (c) 2025 by Dr. Panos Asproulis (p.asproulis@icloud.com).
//!  All Rights Reserved.

pub mod interpolations;
pub mod root_finder;
pub mod search;
pub mod sorting;
pub mod determinants;
pub mod equations_solver;

pub use interpolations::*;
pub use root_finder::*;
pub use search::*;
pub use sorting::*;
pub use determinants::*;
pub use equations_solver::*;
