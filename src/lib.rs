//! # Fast modular arithmetic
//!
//! Provides efficient types for modular arithmetic.
//!
//! Naive modular multiplication typically requires widening the operands
//! followed by an expensive division.
//! This crate avoids division by using [Montgomery] or [Plantard] multiplication,
//! enabling faster modular operations.
//!
//! [Plantard]: https://thomas-plantard.github.io/pdf/Plantard21.pdf
//! [Montgomery]: https://doi.org/10.1090/s0025-5718-1985-0777282-x
//!
//! # Guide
//!
//! Depending on the modulus, you can choose between two types:
//!
//! - [`Residue32`]: for odd moduli up to `2_654_435_769` (~ `2^31.3`)
//! - [`Residue64`]: for any odd modulus that fits in `u64`
//!
//! Since [`Residue32`] is typically faster than [`Residue64`], prefer using it whenever possible.
//!
//! # Example
//!
//! Below is an implementation of a rolling hash algorithm using [`Residue64`].
//! This allows the use of large prime moduli without overflow.
//!
//! ```
#![doc = include_str!("../examples/rolling_hash.rs")]
//! ```
#![warn(missing_docs)]
#![warn(missing_debug_implementations)]
pub mod factorize;
pub mod prime;

mod residue32;
pub use residue32::*;

mod residue64;
pub use residue64::*;
