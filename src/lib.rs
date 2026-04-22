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
//! ## Basic usage
//!
//! Depending on the modulus, you can choose between two types:
//!
//! - [`Residue32`]: for odd moduli up to `2_654_435_769` (~ `2^31.3`)
//! - [`Residue64`]: for any odd modulus that fits in `u64`
//!
//! Since [`Residue32`] is typically faster than [`Residue64`], prefer using it whenever possible.
//!
//! ## Advanced usage
//!
//! `Residue` types hold a reference to their corresponding `Modulus` for convenience.
//! However, when storing many residues with the same modulus, this can be memory-intensive.
//! In such cases, [`Raw32`] or [`Raw64`] can be used instead.
//! The caller is responsible for associating them with the correct modulus.
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

/// small prime numbers less than 64
///
/// `(SELF >> x) & 1 == 1` iff `x` is prime
const PRIME_LT_64: u64 = {
    let mut test = u64::MAX << 2;
    let mut x = 2;
    while x < 64 {
        if (test >> x) & 1 == 1 {
            let mut y = x * x;
            while y < 64 {
                test &= !(1 << y);
                y += x;
            }
        }

        x += 1;
    }
    test
};

/// multiples of 2, 3 or 5.
///
/// `(SELF >> x % 30) & 1 == 1` iff `x` is coprime to 2, 3, and 5
const COPRIME_2_3_5: u32 = {
    let mut table = 0;

    let mut n = 0;
    while n < 30 {
        table |= if n % 2 == 0 || n % 3 == 0 || n % 5 == 0 {
            0 // composite
        } else {
            1 // may be prime
        } << n;
        n += 1;
    }

    table
};
