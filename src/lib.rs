/*!
# Fast Modular Arithmetic

[![crate](https://img.shields.io/crates/v/lib_modulo.svg)](https://crates.io/crates/lib_modulo)
[![documentation](https://docs.rs/lib_modulo/badge.svg)](https://docs.rs/lib_modulo)

High-performance word-size modular arithmetic using Barrett, Montgomery or Plantard multiplication.

## Overview

Naive modular multiplication typically requires widening the operands,followed by an expensive division.
This crate avoids division by using:

- [Barrett multiplication](https://doi.org/10.1007/3-540-47721-7_24)
- [Montgomery multiplication](https://doi.org/10.1090/s0025-5718-1985-0777282-x)
- [Plantard multiplication](https://thomas-plantard.github.io/pdf/Plantard21.pdf)


These techniques significantly improve performance, especially when the modulus is determined at *runtime*.

## Selection guide

| Type           | Modulus             | Notes                 |
|----------------|---------------------|-----------------------|
| `Modulus32`    | odd, ≲ 2^31.3       | fastest               |
| `Modulus32Any` | in `[2, 2^32)`      | supports even moduli  |
| `Modulus64`    | odd, fits in `u64`  | supports large moduli |

## Advanced usage

`Residue` types hold a reference to their corresponding `Modulus` for convenience.
However, when storing many residues sharing the same modulus, this can be memory-intensive.

In such cases, [`Raw32`] or [`Raw64`] can be used instead.
The caller is responsible for associating them with the correct modulus.

# Example

Below is an implementation of a rolling hash algorithm using [`Residue64`].
This allows the use of large prime moduli without overflow.

```
*/
#![doc = include_str!("../examples/rolling_hash.rs")]
//! ```
#![warn(missing_docs, missing_debug_implementations)]
#![warn(clippy::all, clippy::pedantic, clippy::cargo)]
#![deny(rust_2018_idioms)]
#![forbid(unsafe_code)]
pub mod prime;

mod residue32;
pub use residue32::*;

mod residue32any;
pub use residue32any::Modulus32Any;

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
