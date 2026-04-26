/*!
# Fast Modular Arithmetic

[![crate](https://img.shields.io/crates/v/lib_modulo.svg)](https://crates.io/crates/lib_modulo)
[![docs](https://docs.rs/lib_modulo/badge.svg)](https://docs.rs/lib_modulo)
![no_std](https://img.shields.io/badge/no__std-compatible-blue)
![unsafe](https://img.shields.io/badge/unsafe-forbidden-success)

High-performance word-size modular arithmetic using Barrett, Montgomery or Plantard multiplication.

# Overview

Naive modular multiplication typically requires widening the operands, followed by an expensive division.
This crate avoids division by using:

- [Barrett multiplication](https://doi.org/10.1007/3-540-47721-7_24)
- [Montgomery multiplication](https://doi.org/10.1090/s0025-5718-1985-0777282-x)
- [Plantard multiplication](https://thomas-plantard.github.io/pdf/Plantard21.pdf)

These techniques significantly improve performance, especially when the modulus is determined at *runtime*.

# Selection guide

| Type             | Modulus             | Notes                                              |
|------------------|---------------------|----------------------------------------------------|
| [`Modulus32`]    | odd, < 2^31.3...    | fastest                                            |
| [`Modulus32Any`] | in `[2, 2^32)`      | supports even moduli, zero-cost modulus switching  |
| [`Modulus64`]    | odd, fits in `u64`  | supports large moduli                              |

# Advanced usage

`Residue{N}` types hold a reference to their corresponding `Modulus{N}` for convenience.
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
#![forbid(unsafe_code)]
#![no_std]

mod residue32;
pub use residue32::*;

mod residue32any;
pub use residue32any::*;

mod residue64;
pub use residue64::*;

mod prime;
pub use prime::*;
