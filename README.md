# Fast Modular Arithmetic

[![crate](https://img.shields.io/crates/v/lib_modulo.svg)](https://crates.io/crates/lib_modulo)
[![documentation](https://docs.rs/lib_modulo/badge.svg)](https://docs.rs/lib_modulo)

High-performance word-size modular arithmetic using Montgomery and Plantard multiplication.

## Overview

Naive modular multiplication typically requires widening the operands, followed by an expensive division.
This crate avoids division by using:

- Montgomery multiplication
- Plantard multiplication

These techniques significantly improve performance, especially when the modulus is determined at *runtime*.

## Features

- 🚀 Fast modular multiplication without division
- ⚡ Optimized 32-bit and 64-bit implementations
- 🔒 Supports any odd modulus

## Example

```rust
use lib_modulo::Modulus32;

let modulus = Modulus32::new((1 << 10) - 1);

// residue of 1025 modulo 1023
let two = modulus.residue(1025);
assert_eq!(two.get(), 2);

assert_eq!(two.pow(10).get(), 1);
```

See more [examples].

[examples]: https://github.com/qdot3/montgomery/tree/master/examples

## Further reading

Fast modular multiplication algorithm

- [Montgomery multiplication](https://doi.org/10.1090/s0025-5718-1985-0777282-x)
- [Plantard multiplication](https://thomas-plantard.github.io/pdf/Plantard21.pdf)

Fast remainder algorithm

- <https://onlinelibrary.wiley.com/doi/10.1002/spe.2689>
