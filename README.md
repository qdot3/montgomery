# lib-modulo

This crate provides an implementation of Montgomery multiplication along with additional utilities for modular arithmetic.

[Montgomery multiplication](https://doi.org/10.1090/s0025-5718-1985-0777282-x) is an efficient algorithm for modular multiplication, especially when the modulus is determined at runtime.

## Installation

Add this to your `Cargo.toml`.

```toml
lib-modulo = { git = "https://github.com/qdot3/montgomery" }
```
