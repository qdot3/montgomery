# lib-modulo

This crate provides fast modular arithmetics with runtime-specified odd modulus.

## Example

```rust
/// Miller-Rabin primality test
pub fn primality_test(x: u64) -> bool {
    if x == 2 || x == 3 || x == 5 || x == 7 {
        return true;
    } else if x % 2 == 0 || x % 3 == 0 || x % 5 == 0 || x % 7 == 0 {
        return false;
    }

    let (s, d) = {
        let x = x - 1;
        let s = x.trailing_zeros();
        (s - 1, x >> s)
    };

    let modulus = Modulus64::new(x);
    let one = modulus.residue(1);

    'test: for witness in [2, 325, 9375, 28178, 450775, 9780504, 1795265022] {
        let mut mint = modulus.residue(witness);

        if mint.is_zero() {
            continue;
        }

        mint = mint.pow(d);
        if mint == one || mint == -one {
            continue;
        }

        let mut s = s;
        while s > 0 {
            s -= 1;

            mint = mint * mint;
            if mint == -one {
                continue 'test;
            }
        }

        return false;
    }

    true
}
```

## Future Plan

These functions will be added:

- `sqrt()`
- `nth_root()`
- `log()`

## Further reading

- [Montgomery multiplication](https://doi.org/10.1090/s0025-5718-1985-0777282-x)
- [Plantard multiplication](https://thomas-plantard.github.io/pdf/Plantard21.pdf)
