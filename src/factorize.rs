use std::num::NonZero;

use crate::{prime::primality_test, Context64};

// 12 bytes * 6541 ~ 75 KiB
static SMALL_ODD_PRIME_CONTEXT_16: [(u16, u64, u16); 6541] =
    include!("./small_prime_context_u16_raw.rs");

/// Factorize integer and writes prime factors to `factor` in any order.
///
/// This function is probabilistic and may fail.
///
/// # Time complexity
///
/// O(`x`^0.25) expected
///
/// # Example
///
/// ```
/// use lib_modulo::factorize::*;
///
/// let mut factor = Vec::new();
/// // panics if factorization fails
/// assert!(factorize(998_244_353 * 1_000_000_007, &mut factor).is_ok());
///
/// factor.sort_unstable();
/// assert_eq!(factor, vec![998_244_353, 1_000_000_007])
/// ```
pub fn factorize(mut x: u64, factor: &mut Vec<u64>) -> Result<(), ()> {
    if x < 2 {
        return Ok(());
    }
    factor.reserve(64);

    // trial division by small primes less than 2^16
    {
        factor.extend(std::iter::repeat_n(2, x.trailing_zeros() as usize));
        x >>= x.trailing_zeros();
    }
    for &(n, inv_n, r2_mod_n) in SMALL_ODD_PRIME_CONTEXT_16.iter() {
        let ctx = Context64 {
            n: n as u64,
            inv_n,
            r2_mod_n: r2_mod_n as u64,
        };

        while ctx.can_divide(x) {
            x /= ctx.n;
            factor.push(ctx.n);
        }

        if x == 1 {
            return Ok(());
        }
    }

    // find large prime factors (up to 3) by Pollard's rho
    while x > 1 {
        if primality_test(x) {
            factor.push(x);
            return Ok(());
        }

        if let Some(d) = pollard_rho(x) {
            let d = d.get();
            while x % d == 0 {
                x /= d;
                factor.push(d);
            }
        } else {
            return Err(());
        }
    }

    Ok(())
}

/// Find prime factor of `x`.
///
/// This function is probabilistic and may fail.
///
/// # Time complexity
///
/// *O*(p^0.25) where p is a prime factor of `x`
fn pollard_rho(x: u64) -> Option<NonZero<u64>> {
    let ctx = Context64::new(x);
    let one = ctx.modulo(1);

    for c in 1..100 {
        // a = b (mod x) => f(a) = f(b) (mod x)
        let f = |x: u64| ctx.mul_add(x, x, c);

        let mut y0 = ctx.modulo(1);
        let mut y1 = y0;

        let mut prod = one;
        let mut step = 0;
        let mut memo = [[0, 0, one.value]; 1 << 5];

        'cycle_detection: while !prod.is_zero() {
            y0.value = f(y0.value);
            y1.value = f(f(y1.value));
            prod *= y1 - y0;
            step += 1;

            if step % (1 << 5) == 0 {
                memo[(step >> 5) % (1 << 5)] = [y0.value, y1.value, prod.value];
            }
            if step % (1 << 10) == 0 || prod.is_zero() {
                let g = binary_gcd(prod.value, x);

                if g == 1 {
                    continue 'cycle_detection;
                } else if primality_test(g) {
                    return NonZero::new(g);
                }

                for i in 0..memo.len() {
                    let g = binary_gcd(memo[i][2], x);

                    if g != 1 {
                        if primality_test(g) {
                            return NonZero::new(g);
                        }

                        y0.value = memo[i][0];
                        y1.value = memo[i][1];
                        for _ in 0..1 << 5 {
                            let g = binary_gcd((y0 - y1).value, x);

                            if g != 1 {
                                if primality_test(g) {
                                    return NonZero::new(g);
                                } else if g != x {
                                    // FIXME: `x` is composed of at most 3 primes, so return `x/g`
                                    return pollard_rho(g);
                                } else {
                                    break 'cycle_detection;
                                }
                            }

                            y0.value = f(y0.value);
                            y1.value = f(f(y1.value));
                        }
                    }
                }
            }
        }
    }

    None
}

#[inline(always)]
fn binary_gcd(mut a: u64, mut b: u64) -> u64 {
    if b == 0 {
        return a;
    }

    let shift = (a | b).trailing_zeros();
    b >>= b.trailing_zeros();

    while a != 0 {
        a >>= a.trailing_zeros();

        if a < b {
            (a, b) = (b, a)
        }
        a -= b
    }

    b << shift
}

#[cfg(test)]
mod tests {
    use rand::{rng, seq::SliceRandom, Rng};

    use super::*;

    #[test]
    fn random() {
        let mut rng = rng();
        for n in std::iter::repeat_with(|| rng.random_range(1 << 55..=u64::MAX)).take(10_000) {
            let mut factor = Vec::new();

            assert!(factorize(n, &mut factor).is_ok());
            assert_eq!(n, factor.iter().product())
        }
    }

    #[test]
    fn random_square() {
        let mut rng = rng();
        for n in std::iter::repeat_with(|| rng.random_range(1 << 20..1 << 32)).take(5000) {
            let mut factor = Vec::new();

            assert!(factorize(n * n, &mut factor).is_ok());
            assert_eq!(n * n, factor.iter().product())
        }
    }

    #[test]
    fn prime_square() {
        for n in (0..1 << 32)
            .rev()
            .step_by(2)
            .filter(|n| primality_test(*n))
            .take(500)
        {
            let mut factor = Vec::new();

            assert!(factorize(n * n, &mut factor).is_ok());
            assert_eq!(n * n, factor.iter().product())
        }
    }

    // fast since p is relatively small
    #[test]
    fn prime_cube() {
        let p = Vec::from_iter(
            (0..1 << 21)
                .rev()
                .step_by(2)
                .filter(|n| primality_test(*n))
                .take(500),
        );

        for p in p {
            let n = p.pow(3);
            let mut factor = Vec::new();

            assert!(factorize(n, &mut factor).is_ok());
            assert_eq!(n, factor.iter().product())
        }
    }

    #[test]
    fn prime_double() {
        let mut p: Vec<u64> = (0..1 << 32)
            .rev()
            .step_by(2)
            .filter(|n| primality_test(*n))
            .take(500)
            .collect();
        p.shuffle(&mut rng());

        for p in p.windows(2) {
            let n = p[0] * p[1];
            let mut factor = Vec::new();

            assert!(factorize(n, &mut factor).is_ok());
            assert_eq!(n, factor.iter().product())
        }
    }

    #[test]
    fn prime_triple() {
        let mut p: Vec<u64> = (0..1 << 21)
            .rev()
            .step_by(2)
            .filter(|n| primality_test(*n))
            .take(500)
            .collect();
        p.shuffle(&mut rng());

        for p in p.windows(3) {
            let n = p[0] * p[1] * p[2];
            let mut factor = Vec::new();

            assert!(factorize(n, &mut factor).is_ok());
            assert_eq!(n, factor.iter().product())
        }
    }
}
