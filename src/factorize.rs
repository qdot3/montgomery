use crate::{prime::primality_test, Context64};

// 12 bytes * 6541 ~ 75 KiB
static SMALL_ODD_PRIME_CONTEXT_16: [(u16, u64, u16); 6541] =
    include!("./small_prime_context_u16_raw.rs");

/// \[*WIP*\] Factorize integer.
///
/// # Panics
///
/// Panics if you are unlucky. *Do not use this function in any production code.*
///
/// # Example
///
/// ```
/// use lib_modulo::factorize::*;
///
/// let mut factor = Vec::new();
/// factorize(998_244_353 * 1_000_000_007, &mut factor);
///
/// println!("{:?}", factor)
/// ```
pub fn factorize(mut x: u64, factor: &mut Vec<u64>) {
    if x < 2 {
        return;
    }

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

        while ctx.modulo(x).is_zero() {
            x /= ctx.n;
            factor.push(ctx.n);
        }

        if x == 1 {
            return;
        }
    }

    // find large prime factors (up to 3) by Pollard's rho
    while x > 1 {
        if primality_test(x) {
            factor.push(x);
            return;
        }

        // x is a composite number
        let d = pollard_rho(x);
        if primality_test(d) {
            while x % d == 0 {
                x /= d;
                factor.push(d);
            }
        }
    }
}

fn pollard_rho(x: u64) -> u64 {
    // x = Π_i p_i^{k_i} => period T >= lcm {p_i - 1}_i >= max {p_i - 1}_i
    let ctx = Context64::new(x);
    let one = ctx.modulo(1);
    // LCG => Pollard's "o"
    // gcd(a-1,x) = 1 => z = y + b/(a-1), z = a z
    let a1 = ctx.modulo(3);
    let b1 = ctx.modulo(2);
    let mut y1 = b1; // z != 0

    let mut prev = y1;
    let mut prod = one;
    let mut cnt = 0_u64;
    while !y1.is_zero() {
        // 0 -> a * 0 + b = b
        y1 = a1 * y1 + b1;
        prod *= y1;
        cnt += 1;

        if cnt % 512 == 0 || y1.is_zero() {
            let g = binary_gcd(prod.get(), x);

            if g == 1 {
                prev = y1;
                prod = one;
                continue;
            }
            if g > 0 && primality_test(g) {
                return g;
            }

            let mut cnt = 0;
            while {
                prev = a1 * prev + b1;
                let g = binary_gcd(prev.get(), x);

                if g > 1 {
                    if primality_test(g) {
                        return g;
                    } else {
                        return pollard_rho(g);
                    }
                }

                cnt += 1;
                cnt != 0
            } {}
            prod = one;
        }
    }

    panic!("This is a bug in `factorize()` function. Please report following two numbers:\ntarget: {x}\nperiod: {cnt}")
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
    use rand::{random_iter, rng, seq::SliceRandom, Rng};

    use super::*;

    fn base(mut primes: Vec<u64>) {
        assert!(
            primes.iter().all(|p| primality_test(*p)),
            "only prime numbers are accepted"
        );

        let mut factor = Vec::new();
        factorize(
            primes
                .iter()
                .fold(1, |prod, p| prod.checked_mul(*p).unwrap()),
            &mut factor,
        );

        factor.sort_unstable();
        primes.sort_unstable();
        assert_eq!(factor, primes)
    }

    #[test]
    fn large1() {
        base(vec![998_244_353, 1_000_000_007]);
    }

    #[test]
    fn large2() {
        base(vec![982_451_629, 982_451_653]);
    }

    #[test]
    fn large_square() {
        for p in [998_244_353, 1_000_000_007] {
            base(vec![p; 2]);
        }
    }

    #[test]
    fn random_small() {
        for n in random_iter::<u32>().take(10_000) {
            let mut factor = Vec::new();

            factorize(n as u64, &mut factor);
            assert_eq!(n as u64, factor.iter().product())
        }
    }

    #[test]
    fn random_square() {
        let mut rng = rng();
        for n in std::iter::repeat_with(|| rng.random_range(1 << 20..1 << 30)).take(30) {
            let mut factor = Vec::new();

            factorize(n * n, &mut factor);
            assert_eq!(n * n, factor.iter().product())
        }
    }

    // very slow!
    #[test]
    fn prime_square() {
        for n in (0..1 << 32)
            .rev()
            .step_by(2)
            .filter(|n| primality_test(*n))
            .take(1)
        {
            let mut factor = Vec::new();

            factorize(n * n, &mut factor);
            assert_eq!(n * n, factor.iter().product())
        }
    }

    // fast since p is relatively small
    #[test]
    fn triple() {
        let mut p = Vec::from_iter(
            (0..1 << 21)
                .rev()
                .step_by(2)
                .filter(|n| primality_test(*n))
                .take(100),
        );
        p.shuffle(&mut rng());

        for p in p.windows(3) {
            let n = p[0] * p[1] * p[2];
            let mut factor = Vec::new();

            factorize(n, &mut factor);
            assert_eq!(n, factor.iter().product())
        }
    }
}
