use super::Context;

/// `x` is prime iff `(test >> x) & 1 == 1`
const SMALL_PRIME: u64 = {
    let mut test = 0;
    let mut x = 64;
    while x > 0 {
        x -= 1;

        test = (test << 1) | if primality_test_naive(x) { 1 } else { 0 };
    }
    test
};

// from <https://miller-rabin.appspot.com/>
/// x < 350_269_456_337
static SET3: [u64; 3] = [
    4230279247111683200,
    14694767155120705706,
    16641139526367750375,
];
/// x < 2^64
static SET7: [u64; 7] = [2, 325, 9375, 28178, 450775, 9780504, 1795265022];

/// Performs Miller-Rabin primality test.
///
/// # Time complexity
///
/// *O*(log *x*)
pub const fn primality_test(x: u64) -> bool {
    /// remove multiples of 2, 3 or 5
    const MAY_BE_PRIME_30: u32 = {
        let mut table = 0;

        let mut n = 0;
        while n < 30 {
            table |= if n % 2 == 0 || n % 3 == 0 || n % 5 == 0 {
                0 // composite
            } else {
                1 // may be prime
            } << (n % 32);
            n += 1;
        }

        table
    };

    if x < 64 {
        return (SMALL_PRIME >> x) & 1 == 1;
    } else if (MAY_BE_PRIME_30 >> x % 30) & 1 == 0 || x % 7 == 0 {
        return false;
    }

    let (s, d) = {
        let x = x - 1;
        let s = x.trailing_zeros();
        (s - 1, x >> s)
    };

    let ctx = Context::<u64>::new(x);
    let one = ctx.modulo(1).value;
    debug_assert!(one != 0, "gcd(r, x) = 1, x > 1 => r % x != 0");
    // (a - a) r = 0 (mod x)
    let neg_one = x - one;

    let witness = if x < 350_269_456_337 {
        SET3.as_slice()
    } else {
        SET7.as_slice()
    };

    let mut i = 0;
    'test: while i < witness.len() {
        let mut mint = ctx.modulo(witness[i]);
        i += 1;

        if mint.value == 0 {
            continue;
        }

        mint = mint.pow(d);
        if mint.value == one || mint.value == neg_one {
            continue;
        }

        let mut s = s;
        while s > 0 {
            s -= 1;

            mint.value = mint.ctx.mul(mint.value, mint.value);
            if mint.value == neg_one {
                continue 'test;
            }
        }

        return false;
    }

    true
}

const fn primality_test_naive(x: u64) -> bool {
    if x < 2 {
        return false;
    }

    let mut d = 1;
    while d < x.isqrt() {
        d += 1;

        if x % d == 0 {
            return false;
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use rand::{rng, Rng};

    use super::*;

    #[test]
    fn small() {
        for x in 0..500_000 {
            assert_eq!(primality_test(x), primality_test_naive(x), "{x}")
        }
    }

    #[test]
    fn intermediate() {
        let mut rng = rng();

        for x in std::iter::repeat_with(|| rng.random_range(1 << 30..1 << 40)).take(100) {
            assert_eq!(primality_test(x), primality_test_naive(x), "{x}")
        }
    }
}
