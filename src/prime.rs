use super::Context;

/// `x` is prime iff `(test >> x) & 1 == 1`
const SMALL_TEST: u64 = {
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
    if x < 64 {
        return (SMALL_TEST >> x) & 1 == 1;
    } else if x % 2 == 0 || x % 3 == 0 || x % 5 == 0 || x % 7 == 0 {
        return false;
    }

    let (s, d) = {
        let x = x - 1;
        let s = x.trailing_zeros();
        (s - 1, x >> s)
    };

    let ctx = Context::<u64>::new(x);
    let one = ctx.modulo(1).value;
    debug_assert!(one != 0, "gcd(r, n) = 1, r > 1, n > 1 => r % n != 0");
    let neg_one = x - one;

    let mut i = 0;
    let witness = if x < 350_269_456_337 {
        SET3.as_slice()
    } else {
        SET7.as_slice()
    };

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
