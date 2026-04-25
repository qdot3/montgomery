use crate::{Modulus32, Modulus64};

/// Performs deterministic Miller-Rabin primality test.
///
/// # Time complexity
///
/// *O*(log *x*)
pub fn primality_test(x: u64) -> bool {
    if x < 64 {
        return (PRIME_LT_64 >> x) & 1 == 1;
    } else if (COPRIME_2_3_5 >> (x % 30)) & 1 == 0 || x % 7 == 0 {
        return false;
    }

    let (s, d) = {
        let x = x - 1;
        let s = x.trailing_zeros();
        (s - 1, x >> s)
    };

    macro_rules! primality_test_impl {
        ($modulus:expr, $witness:expr, $d:expr, $s:expr) => {
            let modulus = $modulus;
            let witness = $witness;
            let (d, s) = ($d, $s);

            let one = modulus.residue(1);
            let minus_one = -one;

            let mut i = 0;
            'test: while i < witness.len() {
                let mut mint = modulus.residue(witness[i]);
                i += 1;

                if mint.is_zero() {
                    continue;
                }

                mint = mint.pow(d);
                if mint == one || mint == minus_one {
                    continue;
                }

                let mut s = s;
                while s > 0 {
                    s -= 1;

                    mint = mint * mint;
                    if mint == minus_one {
                        continue 'test;
                    }
                }

                return false;
            }
        };
    }

    // witnesses from <https://miller-rabin.appspot.com/>
    if x <= Modulus32::MAX as u64 {
        primality_test_impl!(Modulus32::new(x as u32), [2, 7, 61], d as u32, s as u32);
    } else {
        let witness = if x < 350_269_456_337 {
            static SET3: [u64; 3] = [0x3AB4F88FF0CC7C80, 0xCBEE4CDF120C10AA, 0xE6F1343B0EDCA8E7];
            SET3.as_slice()
        } else if x < 7_999_252_175_582_851 {
            static SET5: [u64; 5] = [
                2,
                0x3C1C7396F6D,
                0x2142E2E3F22DE5C,
                0x297105B6B7B29DD,
                0x370EB221A5F176DD,
            ];
            SET5.as_slice()
        } else {
            static SET7: [u64; 7] = [2, 325, 9375, 28178, 450775, 9780504, 1795265022];
            SET7.as_slice()
        };

        primality_test_impl!(Modulus64::new(x), witness, d, s);
    }

    true
}

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

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

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

    #[test]
    fn small() {
        for x in 0..1 << 15 {
            assert_eq!(primality_test(x), primality_test_naive(x), "{x}")
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 15))]
        #[test]
        fn random1(x in 1u64 << 15..1 << 20) {
            assert_eq!(primality_test(x), primality_test_naive(x), "{x}")
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 15))]
        #[test]
        fn random2(x in 1u64 << 20..1 << 25) {
            assert_eq!(primality_test(x), primality_test_naive(x), "{x}")
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 13))]
        #[test]
        fn random3(x in 1u64 << 25..1 << 30) {
            assert_eq!(primality_test(x), primality_test_naive(x), "{x}")
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 10))]
        #[test]
        fn random4(x in 1u64 << 30..1 << 40) {
            assert_eq!(primality_test(x), primality_test_naive(x), "{x}")
        }
    }

    #[test]
    fn mersenne() {
        for d in 0..=64 {
            // Mersenne prime in `u64`
            let is_prime = [2, 3, 5, 7, 13, 17, 19, 31, 61].contains(&d);
            let x = (1u128 << d) - 1;

            assert_eq!(primality_test(x as u64), is_prime)
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 18))]
        #[test]
        fn composite(x: u32) {
            assert!(!primality_test(x as u64 * (x as u64 + 1)));
        }
    }
}
