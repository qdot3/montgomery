use crate::Modulus64;

/// Performs deterministic Miller-Rabin primality test.
///
/// # Time complexity
///
/// *O*(log *x*)
pub const fn primality_test(x: u64) -> bool {
    if x < 64 {
        return (super::PRIME_LT_64 >> x) & 1 == 1;
    } else if (super::COPRIME_2_3_5 >> (x % 30)) & 1 == 0 || x % 7 == 0 {
        return false;
    }

    let (s, d) = {
        let x = x - 1;
        let s = x.trailing_zeros();
        (s - 1, x >> s)
    };

    let modulus = Modulus64::new(x);
    let one = modulus.residue(1).x;
    // (a - a) r = 0 (mod x), r % x != 0
    let neg_one = x - one;

    // from <https://miller-rabin.appspot.com/>
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

    let mut i = 0;
    'test: while i < witness.len() {
        let mut mint = modulus.residue(witness[i]);
        i += 1;

        if mint.is_zero() {
            continue;
        }

        mint = mint.pow(d);
        if mint.x == one || mint.x == neg_one {
            continue;
        }

        let mut s = s;
        while s > 0 {
            s -= 1;

            mint.x = mint.modulus.mul(mint.x, mint.x);
            if mint.x == neg_one {
                continue 'test;
            }
        }

        return false;
    }

    true
}

// #[test]
// fn f() {
//     use std::io::Write;

//     let mut f = std::fs::File::create("./src/small_prime_context_u16_raw.rs").unwrap();

//     let _ = f.write(b"[\n");
//     for n in (3..1 << 16).step_by(2) {
//         if primality_test(n) {
//             let modulus = Context64::new(n);

//             let _ = f.write(format!("({}, {}, {}),\n", modulus.n, modulus.inv_n, modulus.r2_mod_n).as_bytes());
//         }
//     }
//     let _ = f.write(b"]");
// }

#[cfg(test)]
mod tests {
    use rand::{rng, Rng};

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
