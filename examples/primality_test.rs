use lib_modulo::Modulus64;

fn main() {
    assert!(primality_test(2));
    assert!(!primality_test(2 * 3));

    // 9-th Mersenne number
    assert!(primality_test((1 << 61) - 1));
    // 2^62 - 1 = (2^32 + 1)(2^32 - 1)
    assert!(!primality_test(u64::MAX));
}

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

    for witness in [2, 325, 9375, 28178, 450775, 9780504, 1795265022] {
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
                continue;
            }
        }

        return false;
    }

    true
}
