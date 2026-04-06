use crate::prime::primality_test;

pub fn factorize(mut x: u64, factor: &mut Vec<u64>) {
    if x < 2 {
        return;
    }

    // divide by small prime numbers
    if x & 1 == 1 {
        factor.extend(std::iter::repeat_n(2, x.trailing_zeros() as usize));
        x >>= x.trailing_zeros();
    }
    for d in [3, 5, 7, 11] {
        while x % d == 0 {
            x /= d;
            factor.push(d);
        }
    }

    if primality_test(x) {
        return;
    }

    // Pollard's rho
    unimplemented!("")
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
        (a, b) = if a > b { (a - b, b) } else { (b - a, a) }
    }

    b << shift
}
