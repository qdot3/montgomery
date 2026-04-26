/// A modulus in `[2, 2^32)`, including even values.
///
/// # Fast modular multiplication
///
/// Provides fast modular multiplication using [Barrett multiplication].
/// This works for any modulus (including even ones) and places no restrictions on the operands,
///  but is generally slower than [`Residue32`](crate::Residue32).
///
/// Unlike Montgomery or Plantard methods, this operates directly on standard
/// integer representations (i.e., no transformation is required).
///
/// [Barrett multiplication]: https://doi.org/10.1007/3-540-47721-7_24
///
/// # Usage
///
/// ```
/// use lib_modulo::Modulus32Any;
///
/// let modulus = Modulus32Any::new(14);
/// assert_eq!(modulus.mul(3, 5), 1)
/// ```
#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
pub struct Modulus32Any {
    // Lemire's remainder algorithm (N = 64, L = 32)
    n: u64,
    // magic number for Barrett multiplication and Lemire's remainder algorithm
    // ceil(2^64 / n)
    magic: u64,
}

impl Modulus32Any {
    /// Creates new instance with the given modulus.
    ///
    /// # Panics
    ///
    /// Modulus `n` should be greater than 1.
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus32Any;
    ///
    /// // even number is available
    /// let modulus = Modulus32Any::new(2);
    /// // odd number is also available
    /// let modulus = Modulus32Any::new(3);
    /// // division by 0 is undefined
    /// assert!(std::panic::catch_unwind(|| Modulus32Any::new(0)).is_err());
    /// // meaningless division by 1 is NOT available for performance
    /// assert!(std::panic::catch_unwind(|| Modulus32Any::new(1)).is_err());
    /// ```
    #[inline(always)]
    pub const fn new(n: u32) -> Self {
        assert!(n > 1);

        let n = n as u64;
        let magic = (u64::MAX / n).wrapping_add(1);

        Self { n, magic }
    }

    /// Returns the modulus.
    #[inline(always)]
    pub const fn modulus(&self) -> u32 {
        self.n as u32
    }

    /// Calculates residue of `x` modulo `self`.
    #[inline(always)]
    pub const fn residue32(&self, x: u32) -> u32 {
        let lo = self.magic.wrapping_mul(x as u64);
        ((lo as u128 * self.n as u128) >> 64) as u32
    }

    /// Calculates residue of `x` modulo `self`.
    #[inline(always)]
    pub const fn residue64(&self, x: u64) -> u64 {
        let quot = ((x as u128 * self.magic as u128) >> 64) as u64;
        let (rem, b) = x.overflowing_sub(quot * self.n);
        if b {
            rem.wrapping_add(self.n)
        } else {
            rem
        }
    }

    /// Checks if `x` is divisible by `self`.
    #[inline(always)]
    pub const fn can_divide(&self, x: u32) -> bool {
        // since `self.n` is not 1, `self.magic` never overflow
        self.magic.wrapping_mul(x as u64) < self.magic
    }

    /// Performs `*` operation modulo `self`.
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus32Any;
    ///
    /// // even number is available
    /// let modulus = Modulus32Any::new(1 << 8);
    /// assert_eq!(modulus.mul(u32::MAX, u32::MAX), 1);
    /// ```
    #[inline(always)]
    pub const fn mul(&self, a: u32, b: u32) -> u32 {
        self.residue64(a as u64 * b as u64) as u32
    }

    /// Performs `a * b + c` modulo `self`.
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus32Any;
    ///
    /// let modulus = Modulus32Any::new(2357);
    /// assert_eq!(
    ///     modulus.carrying_mul(123, 456, 789),
    ///     (123 * 456 + 789) % 2357
    /// );
    /// ```
    #[inline(always)]
    pub const fn carrying_mul(&self, a: u32, b: u32, c: u32) -> u32 {
        self.residue64(a as u64 * b as u64 + c as u64) as u32
    }

    /// Performs `a * b + c + d` modulo `self`.
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus32Any;
    ///
    /// // even number is available
    /// let modulus = Modulus32Any::new(123_456);
    /// assert_eq!(
    ///     modulus.carrying_mul_add(u32::MAX, u32::MAX, u32::MAX, u32::MAX),
    ///     (u64::MAX % 123_456) as u32
    /// );
    /// ```
    #[inline(always)]
    pub const fn carrying_mul_add(&self, a: u32, b: u32, c: u32, d: u32) -> u32 {
        self.residue64(a as u64 * b as u64 + c as u64 + d as u64) as u32
    }

    /// Raises `x` to the power of `exp`, using exponentiation by squaring.
    ///
    /// # Time Complexity
    ///
    /// *O*(log `x`)
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus32Any;
    ///
    /// let modulus = Modulus32Any::new(123_456);
    ///
    /// assert_eq!(modulus.pow(123_456 * 100 + 1, 1000), 1)
    /// ```
    #[inline(always)]
    pub const fn pow(&self, mut x: u32, mut exp: u32) -> u32 {
        let mut res = 1;
        while exp > 0 {
            if exp & 1 == 1 {
                res = self.mul(res, x);
            }
            exp >>= 1;
            x = self.mul(x, x);
        }
        res
    }

    /// Calculates the modular inverse of `x`, using extended gcd algorithm.
    ///
    /// Modular inverse can be defined if and only if `x` and the modulus is coprime.
    ///
    /// - `Ok(_)` : the modular inverse.
    /// - `Err(_)`: the GCD of `x` and the `modulus`, where `gcd(0, a)` is defined to be `a`.
    ///
    /// # Time complexity
    ///
    /// *O*(log `x`)
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus32Any;
    ///
    /// let modulus = Modulus32Any::new(3 * 5);
    /// for (a, inv_a) in [(1, 1), (2, 8), (4, 4), (7, 13), (11, 11), (14, 14)] {
    ///     assert_eq!(modulus.mul(a, inv_a), 1);
    ///     assert_eq!(modulus.inv(a), Ok(inv_a));
    /// }
    /// for a in [3, 6, 9, 12] {
    ///     assert_eq!(modulus.inv(a), Err(3));
    /// }
    /// for a in [5, 10] {
    ///     assert_eq!(modulus.inv(a), Err(5));
    /// }
    /// // gcd(0, a) is defined t be `a`
    /// assert_eq!(modulus.inv(0), Err(15));
    /// assert_eq!(modulus.inv(15 * 99), Err(15));
    /// ```
    pub const fn inv(&self, x: u32) -> Result<u32, u32> {
        // invariant: a x0 = x, b x0 = y (mod [y])
        let mut x = self.residue32(x) as i64;
        let mut y = self.n as i64;
        let [mut a, mut b] = [1, 0];

        while x > 0 {
            let (div, rem) = (y / x, y % x);

            (x, y) = (rem, x);
            (a, b) = (b - div * a, a);
        }

        // y = gcd(x0, y0) > 0
        if y != 1 {
            return Err(y as u32);
        }
        if b.is_negative() {
            b += self.n as i64;
        }

        Ok(b as u32)
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;
    use rand::{random_iter, rng};

    use super::Modulus32Any;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 15))]
        #[test]
        fn mul(n in 2..=u32::MAX, a: u32, b: u32) {
            let modulus = Modulus32Any::new(n);
            assert_eq!(
                modulus.mul(a, b),
                (a as u64 * b as u64 % n as u64) as u32,
                "{:?}", modulus
            );
        }
    }

    #[test]
    fn mul_small() {
        let mut rng = rng();
        for n in 2..1 << 8 {
            let modulus = Modulus32Any::new(n);
            for _ in 0..1 << 12 {
                let a = rng.random();
                let b = rng.random();
                assert_eq!(modulus.mul(a, b) as u64, (a as u64 * b as u64 % n as u64),)
            }
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 15))]
        #[test]
        fn residue32(n in 2..=u32::MAX, a: u32) {
            let modulus = Modulus32Any::new(n);
            assert_eq!(modulus.residue32(a), a % n);
        }
    }

    #[test]
    fn residue32_small() {
        for n in 2..1 << 8 {
            let modulus = Modulus32Any::new(n);
            for a in random_iter().take(1 << 12) {
                assert_eq!(modulus.residue32(a), a % n)
            }
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 15))]
        #[test]
        fn residue64(n in 2..=u32::MAX, a: u64) {
            let modulus = Modulus32Any::new(n);
            assert_eq!(modulus.residue64(a), a % n as u64);
        }
    }

    #[test]
    fn residue64_small() {
        for n in 2..1 << 8 {
            let modulus = Modulus32Any::new(n);
            for a in random_iter().take(1 << 12) {
                assert_eq!(modulus.residue64(a), a % n as u64)
            }
        }
    }

    fn binary_gcd(mut a: u32, mut b: u32) -> u32 {
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

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 15))]
        #[test]
        fn inv(n in 2..=u32::MAX, a: u32) {
            let modulus = Modulus32Any::new(n);
            match modulus.inv(a) {
                Ok(inv) => assert_eq!(modulus.mul(a, inv), 1, "!"),
                Err(gcd) => assert_eq!(gcd, binary_gcd(a, n), "?")
            }
        }
    }
}
