use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Factory of [`Residue64`].
///
/// See documentation of [`Residue64`] for details.
#[derive(Debug, Clone, Eq, Hash)]
pub struct Modulus64 {
    // n inv_n = 1 (mod r = 2^32 or 2^64)
    pub(crate) n: u64,
    pub(crate) inv_n: u64,
    pub(crate) r2_mod_n: u64,
}

impl Modulus64 {
    /// Calculates some parameters for Montgomery multiplication.
    ///
    /// # Panics
    ///
    /// - modulus `n` should be an odd number.
    #[inline]
    pub const fn new(n: u64) -> Self {
        assert!(n & 1 == 1, "modulus should be an odd number");

        let inv_n = {
            const TABLE: u32 = {
                // | n     | 1 | 3  | 5  | 7 | 9 | 11 | 13 | 15 |
                // | inv_n | 1 | 11 | 13 | 7 | 9 | 3  | 5  | 15 | <- 4 bits * 8
                let inv_n = [1, 11, 13, 7, 9, 3, 5, 15];

                let mut table = 0;
                let mut i = 0;
                while i < 8 {
                    table |= inv_n[i] << (i * 4);
                    i += 1;
                }

                table
            };
            // n inv_n = 1 (mod 8)
            let mut inv_n = ((TABLE >> (n & 0b1110) * 2) & 0b1111) as u64;

            let mut d = const { u64::BITS.ilog2() - 2 };
            while d > 0 {
                inv_n = inv_n.wrapping_mul((2 as u64).wrapping_sub(n.wrapping_mul(inv_n)));
                d -= 1;
            }
            debug_assert!(n.wrapping_mul(inv_n) == 1);

            inv_n
        };
        let r2_mod_n = ((n as u128).wrapping_neg() % (n as u128)) as u64;

        Self { n, inv_n, r2_mod_n }
    }

    #[inline(always)]
    pub const fn residue(&self, x: u64) -> Residue64<'_> {
        // `x r2 < r n`
        let x = self.mul(x, self.r2_mod_n);

        Residue64 { x, modulus: &self }
    }

    /// Performs Montgomery multiplication.
    ///
    /// if `lhs rhs < n r`, then `result < n`
    #[inline(always)]
    pub(crate) const fn mul(&self, lhs: u64, rhs: u64) -> u64 {
        self.mul_add(lhs, rhs, 0)
    }

    /// Performs `lhs rhs + add`.
    ///
    /// If `lhs rhs + add < n r`, then the result is less than `n`.
    #[inline(always)]
    pub(crate) const fn mul_add(&self, lhs: u64, rhs: u64, add: u64) -> u64 {
        // FIXME: use `a.widening_mul(b)`
        let (x_hi, x_lo) = {
            let x = lhs as u128 * rhs as u128 + add as u128;
            ((x >> u64::BITS) as u64, x as u64)
        };
        // FIXME: use `mul_hi()`
        // y = x n nn = x (mod r) => yl = x_lo
        let y_hi = ((x_lo.wrapping_mul(self.inv_n) as u128 * self.n as u128) >> u64::BITS) as u64;
        // x - y = 0 (mod r), x - y = x (mod n) => z = x inv_r (mod n)
        let (z, b) = x_hi.overflowing_sub(y_hi);

        // x < n r, y < n r => |z| < n
        if b {
            z.wrapping_add(self.n)
        } else {
            z
        }
    }

    /// Checks whether `x` is multiple of `self`.
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus64;
    ///
    /// for n in (1..1 << 10).step_by(2) {
    ///     let modulus  = Modulus64::new(n);
    ///
    ///     (0..1 << 10).for_each(|k| assert!(modulus.can_divide(n * k)));
    /// }
    /// ```
    #[inline]
    pub const fn can_divide(&self, x: u64) -> bool {
        // x < n r
        let x = self.mul(x, 1);
        x == 0
    }
}

impl PartialEq for Modulus64 {
    fn eq(&self, other: &Self) -> bool {
        self.n == other.n
    }
}

/// Residue with odd modulus which is less than 2^64.
///
/// # Usage
///
/// ```
/// use lib_modulo::Modulus64;
///
/// // runtime-specified *odd* modulus
/// let modulus = 5;
///
/// let modulus = Modulus64::new(modulus); // slow
/// let n = modulus.residue(2) * modulus.residue(3); // fast
/// assert_eq!(n.get(), 1);
/// ```
///
/// # Caution
///
/// [`Modulo`] values created from different [`Modulus64`]s can technically interact,
/// but the results will be meaningless.
/// It is recommended to use a block to ensure that each [`Modulus64`] is dropped
/// before another one is introduced.
#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
pub struct Residue64<'a> {
    pub(crate) modulus: &'a Modulus64,
    // x r (mod n)
    pub(crate) x: u64,
}

impl<'a> Residue64<'a> {
    /// Returns value.
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus64;
    ///
    /// let n = 101;
    /// let modulus  = Modulus64::new(n);
    ///
    /// let n = modulus.residue(99);
    ///
    /// assert_eq!(n.get(), 99);
    /// assert_eq!(n.modulus(), 101);
    /// ```
    #[inline(always)]
    pub const fn get(&self) -> u64 {
        self.modulus.mul(self.x, 1)
    }

    /// Returns modulus.
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus64;
    ///
    /// let n = 101;
    ///     let modulus  = Modulus64::new(n);
    ///
    /// let n = modulus.residue(99);
    ///
    /// assert_eq!(n.get(), 99);
    /// assert_eq!(n.modulus(), 101);
    /// ```
    #[inline(always)]
    pub const fn modulus(&self) -> u64 {
        self.modulus.n
    }

    /// Returns `true` if `self` is `0`.
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus64;
    ///
    /// for n in (1..100_000).step_by(2) {
    ///     let modulus = Modulus64::new(n);
    ///     assert!(modulus.residue(0).is_zero());
    /// }
    /// ```
    #[inline(always)]
    pub const fn is_zero(self) -> bool {
        self.x == 0
    }

    /// Raises `self` to the power of `exp`, using exponentiation by squaring.
    ///
    /// # Time complexity
    ///
    /// *O*(log `exp`)
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus64;
    ///
    /// let n = 12_345;
    ///     let modulus  = Modulus64::new(n);
    ///
    /// let mut pow10 = 1;
    /// for i in 0..1_000 {
    ///     assert_eq!(modulus.residue(10).pow(i).get(), pow10);
    ///     pow10 = pow10 * 10 % n;
    /// }
    /// ```
    #[inline]
    pub const fn pow(mut self, mut exp: u64) -> Self {
        // r inv_r = 1 (mod n)
        let mut result = self.modulus.residue(1).x;

        while exp > 0 {
            if exp & 1 == 1 {
                // n < r
                result = self.modulus.mul(result, self.x)
            }

            exp >>= 1;
            // n < r
            self.x = self.modulus.mul(self.x, self.x)
        }
        self.x = result;

        self
    }

    /// Calculates the modular inverse of `self`, using extended binary GCD algorithm.
    ///
    /// Modular inverse can be defined if and only if `self` and the modulus is coprime.
    ///
    /// - `Ok(x)` : `x` is the modular inverse.
    /// - `Err(x)`: `x` is the GCD of `self` and the `modulus`,
    /// where `gcd(0, a) = gcd(a, 0)` is defined to be `a`.
    ///
    /// # Time complexity
    ///
    /// *O*(log `self`)
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus64;
    ///
    /// // 998_244_353 is a prime number, so modular inverse of n exists iff n != 0 (mod 998_244_353)
    /// let modulus = Modulus64::new(998_244_353);
    ///
    /// for n in 1..500_000 {
    ///     let n = modulus.residue(n);
    ///     assert!(n.try_inv().is_ok_and(|i| (i * n).get() == 1));
    /// }
    /// // 0 n = 0 != 1 for any integer n
    /// assert!(modulus.residue(0).try_inv().is_err());
    /// ```
    #[inline]
    pub const fn try_inv(self) -> Result<Self, u64> {
        let mut a = self.get();
        let Self { modulus, .. } = self;

        // performs extended binary gcd
        //
        // invariants: a = [a] x,  b = [a] y (mod n) where [a] is initial value
        let mut b = modulus.n;
        let mut x = modulus.residue(1).x; // 1 r mod n
        let mut y = 0; // 0 r mod n
        let frac_1_2 = modulus.residue(modulus.n.div_ceil(2));

        while a > 0 {
            x = modulus.mul(x, frac_1_2.pow(a.trailing_zeros() as u64).x);
            a >>= a.trailing_zeros();

            if a < b {
                (a, b) = (b, a);
                (x, y) = (y, x);
            }
            a -= b;
            let (diff, b) = x.overflowing_sub(y);
            x = if b {
                diff.wrapping_add(modulus.n)
            } else {
                diff
            };
        }

        // b = gcd([a], [b])
        if b == 1 {
            Ok(Self { x: y, modulus })
        } else {
            Err(b)
        }
    }
}

impl<'a> Add for Residue64<'a> {
    type Output = Self;

    #[inline(always)]
    fn add(mut self, rhs: Self) -> Self {
        let (sum, b) = self.x.overflowing_add(rhs.x);
        self.x = if b || sum >= self.modulus.n {
            sum.wrapping_sub(self.modulus.n)
        } else {
            sum
        };

        self
    }
}

impl<'a> AddAssign for Residue64<'a> {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl<'a> Sub for Residue64<'a> {
    type Output = Self;

    #[inline(always)]
    fn sub(mut self, rhs: Self) -> Self {
        let (diff, b) = self.x.overflowing_sub(rhs.x);
        self.x = if b {
            diff.wrapping_add(self.modulus.n)
        } else {
            diff
        };

        self
    }
}

impl<'a> SubAssign for Residue64<'a> {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}

impl<'a> Mul for Residue64<'a> {
    type Output = Self;

    #[inline(always)]
    fn mul(mut self, rhs: Self) -> Self {
        // n < r
        self.x = self.modulus.mul(self.x, rhs.x);

        self
    }
}

impl<'a> MulAssign for Residue64<'a> {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs
    }
}

impl<'a> Neg for Residue64<'a> {
    type Output = Self;

    #[inline(always)]
    fn neg(mut self) -> Self::Output {
        // (x - x) r = 0 (mod n)
        self.x = if self.x == 0 {
            self.x
        } else {
            self.modulus.n - self.x
        };

        self
    }
}
