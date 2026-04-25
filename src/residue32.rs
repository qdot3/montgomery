use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Factory of [`Residue32`].
///
/// See documentation of [`Residue32`] for details.
#[allow(clippy::derived_hash_with_manual_eq)]
#[derive(Debug, Clone, Hash, Eq)]
pub struct Modulus32 {
    // n inv_n = 1 (mod 2^64)
    n: u64,
    inv_n: u64,
    // 2^128 (mod n) * inv_n
    init: u64,
    // ceil(2^64 / n)
    recip: u64,
}

impl Modulus32 {
    /// Maximum available modulus.
    pub const MAX: u32 = 2_654_435_769;

    /// Creates new context for modular arithmetics.
    ///
    /// # Panics
    ///
    /// - modulus `n` should be an odd integer.
    /// - modulus `n` should be no more than `2_654_435_769`,
    ///   which is the floor of `2^32 / GOLDEN_RATIO`.
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus32;
    ///
    /// // odd integer less than or equal to 2_654_435_769 is allowed.
    /// let modulus = Modulus32::new(Modulus32::MAX);
    /// let modulus = Modulus32::new(3);
    ///
    /// // modulus should be an odd integer!
    /// assert!(std::panic::catch_unwind(|| { Modulus32::new(2); }).is_err())
    /// ```
    #[inline]
    pub const fn new(n: u32) -> Self {
        assert!(
            n & 1 == 1,
            "invalid modulus: modulus should be an odd integer."
        );
        assert!(
            n <= Self::MAX,
            "invalid modulus: modulus should be no more than 2_654_435_769."
        );

        let n = n as u64;

        let inv_n = {
            // 1 * 1 = 3 * 3 = 1 (mod 4)
            let mut inv_n = n & 3;
            // n inv_n = 1 (mod 2^k) => (n inv_n - 1)^2 = 0 (mod 2^{2k})
            // => n inv_n (2 - n inv_n) = 1 (mod 2^{2k})
            let mut i = u64::BITS.ilog2() - 1;
            while i > 0 {
                i -= 1;
                inv_n = inv_n.wrapping_mul(2_u64.wrapping_sub(n.wrapping_mul(inv_n)));
            }
            debug_assert!(n.wrapping_mul(inv_n) == 1);

            inv_n
        };

        let (div, rem) = {
            let denom = n.wrapping_neg();
            (denom / n, denom % n)
        };
        // 2^128 (mod n): magic number for converting integer to Plantard representation.
        let init = rem * rem % n;
        // ceil(2^64 / n): magic number for fast remainder algorithm
        let recip = div.wrapping_add(if rem > 0 { 2 } else { 1 });

        Self {
            n,
            inv_n,
            init: init.wrapping_mul(inv_n),
            recip,
        }
    }

    /// Performs Plantard multiplication, i.e. `x, y -> x y / -2^64 (mod n)`.
    ///
    /// If `x y < self.n`, then returned value is less than `self.n`.
    #[inline(always)]
    const fn mul(&self, x: u64, y: u64) -> u64 {
        // Plantard reduction: <https://thomas-plantard.github.io/pdf/Plantard21.pdf>
        let z = self.inv_n.wrapping_mul(x).wrapping_mul(y) >> 32;
        let z = ((z as u32).wrapping_add(1) as u64 * self.n) >> 32;
        debug_assert!(z < self.n, "this is a bug in lib-modulo");
        z
    }

    /// Calculates the residue of `x` modulo `self`.
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus32;
    ///
    /// let modulus = Modulus32::new(5);
    /// assert_eq!(modulus.residue(8).get(), 3)
    /// ```
    #[inline(always)]
    pub const fn residue(&self, x: u32) -> Residue32<'_> {
        // fast remainder algorithm
        // See <https://onlinelibrary.wiley.com/doi/10.1002/spe.2689> for details
        let x = {
            let lo = self.recip.wrapping_mul(x as u64);
            ((lo as u128 * self.n as u128) >> 64) as u64
        };

        let x = {
            // multiplication by a constant
            let x = self.init.wrapping_mul(x) >> 32;
            ((x as u32).wrapping_add(1) as u64 * self.n) >> 32
        };

        Residue32 { x, modulus: self }
    }

    /// Checks whether `x` is divisible by `self`.
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus32;
    ///
    /// let modulus = Modulus32::new(9);
    /// assert!(modulus.can_divide(18));
    /// assert!(!modulus.can_divide(19));
    /// ```
    #[inline(always)]
    pub const fn can_divide(&self, x: u32) -> bool {
        self.recip.wrapping_mul(x as u64) <= self.recip.wrapping_sub(1)
    }
}

impl PartialEq for Modulus32 {
    fn eq(&self, other: &Self) -> bool {
        // other fields depend on `n`
        self.n == other.n
    }
}

/// A residue with an odd modulus not exceeding `2_654_435_769`.
///
/// # Fast modular multiplication
///
/// [`Residue32`] provides fast modular multiplication using [Plantard multiplication].
/// This method eliminates one multiplication when one of the operands is reused multiple times.
/// As a result, [`Residue32::pow`] and other operations are typically
/// faster than implementations based on [Montgomery multiplication].
///
/// [Plantard multiplication]: https://thomas-plantard.github.io/pdf/Plantard21.pdf
/// [Montgomery multiplication]: https://doi.org/10.1090/s0025-5718-1985-0777282-x
///
/// # Usage
///
/// ```
/// use lib_modulo::Modulus32;
///
/// // set modulus
/// let modulus = Modulus32::new(3);
///
/// // performs modular arithmetics
/// let one = modulus.residue(1);
/// let two = modulus.residue(2);
/// let five = modulus.residue(5);
/// assert_eq!(two * five, one)
/// ```
///
/// Two residues with different modulus can interact, but the result will be meaningless.
/// It is highly recommended to use a block to ensure that [`Modulus32`], therefore [`Residue32`]s, are dropped.
#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
pub struct Residue32<'a> {
    // compare modulus first
    modulus: &'a Modulus32,
    x: u64,
}

impl<'a> Residue32<'a> {
    /// Extract the internal representation of `self`.
    ///
    /// ```
    /// use lib_modulo::{Modulus32, Raw32};
    ///
    /// let modulus = Modulus32::new(1001);
    /// // save memory
    /// let residues: Vec<Raw32> = (1..=1000).map(|x| modulus.residue(x).into_raw()).collect();
    /// ```
    #[inline(always)]
    pub const fn into_raw(self) -> Raw32 {
        Raw32 { x: self.x }
    }

    /// Checks whether `self` is `0`.
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus32;
    ///
    /// let modulus = Modulus32::new(5);
    /// assert!(modulus.residue(10).is_zero())
    /// ```
    #[inline(always)]
    pub const fn is_zero(self) -> bool {
        self.x == 0
    }

    /// Returns the residue.
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus32;
    ///
    /// let modulus = Modulus32::new(7);
    /// assert_eq!(modulus.residue(10).get(), 3)
    /// ```
    #[inline(always)]
    pub const fn get(self) -> u64 {
        self.modulus.mul(self.x, 1)
    }

    /// Returns the modulus.
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus32;
    ///
    /// let modulus = Modulus32::new(11);
    /// assert_eq!(modulus.residue(2).modulus(), 11);
    /// ```
    #[inline(always)]
    pub const fn modulus(&self) -> u64 {
        self.modulus.n
    }

    /// Raises `self` to the power of `exp`, using exponentiation by squaring.
    ///
    /// # Time complexity
    ///
    /// *Θ*(log `exp`)
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus32;
    ///
    /// let modulus = Modulus32::new(1001);
    /// let residue = modulus.residue(2);
    /// for exp in 0..64 {
    ///     assert_eq!(residue.pow(exp).get(), (1 << exp) % 1001)
    /// }
    /// ```
    #[inline(always)]
    pub const fn pow(self, mut exp: u32) -> Self {
        let Self { mut x, modulus } = self;
        // If `n = 1`, then `init = 0`. Otherwise, `n > 1`.
        let mut prod = modulus.residue(1).x;

        while exp > 1 {
            if exp & 1 == 1 {
                // インライン展開されると,掛け算を１回節約できる。
                prod = modulus.mul(prod, x)
            }

            exp >>= 1;
            x = modulus.mul(x, x); // skip last useless one
        }
        if exp != 0 {
            prod = modulus.mul(prod, x);
        }

        Self { x: prod, modulus }
    }

    /// Calculates the modular inverse of `self`, using extended binary GCD algorithm.
    ///
    /// Modular inverse can be defined if and only if `self` and the modulus is coprime.
    ///
    /// - `Ok(x)` : `x` is the modular inverse.
    /// - `Err(x)`: `x` is the GCD of `self` and the `modulus`,
    ///   where `gcd(0, a)` is defined to be `a`.
    ///
    /// # Time complexity
    ///
    /// *O*(log `self`)
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus32;
    ///
    /// let modulus = Modulus32::new(3 * 5);
    ///
    /// let residue = modulus.residue(2);
    /// assert!(residue.inv().is_ok_and(|inv| (inv * residue).get() == 1));
    ///
    /// let residue = modulus.residue(6);
    /// assert!(residue.inv().is_err_and(|gcd| gcd == 3));
    /// ```
    pub const fn inv(self) -> Result<Self, u64> {
        // invariant: [a] x = a, [a] y = b (mod n), where [a] is initial value.
        let mut a = self.get();
        let mut b = self.modulus();
        let Self { modulus, .. } = self;
        let mut x = modulus.residue(1).x;
        let mut y = 0;
        let frac_1_2 = modulus.residue((modulus.n as u32).div_ceil(2));

        while a > 0 {
            x = modulus.mul(x, frac_1_2.pow(a.trailing_zeros()).x);
            a >>= a.trailing_zeros();

            if a < b {
                (a, b) = (b, a);
                (x, y) = (y, x);
            }
            a -= b;
            let (z, b) = x.overflowing_sub(y);
            x = if b { z.wrapping_add(modulus.n) } else { z };
        }

        // b = gcd([a], n)
        if b == 1 {
            Ok(Self { x: y, modulus })
        } else {
            Err(b)
        }
    }
}

impl<'a> Add for Residue32<'a> {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        let (x, b) = self.x.overflowing_add(rhs.x);
        self.x = if b || x >= self.modulus() {
            x.wrapping_sub(self.modulus())
        } else {
            x
        };

        self
    }
}

impl<'a> AddAssign for Residue32<'a> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl<'a> Sub for Residue32<'a> {
    type Output = Self;

    fn sub(mut self, rhs: Self) -> Self::Output {
        let (x, b) = self.x.overflowing_sub(rhs.x);
        self.x = if b { x.wrapping_add(self.modulus()) } else { x };

        self
    }
}

impl<'a> SubAssign for Residue32<'a> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}

impl<'a> Mul for Residue32<'a> {
    type Output = Self;

    fn mul(mut self, rhs: Self) -> Self::Output {
        self.x = self.modulus.mul(self.x, rhs.x);
        self
    }
}

impl<'a> MulAssign for Residue32<'a> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs
    }
}

impl<'a> Neg for Residue32<'a> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.x = if self.x == 0 {
            0
        } else {
            self.modulus() - self.x
        };

        self
    }
}

/// An internal representation of [`Residue32`] without an associated [`Modulus32`].
///
/// Conceptually, [`Residue32`] = [`Raw32`] + [`Modulus32`].
/// [`Raw32`] stores the value part alone, without holding a reference to its modulus.
///
/// This separation is useful for reducing the size of collections of [`Residue32`]
/// and for avoiding self-referential structures when a type needs to contain both
/// a residue and its modulus.
#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
pub struct Raw32 {
    x: u64,
}

impl Raw32 {
    /// Attaches a modulus and returns a [`Residue32`].
    ///
    /// # Caution
    ///
    /// This does not perform validation or reduction.
    /// The caller must ensure the modulus is correct for this value.
    #[inline(always)]
    pub const fn into_residue<'a>(self, modulus: &'a Modulus32) -> Residue32<'a> {
        Residue32 { modulus, x: self.x }
    }
}

impl<'a> From<Residue32<'a>> for Raw32 {
    #[inline(always)]
    fn from(residue: Residue32<'a>) -> Self {
        Self { x: residue.x }
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    use proptest::prelude::*;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 15))]
        #[test]
        fn mul(n in (0..=Modulus32::MAX).prop_map(|n| n | 1), x: u32) {
            let modulus = Modulus32::new(n);

            let res = modulus.residue(x);
            assert_eq!(res.get() as u32, x % n)
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 15))]
        #[test]
        fn pow(n in (0..=Modulus32::MAX as u64).prop_map(|n| n | 1), x in 0u64..1 << 32) {
            let modulus = Modulus32::new(n as u32);

            let res = modulus.residue(x as u32);
            let mut naive = 1;
            for i in 0..100 {
                assert_eq!(res.pow(i).get(), naive, "exp = {i}");
                naive = naive * x % n
            }
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 15))]
        #[test]
        fn divisible(n in (0..=Modulus32::MAX).prop_map(|n| n | 1), x: u32) {
            let modulus = Modulus32::new(n);

            assert_eq!(modulus.can_divide(x), x % n == 0);
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 15))]
        #[test]
        fn divisible_by_1(x: u32) {
            assert!(Modulus32::new(1).can_divide(x))
        }
    }

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

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1 << 15))]
        #[test]
        fn inv(n in (0..=Modulus32::MAX).prop_map(|n| n | 1), x: u32) {
            let modulus = Modulus32::new(n);
            let res = modulus.residue(x);

            match res.inv() {
                Ok(inv) => assert_eq!((inv * res).get(), 1),
                Err(gcd) => {
                    assert!(res.get() % gcd == 0);
                    assert!(res.modulus() % gcd == 0);
                    assert_eq!(binary_gcd(res.get() / gcd, res.modulus() / gcd), 1);
                }
            }
        }
    }
}
