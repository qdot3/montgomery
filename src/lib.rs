use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

pub mod factorize;
pub mod prime;

pub type Context64 = Context<u64>;
pub type Context32 = Context<u32>;

pub type Modulo64<'a> = Modulo<'a, u64>;
pub type Modulo32<'a> = Modulo<'a, u32>;

/// Storage of parameters for Montgomery multiplication.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Context<U> {
    // n inv_n = 1 (mod r = 2^32 or 2^64)
    n: U,
    inv_n: U,
    r2_mod_n: U,
}

/// Modulo with a runtime-specified odd modulus.
///
/// # Usage
///
/// ```
/// use lib_modulo::Context64;
///
/// // runtime-specified *odd* modulus
/// let modulus = 5;
///
/// let ctx = Context64::new(modulus); // slow
/// let n = ctx.modulo(2) * ctx.modulo(3); // fast
/// assert_eq!(n.get(), 1);
/// ```
///
/// # Caution
///
/// [`Modulo`] values created from different [`Context`]s can technically interact,
/// but the results will be meaningless.
/// It is recommended to use a block to ensure that each [`Context`] is dropped
/// before another one is introduced.
#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
pub struct Modulo<'a, U> {
    // x r (mod n)
    value: U,
    ctx: &'a Context<U>,
}

macro_rules! montgomery_impl {
    ( $single:ty, $double:ty ) => {
        impl Context<$single> {
            /// Calculates some parameters for Montgomery multiplication.
            ///
            /// # Panics
            ///
            /// - modulus `n` should be an odd number.
            pub const fn new(n: $single) -> Self {
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
                    let mut inv_n = ((TABLE >> (n & 0b1110) * 2) & 0b1111) as $single;

                    let mut d = const { <$single>::BITS.ilog2() - 2 };
                    while d > 0 {
                        inv_n =
                            inv_n.wrapping_mul((2 as $single).wrapping_sub(n.wrapping_mul(inv_n)));
                        d -= 1;
                    }
                    debug_assert!(n.wrapping_mul(inv_n) == 1);

                    inv_n
                };
                let r2_mod_n = ((n as $double).wrapping_neg() % (n as $double)) as $single;

                Self { n, inv_n, r2_mod_n }
            }

            #[inline(always)]
            pub const fn modulo(&self, x: $single) -> Modulo<'_, $single> {
                // `x r2 < r n`
                let x = self.mul(x, self.r2_mod_n);

                Modulo {
                    value: x,
                    ctx: &self,
                }
            }

            /// Performs Montgomery multiplication.
            ///
            /// if `lhs rhs < n r`, then `result < n`
            #[inline(always)]
            const fn mul(&self, lhs: $single, rhs: $single) -> $single {
                // FIXME: use `a.widening_mul(b)`
                let (x_hi, x_lo) = {
                    let x = lhs as $double * rhs as $double;
                    ((x >> <$single>::BITS) as $single, x as $single)
                };
                // FIXME: use `mul_hi()`
                // y = x n nn = x (mod r) => yl = x_lo
                let y_hi = ((x_lo.wrapping_mul(self.inv_n) as $double * self.n as $double)
                    >> <$single>::BITS) as $single;
                // x - y = 0 (mod r), x - y = x (mod n) => z = x inv_r (mod n)
                let (z, b) = x_hi.overflowing_sub(y_hi);

                // x < n r, y < n r => |z| < n
                if b {
                    z.wrapping_add(self.n)
                } else {
                    z
                }
            }
        }

        impl<'a> Modulo<'a, $single> {
            /// Returns value.
            ///
            /// # Example
            ///
            /// ```
            /// use lib_modulo::Context;
            ///
            /// let n = 101;
            #[doc = concat!("let ctx = Context::<", stringify!($single), ">::new(n);")]
            ///
            /// let n = ctx.modulo(99);
            ///
            /// assert_eq!(n.get(), 99);
            /// assert_eq!(n.modulus(), 101);
            /// ```
            #[inline(always)]
            pub const fn get(&self) -> $single {
                self.ctx.mul(self.value, 1)
            }

            /// Returns modulus.
            ///
            /// # Example
            ///
            /// ```
            /// use lib_modulo::Context;
            ///
            /// let n = 101;
            #[doc = concat!("let ctx = Context::<", stringify!($single), ">::new(n);")]
            ///
            /// let n = ctx.modulo(99);
            ///
            /// assert_eq!(n.get(), 99);
            /// assert_eq!(n.modulus(), 101);
            /// ```
            #[inline(always)]
            pub const fn modulus(&self) -> $single {
                self.ctx.n
            }

            /// Returns `true` if `self` is `0`.
            ///
            /// # Example
            ///
            /// ```
            /// use lib_modulo::Context;
            ///
            /// for n in (1..100_000).step_by(2) {
            #[doc = concat!("    let ctx = Context::<", stringify!($single), ">::new(n);")]
            ///     assert!(ctx.modulo(0).is_zero());
            /// }
            /// ```
            #[inline(always)]
            pub const fn is_zero(self) -> bool {
                self.value == 0
            }

            /// Raises self to the power of `exp`, using exponentiation by squaring.
            ///
            /// # Time complexity
            ///
            /// *O*(log `exp`)
            ///
            /// # Example
            ///
            /// ```
            /// use lib_modulo::Context;
            ///
            /// let n = 12_345;
            #[doc = concat!("let ctx = Context::<", stringify!($single), ">::new(n);")]
            ///
            /// let mut pow10 = 1;
            /// for i in 0..1_000 {
            ///     assert_eq!(ctx.modulo(10).pow(i).get(), pow10);
            ///     pow10 = pow10 * 10 % n;
            /// }
            /// ```
            #[inline]
            pub const fn pow(mut self, mut exp: $single) -> Self {
                // r inv_r = 1 (mod n)
                let mut result = self.ctx.modulo(1).value;

                while exp > 0 {
                    if exp & 1 == 1 {
                        // n < r
                        result = self.ctx.mul(result, self.value)
                    }

                    exp >>= 1;
                    // n < r
                    self.value = self.ctx.mul(self.value, self.value)
                }
                self.value = result;

                self
            }

            /// Calculates modular inverse of `self` by extended binary GCD algorithm.
            ///
            /// # Time Complexity
            ///
            /// *O*(log `self`)
            ///
            /// # Example
            ///
            /// ```
            /// use lib_modulo::Context;
            ///
            /// // 998_244_353 is a prime number, so modular inverse of n exists iff n != 0 (mod 998_244_353)
            #[doc = concat!("let ctx = Context::<", stringify!($single), ">::new(998_244_353);")]
            ///
            /// for n in 1..500_000 {
            ///     let n = ctx.modulo(n);
            ///     assert!(n.inv().is_some_and(|i| (i * n).get() == 1));
            /// }
            /// // 0 n = 0 != 1 for any integer n
            /// assert!(ctx.modulo(0).inv().is_none());
            /// ```
            #[inline]
            pub const fn inv(self) -> Option<Self> {
                let mut a = self.get();
                let Self { ctx, .. } = self;

                // performs extended binary gcd
                //
                // invariants: a = [a] x,  b = [a] y (mod n) where [a] is initial value
                let mut b = ctx.n;
                let mut x = ctx.modulo(1).value; // 1 r mod n
                let mut y = 0; // 0 r mod n
                let frac_1_2 = ctx.modulo(ctx.n.div_ceil(2));

                while a > 0 {
                    x = ctx.mul(x, frac_1_2.pow(a.trailing_zeros() as $single).value);
                    a >>= a.trailing_zeros();

                    if a < b {
                        (a, b) = (b, a);
                        (x, y) = (y, x);
                    }
                    a -= b;
                    let (diff, b) = x.overflowing_sub(y);
                    x = if b { diff.wrapping_add(ctx.n) } else { diff };
                }

                // b = gcd([a], [b])
                if b == 1 {
                    Some(Self { value: y, ctx })
                } else {
                    None
                }
            }
        }

        impl<'a> Add for Modulo<'a, $single> {
            type Output = Self;

            #[inline(always)]
            fn add(mut self, rhs: Self) -> Self {
                let (sum, b) = self.value.overflowing_add(rhs.value);
                self.value = if b || sum >= self.ctx.n {
                    sum.wrapping_sub(self.ctx.n)
                } else {
                    sum
                };

                self
            }
        }

        impl<'a> Sub for Modulo<'a, $single> {
            type Output = Self;

            #[inline(always)]
            fn sub(mut self, rhs: Self) -> Self {
                let (diff, b) = self.value.overflowing_sub(rhs.value);
                self.value = if b {
                    diff.wrapping_add(self.ctx.n)
                } else {
                    diff
                };

                self
            }
        }

        impl<'a> Mul for Modulo<'a, $single> {
            type Output = Self;

            #[inline(always)]
            fn mul(mut self, rhs: Self) -> Self {
                // n < r
                self.value = self.ctx.mul(self.value, rhs.value);

                self
            }
        }

        impl<'a> Neg for Modulo<'a, $single> {
            type Output = Self;

            #[inline(always)]
            fn neg(mut self) -> Self::Output {
                // (x - x) r = 0 (mod n)
                self.value = if self.value == 0 {
                    self.value
                } else {
                    self.ctx.n - self.value
                };

                self
            }
        }
    };
}
montgomery_impl!(u64, u128);
montgomery_impl!(u32, u64);

impl<'a, U> AddAssign for Modulo<'a, U>
where
    Self: Add<Output = Self> + Copy,
{
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl<'a, U> SubAssign for Modulo<'a, U>
where
    Self: Sub<Output = Self> + Copy,
{
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}

impl<'a, U> MulAssign for Modulo<'a, U>
where
    Self: Mul<Output = Self> + Copy,
{
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primality_test_u32() {
        for x in 0..1000_000 {
            assert_eq!(miller_rabin(x), naive(x), "{x}")
        }

        fn naive(x: u32) -> bool {
            x > 1 && (2..=x.isqrt()).all(|d| x % d != 0)
        }

        fn miller_rabin(x: u32) -> bool {
            if x & 1 == 0 {
                return x == 2;
            } else if x == 1 {
                return false;
            }

            let (s, d) = {
                let x = x - 1;
                let s = x.trailing_zeros();
                (s - 1, x >> s)
            };

            let ctx = Context::<u32>::new(x);
            let one = ctx.modulo(1).value;
            let neg_one = (-ctx.modulo(1)).value;

            let mut i = 0;
            let a = [2, 7, 61];
            'test: while i < a.len() {
                let mut mint = ctx.modulo(a[i]);
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
    }
}
