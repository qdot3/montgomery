use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

macro_rules! montgomery_primitive_impl {
    ( $factory:tt, $mint:tt, $small:ty, $large:ty, $mr_set:expr, ) => {
        #[doc = concat!("Factory to produce [", stringify!($mint), "].")]
        #[derive(Debug, Clone, PartialEq, Eq, Hash)]
        pub struct $factory {
            // n * nn = 1 (mod r = 2^32)
            n: $large,
            nn: $large,
            // r^2 (mod n)
            r2: $large,
        }

        impl $factory {
            /// # Panic
            ///
            /// `n` should be a odd number.
            ///
            /// # Example
            ///
            /// ```
            #[doc = concat!("use montgomery_uint::", stringify!($factory), " as Montgomery;")]
            ///
            /// const FACTORY: Montgomery = Montgomery::new(101);
            ///
            /// let modulus = 103;
            /// let factory = Montgomery::new(modulus);
            /// ```
            pub const fn new(n: $small) -> Self {
                assert!(n & 1 == 1, "modulus should be odd number");

                // doubling
                let nn = {
                    // n * nn = 1 (mod 4)
                    let mut nn: $small = [1, 11, 13, 7, 9, 3, 5, 15][n as usize % 16 / 2];
                    // n * nn + 1 = 0 (mod 2^k)
                    // => (n * nn + 1)^2 = 0 (mod 2^2k)
                    // <=> n * (2 nn - n * nn^2) = 1 (mod 2^2k)
                    let mut d = <$small>::BITS.ilog2() - 2;
                    while d > 0 {
                        d -= 1;
                        nn = nn.wrapping_mul((2 as $small).wrapping_sub(n.wrapping_mul(nn)));
                    }
                    debug_assert!(
                        n.wrapping_mul(nn) == 1,
                        "n * nn = 1 (mod 2^32)"
                    );

                    nn as $large
                };

                let n = n as $large;
                let r2 = n.wrapping_neg() % n;

                Self { n, nn, r2 }
            }

            /// # Example
            ///
            /// ```
            #[doc = concat!("use montgomery_uint::", stringify!($factory), " as Montgomery;")]
            ///
            /// let factory = Montgomery::new(99);
            /// let mint_100 = factory.new_mint(4) * factory.new_mint(25);
            /// assert_eq!(mint_100.get(), 1);
            /// ```
            pub const fn new_mint(&self, x: $large) -> $mint<'_> {
                let x = if x >= self.n { x % self.n } else { x } * self.r2;
                let value = self.reduce(x);

                $mint { value, ctx: &self }
            }

            /// Performs Montgomery reduction.
            ///
            /// If `x` < `n r`, then returned value will be less than `n`
            const fn reduce(&self, x: $large) -> $large {
                const MASK: $large = <$small>::MAX as $large;
                const SHIFT: u32 = <$small>::BITS;

                // m = x nn n = x (mod r), m < n r
                let m = ((x & MASK) * self.nn & MASK) * self.n;
                // x - m = 0 (mod r), t = x rr (mod n)
                let (t, b) = (x >> SHIFT).overflowing_sub(m >> SHIFT);

                // |t| < n r
                if b { t.wrapping_add(self.n) } else { t }
            }

            /// Performs deterministic Miller-Rabin primality test.
            ///
            /// # Time Complexity
            ///
            /// *O*(log *x*)
            ///
            /// # Example
            ///
            /// ```
            #[doc = concat!("use montgomery_uint::", stringify!($factory), " as Montgomery;")]
            ///
            /// // composite numbers
            /// for n in [0, 1, 4, 6, 8, 9, !0] {
            ///     assert!(!Montgomery::primality_test(n));
            /// }
            /// // prime numbers
            /// for n in [2, 3, 5, 7, 10_007, 998_244_353] {
            ///     assert!(Montgomery::primality_test(n));
            /// }
            /// ```
            pub const fn primality_test(x: $small) -> bool {
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

                let factory = <$factory>::new(x);
                let one = factory.new_mint(1).value;
                let neg_one = factory.new_mint(x as $large - 1).value;

                let mut i = 0;
                let a = const { $mr_set };
                'test: while i < a.len() {
                    let mut mint = factory.new_mint(a[i]);
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

                        mint.value = mint.ctx.reduce(mint.value * mint.value);
                        if mint.value == neg_one {
                            continue 'test;
                        }
                    }

                    return false;
                }

                true
            }
        }

        /// Unsigned integer in Montgomery representation.
        ///
        #[doc = concat!("Two [", stringify!($mint), "]s with different modulus can interact but the results will be useless. Using block to drop [", stringify!($factory) ,"] is recommended")]
        #[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
        pub struct $mint<'a> {
            // value < ctx.t
            value: $large,
            ctx: &'a $factory,
        }

        impl<'a> $mint<'a> {
            /// Returns modulus.
            ///
            /// # Example
            ///
            /// ```
            #[doc = concat!("use montgomery_uint::", stringify!($factory), " as Montgomery;")]
            ///
            /// let factory = Montgomery::new(99);
            /// for i in 0..1000 {
            ///     assert_eq!(factory.new_mint(i).modulus(), 99);
            /// }
            /// ```
            pub const fn modulus(&self) -> $small {
                self.ctx.n as $small
            }

            /// Returns value.
            ///
            /// # Example
            ///
            /// ```
            #[doc = concat!("use montgomery_uint::", stringify!($factory), " as Montgomery;")]
            ///
            /// let factory = Montgomery::new(11);
            /// let mint = factory.new_mint(11 * 100 + 7);
            ///
            /// assert_eq!(mint.get(), 7);
            /// ```
            pub const fn get(&self) -> $small {
                self.ctx.reduce(self.value) as $small
            }

            /// # Time Complexity
            ///
            /// *O*(log `exp`)
            ///
            /// # Example
            ///
            /// ```
            #[doc = concat!("use montgomery_uint::", stringify!($factory), " as Montgomery;")]
            ///
            /// let n = 10001;
            /// let mut pow2 = 1;
            ///
            /// let factory = Montgomery::new(n);
            /// let mint = factory.new_mint(3);
            ///
            /// for d in 1..1000 {
            ///     pow2 = pow2 * 3 % n;
            ///
            ///     assert_eq!(mint.pow(d).get(), pow2)
            /// }
            /// ```
            pub const fn pow(mut self, mut exp: $small) -> Self {
                let mut res = self.ctx.reduce(self.ctx.r2);

                while exp > 0 {
                    if exp & 1 == 1 {
                        res = self.ctx.reduce(res * self.value);
                    }

                    exp >>= 1;
                    self.value = self.ctx.reduce(self.value * self.value);
                }

                self.value = res;
                self
            }
        }

        impl<'a> Add for $mint<'a> {
            type Output = Self;

            fn add(mut self, rhs: Self) -> Self::Output {
                self.value += rhs.value;
                if self.value >= self.ctx.n {
                    self.value -= self.ctx.n
                }

                self
            }
        }

        impl<'a> AddAssign for $mint<'a> {
            fn add_assign(&mut self, rhs: Self) {
                *self = *self + rhs
            }
        }

        impl<'a> Sub for $mint<'a> {
            type Output = Self;

            fn sub(mut self, rhs: Self) -> Self::Output {
                self.value = self.value.wrapping_sub(rhs.value);
                if self.value >= self.ctx.n {
                    self.value = self.value.wrapping_add(self.ctx.n)
                }

                self
            }
        }

        impl<'a> SubAssign for $mint<'a> {
            fn sub_assign(&mut self, rhs: Self) {
                *self = *self - rhs
            }
        }

        impl<'a> Mul for $mint<'a> {
            type Output = Self;

            fn mul(mut self, rhs: Self) -> Self::Output {
                self.value = self.ctx.reduce(self.value * rhs.value);

                self
            }
        }

        impl<'a> MulAssign for $mint<'a> {
            fn mul_assign(&mut self, rhs: Self) {
                *self = *self * rhs
            }
        }

        impl<'a> Neg for $mint<'a> {
            type Output = Self;

            fn neg(mut self) -> Self::Output {
                if self.value > 0 {
                    self.value = self.ctx.n - self.value
                }

                self
            }
        }
    };
}
montgomery_primitive_impl!(Montgomery32, Mint32, u32, u64, [2, 7, 61],);
montgomery_primitive_impl!(
    Montgomery64,
    Mint64,
    u64,
    u128,
    [2, 325, 9375, 28178, 450775, 9780504, 1795265022],
);
