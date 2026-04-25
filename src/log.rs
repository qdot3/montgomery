    /// Solves discrete logarithm problem and returns the *smallest* solution.
    ///
    /// Consider using [`FxHashMap`] or other fast hash maps.
    ///
    /// [`FxHashMap`]: https://docs.rs/rustc-hash/latest/rustc_hash/type.FxHashMap.html
    ///
    /// # Time complexity
    ///
    /// *O*(√`modulus`)
    ///
    /// # Example
    ///
    /// ```
    /// use lib_modulo::Modulus32;
    /// use std::collections::HashMap;
    ///
    /// let modulus = Modulus32::new(2025);
    /// let mut map = HashMap::new();
    /// let mut offset = 0;
    /// for d in 0..5000 {
    ///     let pow2 = modulus.residue(2).pow(d).get() as u32;
    ///     if pow2 == 1 {
    ///         offset = d;
    ///     }
    ///     assert_eq!(modulus.residue(2).log(pow2, &mut map), Some(d - offset));
    /// }
    /// // Since `5 + 2025 i` is multiple of 5, it is not power of 2, 3, or 7
    /// assert!(modulus.residue(2).log(5, &mut map).is_none());
    /// assert!(modulus.residue(3).log(5, &mut map).is_none());
    /// assert!(modulus.residue(7).log(5, &mut map).is_none());
    /// ```
    pub fn log<S>(self, rhs: u32, map: &mut HashMap<Raw32, u32, S>) -> Option<u32>
    where
        S: BuildHasher,
    {
        if rhs == 1 {
            return Some(0);
        } else if self.is_zero() {
            return None;
        }

        let mut offset = 1;
        let mut gcd = 1;
        let mut factor = self;
        // O(log n)
        while let Err(g) = factor.inv().map_err(|g| g as u32) {
            if g == gcd {
                break;
            }

            offset += 1;
            gcd = g;
            factor *= self;
        }

        if rhs % gcd != 0 {
            return None;
        }

        // solve `x^k = y (mod modulus)` by baby-step giant-step algorithm
        let modulus = Modulus32::new(self.modulus() as u32 / gcd);
        let x = modulus.residue(self.get() as u32);
        let y = modulus.residue(rhs) * modulus.residue(factor.get() as u32).inv().unwrap();

        let sqrt = (modulus.n as u32).isqrt() + 1;
        map.clear();
        map.reserve(sqrt as usize);

        {
            let mut lhs = modulus.residue(1);
            map.insert(lhs.into(), offset);
            for i in offset + 1..offset + sqrt {
                lhs *= x;
                // choose smaller
                map.entry(lhs.into()).or_insert(i);
            }
        }
        {
            if let Some(i) = map.get(&y.into()) {
                return Some(*i);
            }

            let mut rhs = y;
            let inv = x.inv().unwrap().pow(sqrt);
            for j in 1..sqrt {
                rhs *= inv;
                if let Some(i) = map.get(&rhs.into()) {
                    return Some(j * sqrt + i);
                }
            }
        }

        None
    }
