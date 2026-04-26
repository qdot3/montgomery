use lib_modulo::{primality_test, Modulus64, Raw64};
use rand::random_range;

struct RollingHash {
    modulus: Modulus64,
    // avoid self-reference
    base: Raw64,

    history: Vec<Raw64>,
}

impl RollingHash {
    pub fn new(source: &[u8], modulus: u64, base: u64) -> Self {
        let modulus = Modulus64::new(modulus);
        let base = modulus.residue(base);

        let mut history = Vec::with_capacity(source.len() + 1);
        // offset
        history.push(modulus.residue(0).into_raw());

        // hash prefixes
        for &c in source.iter() {
            let last = history.last().copied().unwrap();
            // `last` and `base` share the same modulus
            let next = last * base + modulus.residue(c as u64);
            history.push(next.into_raw());
        }

        Self {
            base: base.into_raw(),
            history,
            modulus,
        }
    }

    pub fn contains(&self, target: &[u8]) -> bool {
        let len = target.len();

        let base = self.base.into_residue(&self.modulus);
        let pow = base.pow(len as u64);

        // hash input
        let target = target.iter().fold(self.modulus.residue(0), |hash, c| {
            hash * base + self.modulus.residue(*c as u64)
        });

        (len..self.history.len()).any(|r| {
            // restore hash of windows
            // all of them share the same modulus
            let sub = self.history[r] - self.history[r - len] * pow;
            sub == target
        })
    }
}

fn main() {
    // generate prime at runtime for safety
    let prime = {
        let mut x = random_range(1 << 63..u64::MAX - (1 << 10)) | 1;
        while !primality_test(x) {
            x += 2
        }
        x
    };

    let rolling_hash = RollingHash::new(
        b"Rust is fast, safe and memory-efficient.",
        prime,
        random_range(2..prime),
    );

    // contains these
    assert!(rolling_hash.contains(b"Rust"));
    assert!(rolling_hash.contains(b"fast"));
    assert!(rolling_hash.contains(b"safe"));
    assert!(rolling_hash.contains(b"memory-efficient"));

    // does not contain these
    assert!(!rolling_hash.contains(b"slow"));
    assert!(!rolling_hash.contains(b"inconvenient"));
    assert!(!rolling_hash.contains(b"compilation is slow")); // 🤔
}
