use lib_modulo::Modulus64;
use rand::random_range;

fn main() {
    let src = b"Rust is fast, safe and memory-efficient.";

    // 9th Mersenne number
    let modulus = Modulus64::new((1 << 61) - 1);
    let base = modulus.residue(random_range(2..(1 << 61) - 1));

    // hash[i] = hash(src[..i])
    let hash = {
        let mut hash = Vec::with_capacity(src.len() + 1);
        hash.push(modulus.residue(0));

        for (i, s) in src.into_iter().enumerate() {
            hash.push(hash[i] * base + modulus.residue(*s as u64));
        }

        hash
    };

    let contains = |key: &[u8]| {
        let hashed_key = key.iter().fold(modulus.residue(0), |hash, s| {
            hash * base + modulus.residue(*s as u64)
        });

        let coef = base.pow(key.len() as u64);
        for i in key.len()..hash.len() {
            if hashed_key == hash[i] - hash[i - key.len()] * coef {
                return true;
            }
        }
        false
    };

    // `src` contains these
    assert!(contains(b"Rust"));
    assert!(contains(b"fast"));
    assert!(contains(b"safe"));
    assert!(contains(b"memory-efficient"));

    // `src` does not contain these
    assert!(!contains(b"slow"));
    assert!(!contains(b"inconvenient"));
    assert!(!contains(b"compilation is slow"));
}
