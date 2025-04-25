use crate::dsa::basics;

pub fn ntt(w: &[i32; 256]) -> [i32; 256] {
    let mut w_hat = *w;
    let mut m = 0;
    let mut len = 128;
    while len >= 1 {
        let mut start = 0;
        while start < 256 {
            m += 1;
            let z = basics::ZETAS[m];
            for j in start..start + len {
                let t = basics::red_mul(z, w_hat[j + len]);
                w_hat[j + len] = (w_hat[j] - t) % basics::Q;
                w_hat[j] = (w_hat[j] + t) % basics::Q;
            }
            start += len << 1;
        }
        len >>= 1;
    }
    // At the end, go through and turn everything +ve, as required by the definition of mod in
    // fips-204
    for x in w_hat.iter_mut() {
        if *x < 0 {
            *x += basics::Q;
        }
    }
    w_hat
}

pub fn inv_ntt(w_hat: &[i32; 256]) -> [i32; 256] {
    let mut w = *w_hat;
    let mut m = 256;
    let mut len = 1;
    while len < 256 {
        let mut start = 0;
        while start < 256 {
            m -= 1;
            let z = -basics::ZETAS[m];
            for j in start..start + len {
                let t = w[j];
                w[j] = (t + w[j + len]) % basics::Q;
                w[j + len] = (t - w[j + len]) % basics::Q;
                w[j + len] = basics::red_mul(z, w[j + len]);
            }
            start += len << 1;
        }
        len <<= 1;
    }
    for x in w.iter_mut() {
        *x = (((*x as i64) * 8347681) % (basics::Q as i64)) as i32;
        // fips-204 wants all +ve answers for its inv_ntt().
        if *x < 0 {
            *x += basics::Q;
        }
    }
    w
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{RngCore, SeedableRng};

    #[test]
    fn test_zero_ntt() {
        let f: [i32; 256] = [0; 256];
        let f_hat = ntt(&f);
        let f_prime = inv_ntt(&f_hat);
        assert_eq!(f, f_prime);
    }

    #[test]
    fn test_ntt_random() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(56);
        for _ in 0..64 {
            let mut f: [i32; 256] = [0; 256];
            for g in f.iter_mut() {
                *g = (rng.next_u32() % (basics::Q as u32)) as i32;
            }
            let f_hat = ntt(&f);
            let f_prime = inv_ntt(&f_hat);
            assert_eq!(f, f_prime);
        }
    }
}
