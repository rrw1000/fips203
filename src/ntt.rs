use crate::basics;
use anyhow::Result;

pub fn ntt(f: &[u32; 256]) -> Result<[u32; 256]> {
    let mut f_hat = *f;
    let mut i = 1;
    for len_bit in (1..=7).rev() {
        let len = 1 << len_bit;
        for start in (0..256).step_by(2 * len) {
            let zeta = basics::zeta_mod(i);
            i += 1;
            for j in start..(start + len) {
                // @todo can probably remove some of the mods here, but log_2 Q = 24, so not many.
                let t = (zeta * f_hat[j + len]) % basics::Q;
                f_hat[j + len] = basics::sub_mod(f_hat[j], t, basics::Q);
                f_hat[j] = (f_hat[j] + t) % basics::Q;
            }
        }
    }
    Ok(f_hat)
}

pub fn inv_ntt(f_hat: &[u32; 256]) -> Result<[u32; 256]> {
    let mut f = *f_hat;
    let mut i = 127;
    for len_bit in 1..=7 {
        let len = 1 << len_bit;
        for start in (0..256).step_by(2 * len) {
            let zeta = basics::zeta_mod(i);
            i -= 1;
            for j in start..(start + len) {
                let t = f[j];
                f[j] = (t + f[j + len]) % basics::Q;
                f[j + len] = (zeta * basics::sub_mod(f[j + len], t, basics::Q)) % basics::Q;
            }
        }
    }
    f.iter_mut().for_each(|x| *x = (*x * 3303) % basics::Q);
    Ok(f)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{RngCore, SeedableRng};

    #[test]
    fn test_zero_ntt() {
        let f: [u32; 256] = [0; 256];
        let f_hat = ntt(&f).unwrap();
        let f_prime = inv_ntt(&f_hat).unwrap();
        assert_eq!(f, f_prime)
    }

    #[test]
    fn test_one_ntt() {
        let mut f: [u32; 256] = [0; 256];
        f[0] = 1;
        let f_hat = ntt(&f).unwrap();
        let f_prime = inv_ntt(&f_hat).unwrap();
        assert_eq!(f, f_prime)
    }

    #[test]
    fn test_trivial_ntt() {
        let mut f: [u32; 256] = [0; 256];
        f[8] = 3312;
        let f_hat = ntt(&f).unwrap();
        let f_prime = inv_ntt(&f_hat).unwrap();
        let f_prime_hat = ntt(&f_prime).unwrap();
        let f_prime_prime = inv_ntt(&f_prime_hat).unwrap();
        //println!("f = {f:?}");
        //println!("f_hat = {f_hat:?}");
        //println!("f_prime =  {f_prime:?}");
        assert_eq!(f, f_prime_prime)
    }

    #[test]
    fn test_ntt_reversibility() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        for _ in 0..64 {
            // Some NTTs twist on transform, so you need to do it twice to get the original.
            let mut f: [u32; 256] = [0; 256];
            for g in f.iter_mut() {
                *g = rng.next_u32() % basics::Q;
            }
            let f_hat = ntt(&f).unwrap();
            let f_prime = inv_ntt(&f_hat).unwrap();
            let f_prime_hat = ntt(&f_prime).unwrap();
            let f_prime_prime = inv_ntt(&f_prime_hat).unwrap();
            assert_eq!(f, f_prime_prime)
        }
    }
}
