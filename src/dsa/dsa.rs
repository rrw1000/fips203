// Utility functions for dsa

use crate::{
    dsa::{basics, convert, encode, sample},
    types::{Bytes, Bytes32, IntRange2Or4, KeyPair},
};
use anyhow::Result;
use sha3::{
    Shake256,
    digest::{ExtendableOutput, Update},
};

pub struct MLDSA {
    // q and zeta are constant - here for completeness
    pub q: i32,
    pub zeta: i32,
    pub d: i32,
    pub tau: i32,
    pub lambda: i32,
    pub gamma_1: i32,
    pub gamma_2: i32,
    pub k: u32,
    pub l: u32,
    pub nu: IntRange2Or4,
    pub beta: i32,
    pub omega: i32,
    pub challenge_entropy: u32,
    pub repetitions_times_100: u32,
}

impl MLDSA {
    pub fn ml_dsa_44() -> Self {
        Self {
            q: basics::Q,
            zeta: basics::ZETA,
            d: 13,
            tau: 39,
            lambda: 128,
            gamma_1: (1 << 17),
            gamma_2: (basics::Q - 1) / 88,
            k: 4,
            l: 4,
            nu: IntRange2Or4::Two,
            beta: 39 * 2, // tau * nu
            omega: 80,
            challenge_entropy: 192,
            repetitions_times_100: 425,
        }
    }

    pub fn key_gen(&self, seed: Bytes32) -> Result<KeyPair> {
        // H is SHAKE256
        let mut seeds: [u8; 128] = [0; 128];
        {
            let mut xof = Shake256::default();
            xof.update(seed.as_bytes());
            let k_arr: [u8; 1] = [self.k as u8; 1];
            xof.update(&k_arr);
            let l_arr: [u8; 1] = [self.l as u8; 1];
            xof.update(&l_arr);
            xof.finalize_xof_into(&mut seeds);
        }
        let p = &seeds[0..32];
        let p_prime = &seeds[32..96];
        let big_k = &seeds[96..128];
        let a_hat = sample::expand_a(p, self.k, self.l)?;
        let (s_1, s_2) = sample::expand_s(p_prime, self.k, self.l, self.nu)?;
        let t = {
            let tmp_0 = s_1.ntt();
            let tmp_1 = a_hat.mul_ntt(&tmp_0);
            let tmp_2 = tmp_1.inv_ntt();
            tmp_2.add(&s_2)?
        };
        let (t_1, t_0) = t.power2_round(self.d);
        let bitlen = (convert::bitlen(basics::Q - 1) as i32) - self.d;
        let pk = encode::pk_encode(&p.try_into()?, &t_1.as_vec(), bitlen as u32)?;
        // H is SHAKE256
        let mut tr: [u8; 64] = [0; 64];
        {
            let mut xof = Shake256::default();
            xof.update(pk.as_bytes());
            xof.finalize_xof_into(&mut tr);
        }
        // @todo pass references here; copying is unnecessary.
        let sk_req = encode::SK {
            p,
            big_k,
            tr: &tr,
            s1: s_1.as_vec(),
            s2: s_2.as_vec(),
            t0: t_0.as_vec(),
        };
        let sk = encode::sk_encode(&sk_req, self.nu, self.d)?;
        Ok(KeyPair {
            public_enc: pk,
            private_dec: sk,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_keygen() {}
}
