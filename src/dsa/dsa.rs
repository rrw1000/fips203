// Utility functions for dsa

use crate::{
    dsa::{basics, convert, encode, matrix, ntt, sample},
    types::{Bytes, IntRange2Or4, KeyPair},
};
use anyhow::{Result, anyhow};
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

    pub fn ml_dsa_65() -> Self {
        Self {
            q: basics::Q,
            zeta: basics::ZETA,
            d: 13,
            tau: 49,
            lambda: 192,
            gamma_1: (1 << 19),
            gamma_2: (basics::Q - 1) / 32,
            k: 6,
            l: 5,
            nu: IntRange2Or4::Four,
            beta: 49 * 4, // tau * nu
            omega: 55,
            challenge_entropy: 225,
            repetitions_times_100: 510,
        }
    }

    pub fn ml_dsa_87() -> Self {
        Self {
            q: basics::Q,
            zeta: basics::ZETA,
            d: 13,
            tau: 60,
            lambda: 256,
            gamma_1: (1 << 19),
            gamma_2: (basics::Q - 1) / 32,
            k: 8,
            l: 7,
            nu: IntRange2Or4::Two,
            beta: 2 * 60, // tau*nu
            omega: 75,
            challenge_entropy: 257,
            repetitions_times_100: 385,
        }
    }

    pub fn key_gen(&self, seed: &[u8; 32]) -> Result<KeyPair> {
        // H is SHAKE256
        let mut seeds: [u8; 128] = [0; 128];
        {
            let mut xof = Shake256::default();
            xof.update(seed);
            let k_arr: [u8; 1] = [self.k as u8; 1];
            xof.update(&k_arr);
            let l_arr: [u8; 1] = [self.l as u8; 1];
            xof.update(&l_arr);
            xof.finalize_xof_into(&mut seeds);
        }
        let p = &seeds[0..32];
        let p_prime = &seeds[32..96];
        println!(
            "rho = {}, p_prime = {}",
            hex::encode(p),
            hex::encode(p_prime)
        );
        let big_k = &seeds[96..128];
        let a_hat = sample::expand_a(p, self.k, self.l)?;
        println!("a_hat after expand: \n {a_hat:?}");
        let (s_1, s_2) = sample::expand_s(p_prime, self.k, self.l, self.nu)?;
        let t = {
            let tmp_0 = s_1.ntt();
            let tmp_1 = a_hat.mul_ntt(&tmp_0);
            let tmp_2 = tmp_1.inv_ntt();
            tmp_2.add(&s_2)?
        };
        let (t_1, t_0) = t.power2_round(self.d);
        let bitlen = (convert::bitlen(basics::Q - 1) as i32) - self.d;
        let pk = encode::pk_encode(&p.try_into()?, t_1.as_vec(), bitlen as u32)?;
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

    pub fn sign_internal(&self, sk: &[u8], m_prime: &[u8], rnd: &[u8; 32]) -> Result<Bytes> {
        let sk_out = encode::sk_decode(sk, self.nu, self.d, self.k, self.l)?;
        println!("s1 = {:?}", matrix::Vector::from(sk_out.s1.clone()));
        let s1_hat = matrix::Vector::from(sk_out.s1).ntt();
        println!("s1_ntt = {s1_hat:?}");
        let s2_hat = matrix::Vector::from(sk_out.s2).ntt();
        println!("s2_ntt = {s2_hat:?}");
        let t0_hat = matrix::Vector::from(sk_out.t0).ntt();
        println!("t0_ntt = {t0_hat:?}");
        let a_hat = sample::expand_a(sk_out.p.as_bytes(), self.k, self.l)?;
        println!("a_hat = {a_hat:?}");
        // Compute the message representative. H is Shake256
        println!("tr = {}", hex::encode(sk_out.tr.as_bytes()));
        let mut mu: [u8; 64] = [0; 64];
        {
            let mut xof = Shake256::default();
            // this bit of the spec makes no sense as there is no canonical
            // bytewise representation for a bitstring.
            // leancrypto just hashes the bytes and we will do the same.
            xof.update(sk_out.tr.as_bytes());
            xof.update(m_prime);
            xof.finalize_xof_into(&mut mu);
        }

        println!("mu = {}", hex::encode(&mu));
        let mut p_prime_prime: [u8; 64] = [0; 64];
        {
            let mut xof = Shake256::default();
            xof.update(sk_out.big_k.as_bytes());
            xof.update(rnd);
            xof.update(&mu);
            xof.finalize_xof_into(&mut p_prime_prime);
        }
        println!("p_prime_prime = {}", hex::encode(&p_prime_prime));
        let mut kappa = 0;
        let mut z = matrix::Vector::default();
        let mut h = matrix::Vector::default();
        let c_len = (self.lambda >> 2) as usize;
        let mut c_tilde = vec![0; c_len];
        let mut use_zh = false;
        while !use_zh {
            let y = sample::expand_mask(&p_prime_prime, kappa, self.gamma_1, self.l)?;
            let w = {
                let tmp_0 = y.ntt();
                let tmp_1 = a_hat.mul_ntt(&tmp_0);
                tmp_1.inv_ntt()
            };
            let w1 = w.high_bits(self.gamma_2);
            {
                let mut xof = Shake256::default();
                xof.update(&mu);
                let w1enc = encode::w1_encode(w1.as_slice(), basics::Q, self.gamma_2)?;
                println!("w1enc = {:?}", hex::encode(w1enc.as_bytes()));
                xof.update(w1enc.as_bytes());
                xof.finalize_xof_into(&mut c_tilde);
            }
            println!("ctilde = {:?}", hex::encode(&c_tilde));
            let c = sample::sample_in_ball(c_tilde.as_slice(), self.tau)?;
            println!("c after SampleInBall {:?}", &c);
            let c_hat = ntt::ntt(&c);
            let cs1 = s1_hat.mul_vec(&c_hat).inv_ntt();
            let cs2 = s2_hat.mul_vec(&c_hat).inv_ntt();
            use_zh = true;
            z = y.add(&cs1)?.reduce();
            println!("z = {:?}", z);
            let r_0 = w.sub(&cs2)?.low_bits(self.gamma_2);

            if z.metric_inf() >= self.gamma_1 - self.beta
                || r_0.metric_inf() >= self.gamma_2 - self.beta
            {
                println!("Case 1 z inf {}", z.metric_inf());
                use_zh = false;
            } else {
                let ct0 = t0_hat.mul_vec(&c_hat).inv_ntt().reduce();
                println!("ct0 {ct0:?}");
                let tmp_0 = ct0.minus();
                let tmp_1 = w.sub(&cs2)?.add(&ct0)?.reduce();
                println!("tmp_1 {tmp_1:?}");
                h = tmp_0.make_hint(&tmp_1, self.gamma_2)?;
                println!("ct0 metric_inf = {ct0:?}");
                println!("hint = {h:?}");
                println!("h.count_ones = {} omega = {}", h.count_ones(), self.omega);
                println!("metric = {}, gamma_2 = {}", ct0.metric_inf(), self.gamma_2);
                if ct0.metric_inf() >= self.gamma_2 || h.count_ones() > self.omega as usize {
                    use_zh = false;
                }
            }
            kappa += self.l as i32;
        }
        z = z.mod_pm();
        let sig_struct = encode::Sig {
            c: c_tilde.as_slice(),
            z: z.as_vec(),
            h: h.as_vec(),
            gamma_1: self.gamma_1,
            lambda: self.lambda,
            l: self.l,
            w: self.omega,
        };
        encode::sig_encode(&sig_struct)
    }

    /// if ctx is longer than 255 bytes , we will return failure.
    /// rnd is randomness. Pass a vector of 0s for the deterministic variant.
    pub fn sign(&self, sk: &[u8], m: &[u8], ctx: &[u8], rnd: &[u8; 32]) -> Result<Bytes> {
        if ctx.len() > 255 {
            return Err(anyhow!("Context string must be <= 255 bytes long"));
        }
        let mut m_prime = Bytes::new();
        m_prime.as_vec_mut().push(0);
        m_prime.as_vec_mut().push(ctx.len() as u8);
        m_prime.append_slice(ctx);
        m_prime.append_slice(m);
        let sigma = self.sign_internal(sk, m_prime.as_bytes(), rnd)?;
        Ok(sigma)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Bytes32;
    use serde::Deserialize;
    use serde_json;
    use std::{fs::File, io::BufReader};

    #[derive(Deserialize, Debug, Default)]
    struct KeyTest {
        pub seed: String,
        pub sk: String,
        pub pk: String,
    }

    #[ignore]
    #[test]
    fn test_key_file_44() {
        let file = File::open("fips-vectors/dsa/dsa_44_keys.json").unwrap();
        let reader = BufReader::new(file);
        let contents: Vec<KeyTest> = serde_json::from_reader(reader).unwrap();
        let ml_dsa = MLDSA::ml_dsa_44();
        for v in contents {
            let m = Bytes32::from_hex(&v.seed).unwrap();
            let pk_expect = Bytes::from_hex(&v.pk).unwrap();
            let sk_expect = Bytes::from_hex(&v.sk).unwrap();
            let keypair = ml_dsa.key_gen(m.as_bytes()).unwrap();
            assert_eq!(pk_expect, keypair.public_enc);
            assert_eq!(sk_expect, keypair.private_dec);
        }
    }

    #[ignore]
    #[test]
    fn test_key_file_65() {
        let file = File::open("fips-vectors/dsa/dsa_65_keys.json").unwrap();
        let reader = BufReader::new(file);
        let contents: Vec<KeyTest> = serde_json::from_reader(reader).unwrap();
        let ml_dsa = MLDSA::ml_dsa_65();
        for v in contents {
            let m = Bytes32::from_hex(&v.seed).unwrap();
            let pk_expect = Bytes::from_hex(&v.pk).unwrap();
            let sk_expect = Bytes::from_hex(&v.sk).unwrap();
            let keypair = ml_dsa.key_gen(m.as_bytes()).unwrap();
            assert_eq!(pk_expect, keypair.public_enc);
            assert_eq!(sk_expect, keypair.private_dec);
        }
    }

    #[ignore]
    #[test]
    fn test_key_file_87() {
        let file = File::open("fips-vectors/dsa/dsa_87_keys.json").unwrap();
        let reader = BufReader::new(file);
        let contents: Vec<KeyTest> = serde_json::from_reader(reader).unwrap();
        let ml_dsa = MLDSA::ml_dsa_87();
        for v in contents {
            let m = Bytes32::from_hex(&v.seed).unwrap();
            let pk_expect = Bytes::from_hex(&v.pk).unwrap();
            let sk_expect = Bytes::from_hex(&v.sk).unwrap();
            let keypair = ml_dsa.key_gen(m.as_bytes()).unwrap();
            assert_eq!(pk_expect, keypair.public_enc);
            assert_eq!(sk_expect, keypair.private_dec);
        }
    }

    // Test vectors raided from leancrypto, which can give us intermediate values.

    #[test]
    fn test_sign() {
        let sk = Bytes::from_hex("e12a6ad723953484a37d8fe9f86572800fd4a19818f516a79fedff38eb90cf61d9519fb38ece5f34f039adef040d35ac5fcf166649a44cbc76e021186e475815415c0187194c169a1445ef48fa9450daee9fcf2f10f4a5a9a7dbd3c3aaa1f2298f9169c4cbfd108b95295089413eca985016d540268f1a6bfae637a4f3ab2e4113140d0c8591613685e1045119a3719aa484121665c0b021a1a8281888890997314c928513c6304c907024880d23449163108c59469222326e60266403b780cab6240a48490c918ddc20900000609132125ac680640209844889a0462ea244299a2608d3c029118708e3c489a4b65052086c20356d53a2701a374049344e1a18692297111ca14850264ec2884d4c2086020802000340043982c1226551a4901cc945d2440014170aca48491c120291143202257050a4418a84315c34499332240ca42011c02c8cb8705246609c006cc0289289b6200203865b140824a36c80146ca188204cb611dc12244b1009d0406004419218194000c30414346108026800a800604612212846d840610ca571224222a4266edc100559482e8280851922061c35201c46881c416461304d22270214b62920254c889444998225cc106a6302061b026c6146300b212a1c08608b362c8a908dd300711bb72de3a24c84304509a76541425094241149c06d92340104c9418ba24c003771838200a1324d0208508b346c09118554384e09842449967024b12961b24c138071223932a43445cbc6088c1041c0260a09186e0046055a367044246c89406e209340428671db080c18a26d12338e8a86645408051109051b238dc9460183169203940001b22c54c21044900949b285a180819224881bc164c92888e100059932098a84648a0602012508da468502a648980802c1009048021060202d5ba22d5b2400a0b2451a0650e0c2204ca084581072c40472d9424111a93090c224a1a065e4b88c220471822808a4008500a9449b18511bc54dd0320683960c18a3048aa04881000250b05002b0914418601c112650304141100993245110302c04065012228250c2050837000cc06490940d4286804ab4115344221bc32d60c420532620c14665581064839884a000914a862c1a414058304cc0248511406ee1404901366d98a210d32246a0280d44080810260892242d58b8659a326c13180822c20c0b358d24462c0b1172912845c83466d0b830dc024408015013b04de142488c089221050442c091c4342a95a7d42f5771f88b4829d85715693d82343af18c894adfa27ed1f0003993a6afc4b7fe46763b64398b348f260951d9ada16767e36fe39e2d39849ceeefc975de8da89dacffa178284544b575e29d8dd1680bd8959a5fbe5eb0e42e786d438db3e28852d0a883dc8415182d41742fe7a6374adb0df91c49d1cf62cebbf6e6ce2980dc2f8d8a0aa29f67a9e8e0edf86d241d7e384e4a5da1cdcdad46954c230d01673c4f4df7ff0e1e70998c682b32fe1169a05be363f8e5f025bdac27474972299d6317e7c6e811e5d6fb1bcd6a8e25310c2b626249b0e220e0120101f153e416ac624d8579a2fdfec5c522e4b1530f2957d3dbce698e751adb28eac0d66a284f6bb8d37092555d98298cbd39375dbb187f8f0436d5e27d73c86e15488689d81ec5ea9f04ca5aebc5b35dc0d8373cd315760eb1eb6ac98abc9b644c83b3fce4682b4f8b811d3bc68b434b3e88d3aacd817bb65719284801788ae74b55edd8dfe4a491d8a0c8ef49b145eb9824cd6ce40973daf14885a1fa3540f311f1e3e04c430a0908664d84d49e4aadffd08ec633229f937c8eb3c47c234ebbce8592775b0899b6fe0b02191813e4d2f48cc4072deff0a6be63ffbd7fa1bfdd6af2e4f7f1227a17f27e5d6d7825558b127dd3bdf3de658bfa6dc3b921e4df892dddacab5cc30be8a92ed1c75245225e5651b3a8b86162073163645bc38a132a4ae25b136b4178ebe5f53c78988b2aa6c13557d8105b22deb1a76fa76b3e97a4abc56dfec25fc6dc8c29b1aa8851bfbc0a831dbc211cc7acbff8f4c5aae6594d2c83c8eb7bc4c060ed8ec3190589d87b5891ea77521f7bbe8579756772399a7ae973ec8456dd3a54db95f286e342f05fc4f3dc8a5054ac89cbd846596dd259cd4f7a5cda835afc7cdcaade4e61523dff522eb67a34ad2d66c5b3a22932bf67c233ac3b570b443e54dc6da7cf2c0b1e71b3e6f4645c61dcaa59045b491aed4e300bb29d68df1d00b1f2eec2d565f50789cf461a094d63f53acb457afce3bdda9469b4f6571715fb3c5eec8073b181c27200dd0ef67a582ecf8bd2e813e60984edee0b789e6cb31ddd25eed29217d1699bf61c404ab78af0c69d63c52bc5e6f921b95bddc7ac1b07136d075186844ec59413945f1657c85a5706e2ca9f90bbc4e95ec42db9f0d1587b5080459a50a75e903f76b56707aad7a75e6cd6a5b683b76d0099b5846b0c99d33f624cd7bc8be9c5c7f88141944b20af30d2b8cfbdc18e10ba44e04021bd69b6eb6fae5f92005e6a14318990014478aa7c630d10a269010672ebcac9d77ed2116236805cfed8911bbe4914e2929b987b8e503fa7d23306897726e45e2986641efe4d841577a97fe46b5004d5669ceea56607cbc296fcba225a2afa45727426d6a4a3e40e104be3a1bd9cbf89495be587d462c0d39d42f6a12b7f89363dd3fc53d7eed988dac2f6dca94feb519afa1070e552ad7382ef8537594a5b312dac4f9396a2bffb62c4f33c656002eaac9c6a61f57ccd2c0fe245be79fa5635a0d62781b5d1f77f174a6c7f8b005fcec878ccf5beb480b0192d33e77e2aea5ca4ddc4f77dab811e2cf8b040d387f77b9d094a1dbb575b7627394c4fab9e3817d4bec68a3cb8b4d20cbd49b6781978fee898ca1be8ee19033e5a5d08b88df40cad57bcb1253df10b29b6061fee175847de7e59a8dd95c68def80de7a9f7d63306b15ed2b8c7853837f7191816d3080888d584fecbdcc0f5924afc4fae94ad3f6c7e606c8f432e492a737067d4c05627f4bea50e5b87f1007ef63b4f023c64f56d09402be829780db8574adcd6b54205040dade7b7eb5134956c56c56c5fece374344f1117bce4feff86b7797bd476ce2b5bf92a3f404e663d4af63887dc28d8efc896d289e1bf91c70c469099a92a66eab3e50cf316cf97cf69c923f27c637756be7e87b706fe0f83c6b69330adc0fe73d75d90c8b3c706b316f7e7e1733538d6dd1c59f3b5372e128f603100b7e3e5619d1b99844af0f39f3ab8304fb507a2a68ff3907d28b2bd2fe720d68bf66cdde5a573436347e755bbc09d66ca16729b7515b330b7beba68457ffda20dd957b197bb5fb6d74e61381ecb88b35272ef91f6befa354dffb4e0d96bc9a82bff16a61258fee2b35c6f63a311a00a769433edeea4be5b6aad9de5442bf9472cd9df4ec436734dfa9e9ddda90fdecb9e1d6dac0f71dd9c2fc6f52ed90fa8c0020617ab6a28f402ee877915e4642a268f917b7bf02bdc17b2c8b2c920a3924d97e3a3a814d2ef541c718b6c43f05ad2344d4a3e2f1911f35442a1c2e53ca057bd620be3a195714056054ded3fbf50bf46d317a203c7f").unwrap();
        let msg = Bytes::from_hex("7f9c2ba4e88f827d616045507605853ed73b8093f6efbc88eb1a6eacfa66ef263cb1eea988004b93103cfb0aeefd2a686e01fa4a58e8a3639ca8a1e3f9ae57e2").unwrap();
        let rnd =
            Bytes32::from_hex("0000000000000000000000000000000000000000000000000000000000000000")
                .unwrap();
        let expected_sig = Bytes::from_hex(
            "de46aafacd4e1599d935a5451c66cbd131ba53065d3796b50912406e1fadb4be6824c65d90e391f35bf63755d3be75b60dca1b8dacbb1bbd73d8669112e33926299e983f7c57e3f0ca9f0ae40dc7a821a179b5f23c956defd454ddee3af4382ede9e45d1ddf86c025cb2cd8dbf90aac0b176890e661abb7a2775f2d26333e421f74c47b18bb11bdd0d2bd4ffa53d69cdf63fe446636e14768474c6617fe7b3d7485914a93390b0d46958b114d92e77615308bbcba7084bdef46341876a635de08256dcbec6f7b1913759da81852b574e2770945ec30d42e5842590f1e7a3377f2a5942de906ca655a2d515120775f1fb806f05124e27d1b0cd0f0e50407d11c9b9808bd150d81eefffc8b60e76ecd4ddf6f47acc6b17537daeedc1862948d694c0acd5f054c877a730968567aca4ed9d83d3dde46b05edc63e08598634e20a8c96020c04ca5b81ecd8e6910b1edfdca8aaad615b5ab2fff53feb9e44af178fb9e924b9ebc0a5c76cc40b03e5a98973db89a1e3647f9ec823a1e458aa0f18827f279829c584e050f902effdeaa1908c49d36e1e89dd34686e8843b4ed154d3163bfb2380a2d2ebc6836e1307c69d3c8b1e1eadc91ef2091dba8e4f25c4e3542d8d2649dc184825925f7d3e4a771e6d48339f496893f17e0cdbadb7c599fdc029f089fc0c385e3aee13d454c52f03b783c1ea56620acf6d89936565e95fbd346e3969427e4cf8497c64a68a23248cef85d6e70a590000421e87c04e93e731f8628cdb359e2086ca5ac140ae6dfc8f452016a11ecfb35611d32f00c8ea8daf6cbd0a49886afbed965b7fb28bfd1b98cdb99bf0f89aa4d95c88d4e7798587159e9fb07b3cc2ee824c65dae57fb169b2cb9ef05203752d61ba148350fbfcebb45ee419a43f042d3249cceb6be339d1ee6b8c9b38c3ab717061be7a957b61eeafc77c24f14173172eb3b855271e50688e292387de1ca1037d92995ad0e658ff99c676fec63a822e079c61b75bd3d0e539451b2f28078698cb7f9aa4b54aa3117c0108eff16a66739ded15f16471462030b9e7829eb49ab42942fcb22e4f46e83e520fc949e802103ef6bdbb4007431cdfdfe138730c97c5b80d3961dd80de656bb33e43d6f33cf3e9a9050b585c9c8435cf9dfc7086ceebed5dc48ac8b48780ce1a8f02cb3b93b4a21c6e1cec3dc8629bfa0ea289478bad20bbd1f74d3e4c7d8644ca2f3eb4e3f451e020af2622a5bbd309c463a9c98d534a96222a2cff2c98b1da51bbcf319da1ed2e58870dc66c66bc7d7b6106390776dc128a9e36114b6bf4e67a5a8fb25a36aecbda605000deb33cb846827d8ff08159498d2d43d23ce3ddd79010033870d4b697ab77d6968e8e3ecfea9a3ba8da302758033fd536ffd5da1faa41d455bbfb2e26f40c9a13f0646ce07d5bbee86ac7f9288324224af7c28e497e239ddadadf8e4a9184245c0a773d32414b1921abbfeba390f567c041e49e3a62dd435c9646b0427805b5eefff11e4faa46c6e3d6bc4dbf28524c6a9bfaa39d59451fcfbe11d509140bb1b2d4e8a6b3138411f3907691a67707297a8d97017152d3e4ecf2fa133200ccc694eb3ac853f442d4e43791679d39aa5f03f6755ab962c39946e6126c8ba156b9ebab512f1e6484b2967475699ac977774c49ca89de2d0999d621cb6f3c0fcba0e30d520c7fb203243a8e754a6850cdeb3eea327c7a256a06e38c44514b7405494142a6a6582e696f9eebd4edc001cfc36a0e003cd11898044cad8e36a71c9602bf0f35958fdb0c6831d061dce2c7ea2be4fa98e68b3c026ca46417634830b28d3c44933161cb5a6497b1bf1e5fdd9f16e6537028ea4251653a2110d17edfe583141f8cf71dc7622c65b8d0e85dac76e23be0e07f119cbd2bbf637d5158571423dd03379ceb17e58ec0ed6b405a74fa9c13b1bcfffddb8dfba4bb01d4def9b12d64cce41192267d5dbc0c96eb1638e4ae6f4cac54a3b384ea486ad1a0667b9f8c3556b7eddbb25a45d1f60bfa0bdf8f0880b10b4d0367153c1679551e81faa4388fc4cb85d932c2f29b5d47c60faf1899912ccaf6780ad924381bcea3750c0a3da9214149cadf042d0940e83a45c73f7a3a7f42c5d3f95dd7b86edffb6272f42502fb40d7a6dff43fa6dde9fe28bc820dfb2d67af9ef3e415e384649312212978188a3a2641eb6a457803299d4d94a506a9d68d5b5a1a27766e860881d6b13d66bae39e58a9d3a2ec9228c848fa36bcc39818a5068119ff83372bbe42c0eef583a503ff6364fd0ceef1fed22900642229d87891ceb3870ebe853952779faee112d2960c8817423e678285d04b4742782b8c3626368e7790ef38e24844b906b58f76600c6ac1ad529f3216dd508964ce36b969eec45957463f9a9e6ffb4b90464929d9d82e125a60f4fd97d8e6bfcf32d7bcb3e73bf687c0293b81eec9926aae527ca1d7034cb21b7e4a6e7a0a58662f30526daae35cc87782ef508b32a129368e24e95caa2226803c977baf6123e5d0671b54cab76985f9697093d3d7105f1f7299293ca6bcd1b444190761d1729b8f8f7b53fac28a2e013a209815fb8028a35e302149e70e7390dcbe15154e9a212ecf6e42c1959652ebcada4ec6793473ec7043160768819b180809dd32cf4a12cb6401dffaeb3c88a94ebb88f9d11e1b789ebd942e8357a21c393364bca100def90fe0cda79f714de995709fab118f22311d02cef470b22a87b3cc473cf778bc1b0e225607177ad166e569630d4dc4e97d7854f3d9742adc54fc7756faaf86de8ef6dd4563f6857029db1f9fe5c102fc2420a1b8f688d73f37d5bbb6395ea3bfcad7f8dc374a6333de8fe9866d8b27c9354ae8f8babf44ffbd0b1e4dc4e6b91eebe9507934690b24237e9007057b2ef4a789cc11e5e744a440cf8bb7d91e32d091187f088f52b6070aa30e78acaeef6b458c577c21bb60dcf9867cf68a89f2e88405b37d302eca56774826a4daff6dcc2375e1cbc29769f8a5e0017ced02c6106f83c5ced56e3121ab1bf59781ded9d29b93f729b05b1a472d7a2241269a486c83e21eb49af88ab2a4e9cf0f55735ffa75ddca8963c3ae6c737e6bb916958a500ceb2db6f016ba603497bf8613f718de1052625ed590db96f0de2241b4aeeb38d3341a9067ac2cffd537f9855563446712916bcd905dfc9dafa52b7354aab4c18e8d33c953b81f499e46f88497879c1a9c46cd1412c4d2201d70d3f0c7f519414fa8e88f8b9a066683f1e187ab08cd4f0da128fad56b1908de897fe1c576b34ff745f1d50f3e7975a1617303f41424b4e57606e76777b7c96a5c1f508101519214d4e5c8288a3cbdee4f8f9092e3342474b5e616ea3b2b4b5b8cbd9dae91e293e465a6b79858c989fa6b3b5bbcedcdfff000000000000000013233548"
        ).unwrap();
        let ml_dsa = MLDSA::ml_dsa_44();
        let ctx: [u8; 0] = [];
        let result = ml_dsa
            .sign(sk.as_bytes(), msg.as_bytes(), &ctx, rnd.as_bytes())
            .unwrap();
        println!("result = {:?}", hex::encode(result.as_bytes()));
        assert_eq!(expected_sig, result);
    }

    #[test]
    fn test_keygen() {
        let m =
            Bytes32::from_hex("0000000000000000000000000000000000000000000000000000000000000000")
                .unwrap();
        let pk =
            Bytes::from_hex("ba71f9f64e11baeb58fa9c6fbb6e14e61f18643dab495b47539a9166ca0198131c44f826bbd56e34e55db5e5e2d733485e39ea260fc6000c5ea4ba80d3455cde53b46f34482aedfd5450fc2e1ba4f25d15f9c144242fb39bb52287189030c50498e1717b7c758b190a6748ea9aa3f7acaaf2c7cb526ed717c9f79aeb84214fa5cd8ded92a0c3fa1558810f12c7050a367708d196cd24e5af974904aed8e4ce8872e8696b0b7bca50e452cd7d30ea9a4adac0311d672c6bde8496240b07431463708895cd9bafc31632d7397649388fdafcbf7d305a3de9a495eca7433a8f83ba0f0b25c413c6e39c96eb7d691b34d37ce37f1eead1cf217e25ef34eecf3f7c60f84b8edfdde8405d4f832576c61ef98e0a2f28da187700953924f686b94614705bcf53d33fedd4348edddbdf28b5065e1f20775043e85cf931f829179363a1a7e7404a838ec00086b0976386fe637c98244757e3f769ddd4467471bfad670f9a05f8246ee50a7b1eaf87fc4069c3ae2aa2033258117792f0bcd49e083fd1bc7496abff29cc94e4868b21214ed316525399a610fbdd4a80e7c80715f29578e2a84bb40bdddbd9f47a11b6e7da118a1b658d359e8aef55eb46b5376b5b655979984a922beebfc59bcd600d5309dccd72dbf0787db8ba757b537c1eafd5c0f50ea4bc9583549e2829a42c28cac248c96d78124c47159b18aedd754aba17b19d430fb78f633ea9d26f54a9bd50f8d8f6b73594f828976e7ea09c53bbb9f11a56c9507fb89b9a5ebc037a37267a95f85b8d64ca97192b10a66f417b3f61fe9ca57130a48fd925eae2ab5502d571c8a51903c1d398f4c1f76a7e11743976afdbc697f23094a3cd761ff9685de32e09fb3c28add453490300bc7c89dc01780096071722945775f264e1b0623bcf4619c712c838761205d87691b75ef360196cbb9e9b92a0d4c4ed62326e5024d77510b8ee2c7426cc22eae209dc9f13bde6bf08f5e7181bd3b459450b451a51539a715c21d67dd330eb5970db00d9edbfb2822b036fa13bafeb86d8dc78866e3f8d43e53d78cca5595a6faf886b5dc112f1cf4adcfa875800d90b48883af97316fe1506873fc157e570eacbfd222868d14234101966afb6bf9940829253a953ada89fc756b6a849f70acb9838e69faa50bba75e3e89c2adb57e86d088ab9b04a28e670709172243ec5e0008a5ceaf3f8722f487302596ffd755ad1b82a49c34b3469515b46aa290cd86ee38ea7a9be3f103610335b531cca333ddfe32b14510f4b07ef95fc6684e8c454a92c10dbb5d59c7a7c63fb305fe881967d99e669eb632840582560bb403431d40f75a4954908482278292821f4ea91e42e78fa48caee3c836146dcfd738d117e92e9a15137d28e8e6a4b4622650cb413504cb3a335d44beec5746c1c294b1e8cb99cb608d928f8ce3563632c521f23d13c61a8f61c01df8c96c7360db4f3c68aa5d2fdd342a62ff3459c116389421ab43e8584c45882b50e6e4e96db6f0b8fde890d5dbfadcd88690b449e64240ddb2023747f308363e301aa77757169fc6150628d5920b5aa1ab1c8cbf44cb00e025d7879d72b479e3af5311c785725590da9c89b9fc3b8450769554eb44d203eba2bbaef9cad2237011c2ea44eff00f299a48ffe28ca93ddf85f76608242ef8d6cc24610a1e2078fcac4f9385c314905ecaa82e553916d94d1a7c1ec652aa08897083daa2ebb1775fbc471ae27777d7904ea9f1b92bcac3d8a3158426087b645b1108f0d65fec93789c053743ca14fd63d05e98b652df2b9c2ff9ce05f1940703ffb273f80e0e2732eca9960d981b4cfd3b7bb8045b3c3830546b9dd8db0d").unwrap();
        let sk = Bytes::from_hex("ba71f9f64e11baeb58fa9c6fbb6e14e61f18643dab495b47539a9166ca019813c5c1610a40774f0eba3334f8b5be56e87871b3c3a772c0720fa37666ae1735fde6bc38a1c35f8cf08e440924c9037197bb87fdc4646b86da5a0589a326cc0c0d950ff8b5a9ea4135eab8a93f80f0927e124046e25b2366aea25a6d1d0fef982104b809da1248e1404c11038599b04d61040ac0342823168008124511490ca0960dc2a8015aa8650b3025d494448b388c1024412240698a04040c2690e33469e22468d244000b8490d9423021918d09996420396082362190b84c02126e24b970842201e1a86444260e01170e82406c9b386d0a9350522225d2b82c1a2570901222d1b28022465201348040342141066041328209459052322983b8449312210ac18960b2719c906990829089124940408a5c104801879103b040d84405e4b0290a27901cc95058028c0c497249846d19a785e30030cbb221e208104ca66411436852946909393291a26dca1484a39211db942c01c66414248c1a15310a230ed81420a4268d434090c14286a410850c322dd028601bb50c52446121152421c604d83431d91271a23452033089543210839408a2a470d2448e84028d9146620b000ea0807022b50d0bc90100264de24480932065c9968ddaa44c42a8001a276cc3944da0302cdaa80423a38c63040e59c80020a96cc4065082002e5cb2446146510c1880183882022590a2b26021a04020492a52288693a6289b947102b96d49a02c09252818930910494c13c991a3882d50123163b400a4004e04b6681214260017009b920c51a86511448d44426900a56d04a9694394694cb86052262da13211cab0111b81901211644286899c1251dcc0411347689b2092d1269123965013b34c12456921a1050c466599383099262e19157211398d89109250281148a64952b62dc9304dd9329191a64414340e13394198060c22256c14477021b760239825e4246ed1b64451442c414285e04442894242c2242cccc66813126998304923286a14c81160923011265292b0005c842c4a14290a9540e3c4050cc889a126624410049ac8700c936090b63012177123184ed410858420825c144890b251a018929394801a4972c00206e1a68ca1804d4a260818076242124010416444424e9486409a844121c57101962919404e94362224c32022887184204c642069db842424254621c4718a14640120654c84244424308a2270d386300bb46c21230a8bb84d5cf7d7e689302bedb1c5867e7d269b1cdb07f825641082e19a8da2f93077e8b1fc3d4e6b2d3258336b4f9c6455153ac040a847fb647fbb6b552a400071fe1772485b7a9d1f0d147bf3388c565471e4e62cc3ce0d0c0fc360df9289ed9918376b8b8b93145047f8fea2986007c2aa89922f69eb475b597b2bba237b9c842e3ff1d325e82a1f23e94989d006bc7ce4946f2e8b77e10848463c47fe7b209e2a617ddd41796ae6145e709cda9406f2261257c213b4b30da30ac25b0d06cf79a812c5fcb0ef11d9fedfe0994afe3b69b06a2916cf692b9da76028e5f3a04879e696d21f735c378315364db0a4e0ab6b53d31efaf30d65e37a1b6a77046f04c64ba1072a9780e0c566c94339a4d19d0068c57d6e6f0b512db7134a950eaf4f7b01a5fdd065b91bfa29e4423679cde74bc6a8f1c84c4df78387231dc85ce32670445903c4bebee3f50c43e50449496911aa93e7e395787414d31768d9912520f83c02ff01124dcf0e125fafd5b9d7e7dda4f5b50c70aebb8599a2e4476a0de531b0402672df7575142d8601605c94017923f64ac577c4bed8d8e89a74ca9f3819cbf142a72debe77c4efb7127e2d8c1b7bfb64286c0bd52233f43c67d5717f97ad82854873ddc7f71d656aaa6ef707060af280b9f454b4fedb4776e83b2fdba20a45aefeb549a1ed0382021893ca9a6e74ccc30a2553937ccef343899b502cf46ddb8dd1d95fefb60c9b20469a1503b2a687587830d33cee9a72d798fcf4a9b452c8549f559c5d9fc6bfe083f446c2d903981d9f26492483ab452ea5bb1008ffeac975da02759593e7e066361073a83b27b531a3d0dda517ca990ea3235d1d7b5e09da5f02dc1525b1da685965b54fc2a3a73a1790e0efb69e70a78fa550344ea8c753dbf18639baa8cb1259aa74f68f92aba8007c618ccb6f5069ff46b9751bbfff37df321360f0f5c0e7f5626dd129ae3ae2a7c56cdb611eda4c98fec83163cd5116878c1a93ebaa26db405eaf4a7aba277837de9a51504707624ef2e1bbbca292411167f2e3d390c0e51f84a2f138390e33f85835d38a94dbbe71e6c821e86b11ffd89eff4bfe208d6005d28f704baead1f25de0eb241b18fc7fa0ddd90dc139be7fcbeb9730fae4b5d17270ce4c670c42570a9cf25bc4fae5cd31e5d55ad0226a94be52948c6702a986a0adbfcd3ac482bb12abbb79a2f6602842153b2f82a3b3cd1688e74d36534bff8c48d3c451eb2c5f98feb9e7864d60af96e83b2162467482f058639c86a785a9a1d4b69bc30e77a64c3bbcd7deb4e3d30f1a6721203d87a88ab85e027a9742fc688f0adf15728e597e910cfe5df33c56a136ef39c7ca5d650c2b9f901c9b89e1e093549361f303be8839d1454cceb5fbc4435fa0dab58a8fc285360eea491ca077961c4aaa3e96de9971b94fdfa5207ccf0d9dab2c4896f07eb6771a383c6512f41ea28deee407fdae3c574f5d416a897a27ef7cf596f0432d624a2c4eace52f3cbf2c6331b80c9c9165bf1334246932024ec0be44b32136b4e434279135850364c757f1dcfa6385e2563312c5f553f0c844eabb7911cee760caeb3e193bf3a9c3811487239ad2e01478f46e418a5de56b7f1755ba68f9a374613b5de2ed26c580c772dbdbfab1f7e3f57d94f84e30deb29d70a91df288fc43a276dfed58e2b0db5383e532b6eedfb392e43dc3da7201a068f5231ee522098d6859b2d56463a8917b3c2561657966dbc47856b6ffc82bcc379ffd08b259f3d9d7873ba8fcbe4c9413b601159160701df00470b149bdf32f4d3cfcfb9debc77241717d13067aaed23c7a2651185169f1267061fb6b30e4fea73f66f4f92756ac2623418af8b2a398711b7c6807b43425e1d99bfdcd5df53195287906a332f59971a0c343975fc320ad137c9e34ce7ce85520b26ca197a1fa2df2ecd4e3fa833b3bd2c24482804252cf1df6adc6398f35e98ab18710407680c9c1dbac8c7edc8646b97082e2e121faa7fac21e4a338384cb92202c61bd126c5ddd458b327a18bd716f142ca5cdb3e49d7eb8d962b5b85a88f799b69a6a66c7bd629f56b43c0290629b5e274cdec7a07229e7939a77d32e8ef730fccead9c4e0677a83a0330ab765d336dd2aa155dcd2ac7f31529774f4936b05d0b14b48faa1e8cd45056e56c139b17f890715ad63d6c4a9f2d976c8b635bdfe58602816f612c6e4b225367cb9a7bb79c018f1b8c5315180aadbe3ab75ac356206fe27c12df3b569784e3a538fb052418266e72db400d6f32c4297f34f9f1af186c3765655f11b3e5a3c8049b7df14011ff215fbf17bf89ee976cf0dbab6270104e7e319d1f64c59e209e3582").unwrap();

        let ml_dsa = MLDSA::ml_dsa_44();
        let keypair = ml_dsa.key_gen(m.as_bytes()).unwrap();
        println!("public_enc {}", hex::encode(keypair.public_enc.as_bytes()));
        println!(
            "private_dec {}",
            hex::encode(keypair.private_dec.as_bytes())
        );
        assert_eq!(keypair.public_enc, pk);
        assert_eq!(keypair.private_dec, sk);
    }
}
