use crate::{
    format,
    kem::{basics, matrix, ntt, sample},
    types::{Bytes, Bytes32, IntRange2To3, KeyPair},
};
use anyhow::{Result, anyhow};

/// A parameter set; intentionally no Default() - you need to select your encryption strength explicitly.
#[derive(Eq, PartialEq, Clone)]
pub struct ParamSet {
    // n = 256, q = 3329 always.
    pub k: u32,
    pub n1: IntRange2To3,
    pub n2: IntRange2To3,
    pub du: u8,
    pub dv: u8,
    pub required_rbg_bits: u32,
}

impl ParamSet {
    pub fn ml_kem_512() -> ParamSet {
        ParamSet {
            k: 2,
            n1: IntRange2To3::Three,
            n2: IntRange2To3::Two,
            du: 10,
            dv: 4,
            required_rbg_bits: 128,
        }
    }

    pub fn ml_kem_768() -> ParamSet {
        ParamSet {
            k: 3,
            n1: IntRange2To3::Two,
            n2: IntRange2To3::Two,
            du: 10,
            dv: 4,
            required_rbg_bits: 192,
        }
    }

    pub fn ml_kem_1024() -> ParamSet {
        ParamSet {
            k: 4,
            n1: IntRange2To3::Two,
            n2: IntRange2To3::Two,
            du: 11,
            dv: 5,
            required_rbg_bits: 256,
        }
    }

    /// Purely internal keygen function. Returns (ek,dk)
    fn k_pke_keygen(&self, d: &[u8; 32]) -> Result<(Bytes, Bytes)> {
        let mut random_bytes = Bytes::new();
        random_bytes.as_vec_mut().extend_from_slice(d);
        // We know it fits because max k == 4
        random_bytes.as_vec_mut().push(self.k as u8);
        let (p, sigma) = basics::g(&random_bytes.as_bytes())?;
        let mut big_n = 0;
        let mut a_matrix = matrix::SquareMatrix::new(self.k);
        for i in 0..self.k {
            for j in 0..self.k {
                let mut b = Bytes::new();
                b.as_vec_mut().extend_from_slice(p.as_bytes());
                b.as_vec_mut().push(j as u8);
                b.as_vec_mut().push(i as u8);

                a_matrix.set(i, j, &sample::sample_ntt(&b)?);
            }
        }
        let mut s = matrix::Vector::new(self.k);
        for i in 0..self.k {
            s.set(
                i,
                &sample::sample_poly_cbd(
                    &basics::prf(self.n1, &sigma.as_bytes(), big_n).as_bytes(),
                    self.n1,
                )?,
            );
            big_n += 1;
        }
        let mut e = matrix::Vector::new(self.k);
        for i in 0..self.k {
            e.set(
                i,
                &sample::sample_poly_cbd(
                    &basics::prf(self.n1, &sigma.as_bytes(), big_n).as_bytes(),
                    self.n1,
                )?,
            );
            big_n += 1;
        }
        let s_hat = s.ntt()?;
        let e_hat = e.ntt()?;

        let t_hat = a_matrix.compose_hat(&s_hat)?.add(&e_hat)?;
        let mut ek_pke = Bytes::new();
        let mut dk_pke = Bytes::new();
        for i in 0..self.k {
            let e_value = basics::byte_encode(&t_hat.at(i), 12)?;
            ek_pke.accumulate(e_value);
            let d_value = basics::byte_encode(&s_hat.at(i), 12)?;
            dk_pke.accumulate(d_value);
        }
        ek_pke.as_vec_mut().extend(p.as_bytes());
        Ok((ek_pke, dk_pke))
    }

    fn k_pke_encrypt(&self, ek: &[u8], m: &[u8], r: &[u8; 32]) -> Result<Bytes> {
        let mut big_n = 0;
        let mut ek_pos = 0;
        let mut t_hat = matrix::Vector::new(self.k);
        for i in 0..self.k {
            // @todo clean this up - needless copying.
            let b_val = &ek[ek_pos..ek_pos + 384];
            t_hat.set(i, &basics::byte_decode(&b_val, 12)?);
            ek_pos += 384;
        }
        let p = Bytes::from_bytes(&ek[ek_pos..ek_pos + 32]);

        let mut a_hat_matrix = matrix::SquareMatrix::new(self.k);
        for i in 0..self.k {
            for j in 0..self.k {
                let mut b = Bytes::new();
                b.as_vec_mut().extend_from_slice(p.as_bytes());
                b.as_vec_mut().push(j as u8);
                b.as_vec_mut().push(i as u8);
                a_hat_matrix.set(i, j, &sample::sample_ntt(&b)?);
            }
        }
        let mut y = matrix::Vector::new(self.k);
        for i in 0..self.k {
            y.set(
                i,
                &sample::sample_poly_cbd(&basics::prf(self.n1, r, big_n).as_bytes(), self.n1)?,
            );
            big_n += 1;
        }
        let mut e1 = matrix::Vector::new(self.k);
        for i in 0..self.k {
            e1.set(
                i,
                &sample::sample_poly_cbd(&basics::prf(self.n2, r, big_n).as_bytes(), self.n2)?,
            );
            big_n += 1;
        }
        let e2 = sample::sample_poly_cbd(&basics::prf(self.n2, r, big_n).as_bytes(), self.n2)?;
        let y_hat = y.ntt()?;
        let mut u = a_hat_matrix.compose_transpose_hat(&y_hat)?;
        u = u.inv_ntt()?;
        u.accumulate(&e1)?;
        let mu = basics::decompress_poly(basics::byte_decode(m, 1)?, 1);
        let mut v = t_hat.compose_transpose(&y_hat)?;
        v = ntt::inv_ntt(&v)?;
        format::accumulate_vec(&mut v, &e2, basics::Q);
        format::accumulate_vec(&mut v, &mu, basics::Q);
        let mut c1 = Bytes::new();
        for i in 0..self.k {
            let val = basics::byte_encode(&basics::compress_poly(u.at(i), self.du), self.du)?;
            c1.accumulate(val);
        }
        let c2 = basics::byte_encode(&basics::compress_poly(v, self.dv), self.dv)?;
        c1.accumulate(c2);
        Ok(c1)
    }

    fn k_pke_decrypt(&self, dk: &[u8], c: &[u8]) -> Result<Bytes32> {
        let expected_dk_len = (384 * self.k) as usize;
        let expected_c_len = 32 * ((self.du as usize) * (self.k as usize) + (self.dv as usize));
        if dk.len() != expected_dk_len {
            return Err(anyhow!(
                "Expected decryption key to be {expected_dk_len} bytes, got {}",
                dk.len()
            ));
        }
        if c.len() != expected_c_len {
            return Err(anyhow!(
                "Expected ciphertext to be {expected_c_len} bytes, got {}",
                c.len()
            ));
        }
        let (u_prime, tail) = matrix::Vector::decompress_from(c, self.k, self.du)?;
        let c2 = &c[tail..];
        let v_prime = basics::decompress_poly(basics::byte_decode(&c2, self.dv)?, self.dv);
        let s_hat = matrix::Vector::decode_from(dk, 12, self.k)?;
        let ntt_u_prime = u_prime.ntt()?;
        let tmp_0 = s_hat.compose_transpose(&ntt_u_prime)?;
        let inv_tmp_0 = ntt::inv_ntt(&tmp_0)?;
        let mut w = v_prime;
        format::subtract_vec(&mut w, &inv_tmp_0, basics::Q);
        let m = basics::byte_encode(&basics::compress_poly(w, 1), 1)?;
        Bytes32::try_from(m.as_bytes())
    }

    /// Alg 16 ; returns (ek,dk)
    fn keygen_internal(&self, d: &[u8; 32], v: &[u8; 32]) -> Result<(Bytes, Bytes)> {
        let (ek_pke, dk_pke) = self.k_pke_keygen(d)?;
        let mut dk = Bytes::new();
        dk.accumulate(dk_pke);
        dk.accumulate(ek_pke.clone());
        dk.accumulate_32(basics::h(&ek_pke.as_bytes())?);
        dk.append_slice(v);
        Ok((ek_pke, dk))
    }

    /// Alg 17: encaps_internal
    /// Returns (key, ciphertext)
    fn encaps_internal(&self, ek: &[u8], m: &[u8; 32]) -> Result<(Bytes32, Bytes)> {
        let mut mek = Bytes::from(m);
        mek.accumulate_32(basics::h(ek)?);
        let (k, r) = basics::g(&mek.as_bytes())?;
        let c = self.k_pke_encrypt(ek, m, &r.as_bytes())?;
        Ok((k, c))
    }

    /// Alg 18: decaps_internal
    /// Returns the shared secret key
    fn decaps_internal(&self, dk: &[u8], ct: &[u8]) -> Result<Bytes32> {
        // Write these out to avoid value dependencies on the indices.
        let k = self.k as usize;
        let dk_pke = &dk[0..384 * k];
        let ek_pke = &dk[384 * k..768 * k + 32];
        let h = &dk[768 * k + 32..768 * k + 64];
        let z = &dk[768 * k + 64..768 * k + 96];
        let m_prime = self.k_pke_decrypt(&dk_pke, &ct)?;
        let mut m_h = Bytes::new();
        m_h.accumulate_32(m_prime);
        m_h.append_slice(h);
        let (mut k_prime, r_prime) = basics::g(&m_h.as_bytes())?;
        let mut z_c = Bytes::from(z);
        z_c.append_slice(ct);
        let k_bar = basics::j(&z_c.as_bytes())?;
        let c_prime = self.k_pke_encrypt(&ek_pke, m_prime.as_slice(), &r_prime.as_bytes())?;
        if c_prime.as_bytes() != ct {
            // Implicit rejection.
            k_prime = k_bar;
        }
        Ok(k_prime)
    }

    /// Returns: (key, ciphertext)
    pub fn encaps(&self, ek: &[u8], m: &[u8; 32]) -> Result<(Bytes32, Bytes)> {
        // Type check
        let key_len = (384 * self.k) as usize;
        if ek.len() != key_len + 32 {
            return Err(anyhow!(
                "Expected {} EK bytes, got {}",
                key_len + 32,
                ek.len()
            ));
        }
        let mut ptr: usize = 0;
        for i in 0..self.k {
            let the_slice = &ek[ptr..ptr + 384];
            let dec = basics::byte_decode(the_slice, 12)?;
            let enc = basics::byte_encode(&dec, 12)?;
            if the_slice != enc.as_bytes() {
                return Err(anyhow!("Subkey {i} mismatch {the_slice:?} != {enc:?}"));
            }
            ptr += 384;
        }
        self.encaps_internal(ek, m)
    }

    pub fn decaps(&self, d: &[u8], c: &[u8]) -> Result<Bytes32> {
        let c_len = 32 * ((self.du as usize) * (self.k as usize) + (self.dv as usize));
        if c.len() != c_len {
            return Err(anyhow!(
                "Expected ciphertext length {c_len}, got {}",
                c.len()
            ));
        }
        let d_len = (self.k as usize) * 768 + 96;
        if d.len() != d_len {
            return Err(anyhow!("Expected d_k len {d_len}, got {}", d.len()));
        }
        let test = basics::h(&d[384 * (self.k as usize)..768 * (self.k as usize) + 32])?;
        let expect = &d[768 * (self.k as usize) + 32..768 * (self.k as usize) + 64];
        if test.as_bytes() != expect {
            return Err(anyhow!("Hash test failed"));
        }
        self.decaps_internal(d, c)
    }

    // Somewhat pointless, but the standard does it, so we should too :-)
    pub fn keygen(&self, d: &[u8; 32], z: &[u8; 32]) -> Result<KeyPair> {
        let (ek, dk) = self.keygen_internal(d, z)?;
        Ok(KeyPair {
            public_enc: ek,
            private_dec: dk,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils;

    struct TestSet {
        d: Bytes32,
        z: Bytes32,
        m: Bytes32,
        ek: Bytes,
        dk: Bytes,
        secret: Bytes32,
        ct: Bytes,
    }

    impl TestSet {
        pub fn from_hex(
            d: &str,
            v: &str,
            ek: &str,
            dk: &str,
            m: &str,
            secret: &str,
            ct: &str,
        ) -> TestSet {
            TestSet {
                d: Bytes32::from_hex(d).unwrap(),
                z: Bytes32::from_hex(v).unwrap(),
                m: Bytes32::from_hex(m).unwrap(),
                ek: Bytes::from_hex(ek).unwrap(),
                dk: Bytes::from_hex(dk).unwrap(),
                secret: Bytes32::from_hex(secret).unwrap(),
                ct: Bytes::from_hex(ct).unwrap(),
            }
        }
    }

    #[test]
    fn test_ref_512() {
        // This is the first vector from the kyber test repo, used so that I can check the intermediate values.
        let t = TestSet::from_hex(
            "7f9c2ba4e88f827d616045507605853ed73b8093f6efbc88eb1a6eacfa66ef26",
            "3cb1eea988004b93103cfb0aeefd2a686e01fa4a58e8a3639ca8a1e3f9ae57e2",
            "c29ac66c84bee3f129508c2b8c790c99a5ca41e5707e9b8c75c04d7ea8a481980a358b066a4a7e28d15d10374f75c33da6029cd490746fb55f5feab2ce823dec1830c0c18ef89ba3cc3bce4d252a07a40d401a4c8d273618b595db21ca4a959eb4d133556801b2b78254b2b5955c0400962dbb5a487aaa24a430b614e2af9338afd4b0339bd830a9cb761db1a83220352dd523b8a583022ca13e246506548c57ba3850ac5c50d864e2d48b89694277cc7a8174ad7b1ba2ee621d57d29b23c22df831872051145d2535dbb025e29c6da38cb475958e2a8808a431bed5b486821a1513a729cd979b6c8b382503cf53337616e1059ddc59977219d2f17b32acabbe86c468c1267b2b862ac565a7e232266c83fc700cda19b616377389428e62614cd6ac9f08fb12d21827feb324b14bc16cc0c20904b3973b9dfe60af577776449ac34eb00989e5a876143bc9b930a3c2a5bb861208f225dd97625bb36b2d8810d0452e458968b90c6a6d97c2dfca7341d53c741a60e2b91385c86bb3b1be5f10708c546ae6318209ba268dd2b4e3e4722761491ae215600248e7541e6080b865f44a169c02a5b374f1fa5502333753979c9da75452c4814a8b8fd682b8e0370699296bc2f42d38e5a2de92c5f5f61f24b65e1b9b75249b8759708706640000b490937683e1ccbbf39a0baecb1edd62826721524bcc4d3bf06f8e5ab368aa426eb42c701839e8bbcc69dbccd2757f84c60a409b6629d98dc92947367c52f0698c3c41c01110bafbeb0621cc175170375f8250c0f62406e32e10725b47e3c165845a7ada8fb2c55555f5cc8602246fe06f875a6eb5c9459dd19292428492a76738b8057e84c8e11436c4eb788e494e3952657ed331a2e4a47d75492365c473f484200502868c733bb9b2b3e44af2c98adfb91a4e1a144eec5a32d24f13c78775f206d7a94dee12007f778e484a1430666bfeda6dbb31b757d72b5c79a469522fe6599abb45b52cc6403c4970f65684bb770ca78893ba439057013655309cd84c524c9ca558e3612041b7b2026e2187598afb46f1d4ca85096dbc9bcc1c25779dfb607052e11649bb7f5f7268f979c4d8140afe6ce53830f38602290d751427f07b27",
            "7cf108c75a4d3592053d0ca79ce527ec1734f7023656e26253e2bfc68960dea73d28d7a9821597b48b4504837e27132c8ebc48505303aeb568f9d1928f7244cf98b4a88843c5db69845abfc4e40683dca3bde694445c63cc512bd2e4608e91a5697738bbca3b09a98da49757ad65c085125d11d01178d15acb251f417cbeff2a2ac1b05dc70839ffb297be205152e98371745de76243d3a302dfda7f16f5a8f0d2b73021a7f5c0490ff6b5a4125843889364288bede430835cce5bf8b86377132aea3ced911311a43b87e10ae471b49f1a31b636ce5d7c415d7c17539575724c3166ab23ca183f9962a0b17190f3dc204dc2551eec234b7c4e298a4e60f15d4e40bcaf1ba0eaa806e214bda0a0b0e0f51b34745e148290e222a4496c946766029cd82964429fd85a30ee3c0c6dea4e6bc840237c38d8e9b6a893702444bfeef783a104a8b46b8309a914cbaaa8cec9c08bb54a3cfc1a22c5a9d1b87685a30417d846d5e22242e8c0c8316cb7f15832b70b4bbc453bb67485e34c87d9a2d4b294d257a98d864afab445aea6bf4830404d170ca7484cc799b73ae3ce375c3a13d1318bd29158620328c8bc22d812d269c8ec133ce4f977bd95300054b5ef246bb086a6d3340575464183e10db1a126f81a43d2f800e962568de3411515b90dc1256446cf8fc33149e66a8657c221335eb57a8ca6262b06981a032823222c9da22835271033212b4e9d516132a428709b5882ab8764771832968dc73c4cd241adb145b5e654b38bab03e3524f113525dd0a63e8f195fceb26b8b03063c4863459cc10c5884cac016eebb7afa2519b27293b845e7e893e70bb21d0e004f1d05015dc71fcf6b150c750aa98c96f89c8da199eefc1712c409bc833a84d404ad72c873bf910c9190645b94dd2a31347596029d22abd071f5387588a9917fee4b33d1bbcfa1858390662b656c09ff824ed177abf44a118b7c78bba96b9f9b5a3a54e597829be686868e5463f99be58da0be7886501730e01e26aad0535896616d98125c38b3eb778114176b1bd498ebd92b5bd0a98a0a39f77b599c3e63e66fb62167b06c29ac66c84bee3f129508c2b8c790c99a5ca41e5707e9b8c75c04d7ea8a481980a358b066a4a7e28d15d10374f75c33da6029cd490746fb55f5feab2ce823dec1830c0c18ef89ba3cc3bce4d252a07a40d401a4c8d273618b595db21ca4a959eb4d133556801b2b78254b2b5955c0400962dbb5a487aaa24a430b614e2af9338afd4b0339bd830a9cb761db1a83220352dd523b8a583022ca13e246506548c57ba3850ac5c50d864e2d48b89694277cc7a8174ad7b1ba2ee621d57d29b23c22df831872051145d2535dbb025e29c6da38cb475958e2a8808a431bed5b486821a1513a729cd979b6c8b382503cf53337616e1059ddc59977219d2f17b32acabbe86c468c1267b2b862ac565a7e232266c83fc700cda19b616377389428e62614cd6ac9f08fb12d21827feb324b14bc16cc0c20904b3973b9dfe60af577776449ac34eb00989e5a876143bc9b930a3c2a5bb861208f225dd97625bb36b2d8810d0452e458968b90c6a6d97c2dfca7341d53c741a60e2b91385c86bb3b1be5f10708c546ae6318209ba268dd2b4e3e4722761491ae215600248e7541e6080b865f44a169c02a5b374f1fa5502333753979c9da75452c4814a8b8fd682b8e0370699296bc2f42d38e5a2de92c5f5f61f24b65e1b9b75249b8759708706640000b490937683e1ccbbf39a0baecb1edd62826721524bcc4d3bf06f8e5ab368aa426eb42c701839e8bbcc69dbccd2757f84c60a409b6629d98dc92947367c52f0698c3c41c01110bafbeb0621cc175170375f8250c0f62406e32e10725b47e3c165845a7ada8fb2c55555f5cc8602246fe06f875a6eb5c9459dd19292428492a76738b8057e84c8e11436c4eb788e494e3952657ed331a2e4a47d75492365c473f484200502868c733bb9b2b3e44af2c98adfb91a4e1a144eec5a32d24f13c78775f206d7a94dee12007f778e484a1430666bfeda6dbb31b757d72b5c79a469522fe6599abb45b52cc6403c4970f65684bb770ca78893ba439057013655309cd84c524c9ca558e3612041b7b2026e2187598afb46f1d4ca85096dbc9bcc1c25779dfb607052e11649bb7f5f7268f979c4d8140afe6ce53830f38602290d751427f07b27cda93dec4c4dc4d8484457fd882399c4b918c49fa8389a1dfa8c9f92f39b00cf3cb1eea988004b93103cfb0aeefd2a686e01fa4a58e8a3639ca8a1e3f9ae57e2",
            "35b8cc873c23dc62b8d260169afa2f75ab916a58d974918835d25e6a435085b2",
            "221d7d86011659313c83ce3fd0ab26797ef217e11d1f0bc76e7952fbe52a0a58",
            "5a645120b878936d202efc4851f38e6bb6573c3b14b0b9bb44bf372d8b1aa8034a9f1a1584076f0a38e89a9d49a50b792ace7584981be8e239272deef914418fefe2dad97dc0ec20cfe8a9599b9bbe3ecce91f97e10cd9ef2c4950e3ea3c46fe481eb0d24878c4624ad344f0dc9863e7d170937a8cecc6f7d00f9565529d572959cc49d0f7042ff43b7d1d71efd22f2654e14e78c31f34a26ae53b067ae0380a65a732459503da5e9406d50a70e3d5ebdbf3c9c01cad1cc001ebe69e6cff20e64ce5c802b691587e404cf6efa2799a2ffd353492a75f0e2ea52a974e0545a086a3bd14b69045238140a7200b10c3276cc6b2a67c173f7c1ad64545adb8ebdf7835e9aa1f54f6891369988f3625c45f2fec8d9a07b911d32ba69d9ff5d74f10a6808b25a3c81709945bf213c3450d74481f065042186b0d36fb55271162dfaf41e516a408f83ccabe8a0ba7effa16f88f6d7dbdb64e608c8f18d686c7e5d548d737116fa562dc76e7994a86374a9c85b8b17c4f025fa23a4de1a997a87e5a65f5c5386772491fd8d10731f5f5aa60366ffe3fd209cb7b7a8615320ad0728f41e812bd88d2d6104753917e89e1ca0f10177cf5dab040e466908b27446215709b0912972c428c4d9aad9432d9a159069c96154001cb0be4de597e9871b04bddaa4539f838bc12ab0a3ea7e8c8481bcdfcc1834369fcf061f7c599efdc4f6c434102991446aea12881e163fd4ee6c458b82e42759f8b11b0612c12d5a777acf4c7cdd26fb7da0b9098dc4af94704daa529945ab169cdb22d3966fcd26950e2418cac9bc7dc32c4a604f368f0f8a9c7ddce8b5e476b26b33116d607df1b49c205ada0d2ea5a5a64eecb22549ddf18a0daed2e5d44cb6174b9781236eee11f95ab0c45836bcafc73af3bec11440bc1c605669eb019cfac0097943cf29bbffce0f823293da623e5fa6d2a7ee0c7b4507596a62ef46bbe4e1b63cc96ba9878a7b39f84b59dd336f1659a24cbb33015d515e9e80e3e7902b1d583f8ee97153cd8ba1fadeb9ee7e2f2c23dfee85a50d5c3554b79922a4537d6dc4f09418dfdd744596cbf68",
        );
        let ml_kem_512 = ParamSet::ml_kem_512();
        let r = ml_kem_512.keygen(&t.d.as_bytes(), &t.z.as_bytes()).unwrap();
        assert_eq!(t.ek, r.public_enc);
        assert_eq!(t.dk, r.private_dec);
        let (key, ct) = ml_kem_512
            .encaps(&t.ek.as_bytes(), &t.m.as_bytes())
            .unwrap();
        assert_eq!(t.secret, key);
        assert_eq!(t.ct, ct);
        let decaps_key = ml_kem_512
            .decaps(&r.private_dec.as_bytes(), &ct.as_bytes())
            .unwrap();
        assert_eq!(key, decaps_key);
    }

    // I did think about unifying these, or at least making them common code
    // @todo really should in the future.

    #[ignore]
    #[test]
    fn vectors_test_768() {
        // Reads the file produced by the reference code's `test_vectors768` program and tests it against our code; see README for URLs.
        let vec = test_utils::read_test_vectors("fips-vectors/kem/test_vectors768.txt").unwrap();
        for (idx, t) in vec.iter().enumerate() {
            eprintln!("768: Test case {idx}");
            let ml_kem_768 = ParamSet::ml_kem_768();
            let r = ml_kem_768.keygen(&t.d.as_bytes(), &t.z.as_bytes()).unwrap();
            assert_eq!(t.ek, r.public_enc);
            assert_eq!(t.dk, r.private_dec);
            let (key, ct) = ml_kem_768
                .encaps(&t.ek.as_bytes(), &t.m.as_bytes())
                .unwrap();
            assert_eq!(t.secret, key);
            assert_eq!(t.ct, ct);
            let decaps_key = ml_kem_768
                .decaps(&r.private_dec.as_bytes(), &ct.as_bytes())
                .unwrap();
            assert_eq!(key, decaps_key);
            // Corrupt the ciphertext and check that we get an implicit reject.
            let mut bad_ct = ct.clone();
            bad_ct.as_vec_mut()[0] ^= 1;
            let bad_decaps_key = ml_kem_768
                .decaps(&r.private_dec.as_bytes(), &bad_ct.as_bytes())
                .unwrap();
            assert_ne!(key, bad_decaps_key);
            let known_bad_decaps_key = ml_kem_768
                .decaps(&r.private_dec.as_bytes(), &t.implicit_reject_ct.as_bytes())
                .unwrap();
            assert_eq!(t.implicit_reject, known_bad_decaps_key);
        }
    }

    #[ignore]
    #[test]
    fn vectors_test_1024() {
        // Reads the file produced by the reference code's `test_vectors768` program and tests it against our code; see README for URLs.
        let vec = test_utils::read_test_vectors("fips-vectors/kem/test_vectors1024.txt").unwrap();
        for (idx, t) in vec.iter().enumerate() {
            eprintln!("1024: Test case {idx}");
            let ml_kem_1024 = ParamSet::ml_kem_1024();
            let r = ml_kem_1024
                .keygen(&t.d.as_bytes(), &t.z.as_bytes())
                .unwrap();
            assert_eq!(t.ek, r.public_enc);
            assert_eq!(t.dk, r.private_dec);
            let (key, ct) = ml_kem_1024
                .encaps(&t.ek.as_bytes(), &t.m.as_bytes())
                .unwrap();
            assert_eq!(t.secret, key);
            assert_eq!(t.ct, ct);
            let decaps_key = ml_kem_1024
                .decaps(&r.private_dec.as_bytes(), &ct.as_bytes())
                .unwrap();
            assert_eq!(key, decaps_key);
            // Corrupt the ciphertext and check that we get an implicit reject.
            let mut bad_ct = ct.clone();
            bad_ct.as_vec_mut()[0] ^= 1;
            let bad_decaps_key = ml_kem_1024
                .decaps(&r.private_dec.as_bytes(), &bad_ct.as_bytes())
                .unwrap();
            assert_ne!(key, bad_decaps_key);
            let known_bad_decaps_key = ml_kem_1024
                .decaps(&r.private_dec.as_bytes(), &t.implicit_reject_ct.as_bytes())
                .unwrap();
            assert_eq!(t.implicit_reject, known_bad_decaps_key);
        }
    }

    #[ignore]
    #[test]
    fn vectors_test_512() {
        // Reads the file produced by the reference code's `test_vectors512` program and tests it against our code; see README for URLs.
        let vec = test_utils::read_test_vectors("fips-vectors/kem/test_vectors512.txt").unwrap();
        for (idx, t) in vec.iter().enumerate() {
            eprintln!("512: Test case {idx}");
            let ml_kem_512 = ParamSet::ml_kem_512();
            let r = ml_kem_512.keygen(&t.d.as_bytes(), &t.z.as_bytes()).unwrap();
            assert_eq!(t.ek, r.public_enc);
            assert_eq!(t.dk, r.private_dec);
            let (key, ct) = ml_kem_512
                .encaps(&t.ek.as_bytes(), &t.m.as_bytes())
                .unwrap();
            assert_eq!(t.secret, key);
            assert_eq!(t.ct, ct);
            let decaps_key = ml_kem_512
                .decaps(&r.private_dec.as_bytes(), &ct.as_bytes())
                .unwrap();
            assert_eq!(key, decaps_key);
            // Corrupt the ciphertext and check that we get an implicit reject.
            let mut bad_ct = ct.clone();
            bad_ct.as_vec_mut()[0] ^= 1;
            let bad_decaps_key = ml_kem_512
                .decaps(&r.private_dec.as_bytes(), &bad_ct.as_bytes())
                .unwrap();
            assert_ne!(key, bad_decaps_key);
            let known_bad_decaps_key = ml_kem_512
                .decaps(&r.private_dec.as_bytes(), &t.implicit_reject_ct.as_bytes())
                .unwrap();
            assert_eq!(t.implicit_reject, known_bad_decaps_key);
        }
    }

    // One of the official test vectors, to make sure we're all
    // implementing the same FIPS203 .. from 1024 this time, for variety.
    #[test]
    fn test_ml_decap_1024() {
        let dk = Bytes::from_hex("13E490ADD15409C73E2BC9B0496A9513FC89D3857C5E0A9661E8B8CDA9B3B0B81778863940AA9A970123BE3ACF07E97BCC397BF9E6C245029ED518BA64A20AA1E8154C3C736A1C17B43770AE695C38E96B6EE270DE3C81228497E74206272124E3112FCF7ABE938B9B8B51B874FC4362066559C6AA41E462F912026726377BF6681431BCC5D6AE9567230A5681E3192AA6586031DA7A41891C0FE61F9041C475D33570F2AEBFE03910C43B28E549E5C1C21BE2146352CAE4A66FF8D9ABD769963D93BB6F7C8E1856B2DB4311468645C3D8A824720D2FB93CA547A0A7BB15DB3B6A6E04367837AA378734A88064E033C9B8295099F336B1A80493FB8129F79640D263E943B6E535AAD7459F657A59A1CB12B19C8121F15310357D9125825D745F3C6B788DCA96F743C0F7DA9C033A6F2C8A9F998277F1B053756166DB866D2C2AC427E85AF361C7DED10D292A5F876C8B228A655871BA227C1A56B98B8B962F233C3A570A98A01600296249E64ABD84C70CFE91050CBA61D3588406E7A7A1CB199D07862D65C9F617433D173E557B385FD86EB3A642C75494E2273FC7B904E66B4E58B97CD82549A4635070224E5E039D76E7CC055872A356529E39A27E160D64570B5C92A9B9513A6670A227F7C9CEB3450F55B2AC50C1155AA1611A5DED59A069F43B9E3A9E4DC70F4D1164C830AD7BAB499491115D9201518B60CE8175F43327F8756E3DE2B34F653FC4933381F64EF9926D30686FF981C860A517040BB9EE140D0B541FFCE389F1F51DBB4CB51AC393DFEA5B6A55500D3C305EF88112985E623205516807AB1664CC8ABCC107583B8BB3E2E585F95B89AEDB45E5278C4D7B0960355307453E4A954B51B0A6E817CC16AB1D47839095798E6F7C53EB9BBDFDD0B5005C244774BF3AFA35810226E5D4C7C194C285CB8BD6DC2800963D2C1C7204B43B7290ABE6C17AED934D50F538AEDA48F28527DEF27C331B9F122C2941B884C2BC3E950C1A15A6AD5CD17623339C09872FB899245FB18EF879C5B53602C9436A15B79324753A9090851C3184F699CF75E0648AE1849E401385706521219DAFC8326CE8968BBCA21FC059E62963A532A096DA91457427C2A65B45C22E06038429B52B5F732C87C053E2AB38E67457F3EB967F402B7E0762AE4507AFE1196F22AE9E6797BE959D69DBA5B6FA6493E84EE8B3807FBABAD2C07180D8877F17CA92E08E0356CDCCFA9B9ED6268D19196E0BBD14693BFE72510F6151A9407338128449D3806C230B617525E9E01290AA21EB8BC8C70AC60028492A79B354407F5B8C786239368800C1B2D412DC6B14FA9226D1E765BE1780C0C346BCE1147D692A0BF195A1A2CACB75C1076B59222B214F682B8DEB7B7F51BF7BD9C775AA5474C4CE52963C546170B5539DC39C6ECB963D9C260D7253CF99F404DD49BFB7B4C3A90726B3692704910CFFA2AB8AE528D2A944E80735FE25917F960FC9D1014AAC200380A3591638ECBB44F5975708F584C0C392F2A20D2D1B698706853E28117C9139F303AF622286BE0176AB91812C962BA1C311AA9AA8F334C387424177EA6564EB318F832A68C11DA2418FB4BB6B74F1736955CAFE3170CB467F19F5CBEA7A8CC93692A3673CCD8450AA85860AA04BD4EB8008E63EDAE0A32C4413B3242AA9EB1E94371D8A02BA93173BBC93840EC2347560AED85A97D3DBBF0D7CA052635C794B824BD66798E01E70929FB1CBB81E865E0819749B92063EA99DED733751024A143339ACBA3BC9B065A19006E28B522D2A10B22B325E723E6804A2F19270C5744B1754A13B2128F1C3BE87D3968C8213E5B55442E6749FB21C49396B0B427C030A3D9F397525D8B8953B945C132E5B506B5130708A55C41D8B9F7C9B0781D5340D029E538378640A1A57340DE7697A8CD37B58C64FBD7518F55AA276650C94B40C71B435FC83927740222BA32BD8947B7FE3277D9301B1524ECB735EA358AB6BC41682728DA25B70D1AA2E35A4606BA27351A54C88353A3D455651C3118AA578C4B19C557A6DC1B290AB2C97CCF96665C10FBB67A935477168E0185ACC1155A0132011C4F2523732373BD5425AFA7C061BC66C7883803039296C989224A0BDE9BC4B8B6698EABC6C40A10A842C697B3873D84BBCE048B69DC705E6A7BB6C813BC1491219D429AEC6B9559BC714E17C27600E225AB2398223B9247952B73562B78E6F690561E62951B758AD815E3AC031933CBA45749EF836012456C3E4677347DB502A564A41295FF44449E4E5A53A60A0564A8FA062BCB071C42858A9FEB134BEC8281A57525A6A1D4008483BFBAAFEA86C3BD96BEBE0B842652553542406EC555F8C91AB1953B036A0EE5401C7B9727B7150B2FA50C328B7D5BA337457BB3A7936C0F2443EF816E9B881A7D15C7E6B5FDED996F7ABBBEB1435F70AAEB48B309422C95CD48B66E033BB8212DB82C39F0736E55C7722567CBDF89F1E682343F63ABE16ABF8742CBDBCAB3AF5845269776CEC36920B5B442C986AF62CD5124F10E61B9DB446079AC0A93625A3C15A9709CF01130CB4D721E32B03675B3073B5B3B6E387462399DB7B1F0D18311DEB57C2A33834F8B6BE3574F358A75E6905C0F68A391822E48762E126B62A61194A5707A01083060BA8CDD5376302794565BF647C2F506958C9B4AA0A539EFB859F15142DB0144558D5237067C67340A6BE56544B174DE7998E25E281D2621951ECB4DF063CAAE648DA5B87A69B3D06E779C38C9F8B443BA627398A5301156BA76C27414DD30B29393B77481153DA6AEA138C39AA9049219E4D6C73C2C28EE8F57BBCA2B43B3247E49736261CCA9C52987B0BB7242013BE083B5F770C97E78602A7B0FECB8DF584B86B5C1866FC9C1841C3D3DC17D45BC0B76A6792CB83B3983B553876CCB63DDC1720F782AE08E1A7BA329464E0B8345127867B34EAAC7911AA8DE5BC91FF7479F46105C9E0A5318A780C3156C7C40B2FB7ACD47969E0F7094966BAA7C5C211666B4784B5B10AB1C24677F863615763678BE15A7FBAA38D238253B92C96A76071B50B0CD0222FB4BD57D052CA695E8B3AAF2ED35C131C7EDECC04DB9C25A8EA266B34A4C2B446EC23AB0385501DD4516AB286965BB4322C854AF7CABCF2144B36881899AECC3CCA1D1868AEC77F63AB0EAA85AF792972231BC381EA561F4BA5F455128E90842E44C932F43EE27A0C495A6C784707FA4950C9715D2201AE5A04A0D8820576E322FED4AA20443844E727DD600D58E27691A879124570C858C050E076D003CFA28A149C517E67B58C9F53762149265B9AB37FA3BB32A14911A02FE91A4F8F70329C49A010ACCC0C1259A94A458AC0A99F219C6BBB0141393E09607B94B6138862AE0E809CB5910A48D7C92035496EB2731D70C90F866CE0988896E5944CF5CA000911B817621008244848907E9C75A0241519B08C13282A76B07441953BBCE7C9A55CB2E0442D20415C5455B9B5B6B98F34AC68C33ED0849D826C051A1B1221D08520E8A8285C983FF9025773B6265BBDAB802A03C68F21224A68A7A7FAD62D22AAC1E8E7A6AF04A4C05732EE522CC55587A0A55F1A145DD34724B062BE5EDC7E716360CE0341CB9B7751917207E228C5E452C1CB66043116CA375BEB5C9981219BA8551762B55BF562215F840202E6535903286E31AAE3085CC0A3BF92CC8C3C943DBA7876A5A8CB909CA91746BCCE402712705E8CB84890C983A62534E064A53D013188FB38D4EA69B07B75F04B315CA86E7E9A9E94089DCB9A0A7026CBA238010B48ABAEA51483D255ED7A5700C8C4CFA5BEE4189302C12DE462613BB74B776271D174CC7F3642F7F158F5E70EFDC62E48E36D215264F7227F57006C2FE2B6C4D02008CC85F5D21FB146714424457EFC1F826942B0824FFB30812D5570EBA48FCD365566E09D12E75D4593803F2BA841082FE65501FB2271A902286E3B20F7637F1DA87A779918EC022C85B1242ED78E3A79A3223B5D73C4582E430D46696EB16A3A22BAAD70F88CEC827A65CB3228B44638E0BDC7D361AD40405FB22A3599673B407D8B7B5BB24281042B834F140D3D5C8D3AB2029582C46388A33D370692E2677DA34ACD4A748C69130B8886AFA88BEDE95DD27122A4F4827172029B00708E701F0763B2B5F388FDC82FFAF57AF804725915C432A2B0D671673D41C3BB04C868B35C700045B9A056B7B2359772670C9AC5840C5BB0DC2B9FAC5B4ED12A20251A91085C87A8B24DBCC165D94307965EDE4B8FCE38158BE4BDF3473D160184F745053BD99410091C15BC1288F1130EA2277054CA106438EDBB4F3CF893B5A9007B8C192863BA52A17E39EFDC7060E0C8316936A289684F107DA9B710FC15EE55A9D7439071F7D26129D71267870D23C3FCBA53A325F5F4EFDD0818753D68930C0C2797E74DD2FA9550567AC47E70FE9D84F482BDAB1A8CC149AA7889E4ABF8ED76F68D624DE6").unwrap();
        let c = Bytes::from_hex("FC555C15C5FF48C0FA546CAB63C5573BBABF962AE5820147BBB59C4A62468EDBF8EDB240934C6DF94398D1871045ADD21466D0B0FBF74CBF75AB3F8E3C0EC6C6E1218EF1492430647F773EA0A0EABF330A7D4096F3BFD3B1BF6E8CC83996A9A979744FF0FD36F3E18F673DF7E25047BE569A2932C0249C08C5092D44F5C7DE206B2CF5F424B0DEA2086FBF9C1AD2D874774306CC6B072FF4ED2C379A97274567681AE91D77D34C341E7D34EFBBFE3BDDA51ADAB6A52747EEB0D7595031F9ACD419C19CE035964B6CDC0ABEE8D4ACEECF063203E3D5DF66A6A6991CAAEAC7F2C42183D0922587EF583FCEBBE45036CB58D74FF61BB5183DC8E03B0E57679E774B55DABEA668AF0F1C38F0123A1E54AD1381883ED2F9D9A0F5A5269C9EA0A63BC30955EE75FBC1AA2AD2478A1973B634C6324CEF49B4F7F4464EC46CCDB27A52672BAC50F235D20DB9967BCA8D78309EA629E4B6AC3B2BE473A8CBD0299663CD02E03F085CCA8454EEE972A95C5814A77F0707219C1F9B7D69493E158F0AFC377ABE55574A9E8D9AF986BBE0C56E7ED0AE054112CF065286434FC1FC4B539984BBA048E9C26090B54F3783979118171DEF512C1A7E12FB5BF97948BED2A0A64687F39D73E7AA0C6A6F9341850106FB4FF25C714EACCF9779C4F9BFF0438F2C5BC974DD1FCE177C78DA9F6DCD44D4F355903D993230F9FC94B85045A72C94022C09E6D1B1AC7838B1A59C4455432FEEB986A8997EC4F794617EAF77CC9A50FEEB6735F5E0DD59BBFA518D437438DCA65A1F8609847A2DA6E2ECE376748BD022894C7CA17B4E84624F4EDDB69D699DF646047D8D9CC09EE6C0BAC2751BEA96EEDD77C0593CB9F8172CB0E004334E9593919F498234F2F9B6473169A57CA797971BDE367B93BD08ABD942D3EF9EC80C135DFDC31EEFDAB028B382BB36C0467283F0FD9A2F5B50DE3E4025646CEB96667610BE2B1BA073B4D14E7227DB6A6027993DADB6243764E599DC658AE96667F2726D40E4F98DAA82396D2EC98B08323B79EDD50CD56D319AC72BAA082802400795F0E93A240E3934ACA381871356279488BFF55579F04459D9742F022BFBB59BB5F8438069FB62A40FFA9EA8620899E94E85E05E0096063DD4A578D3A8212AAA52003428CBB9543D7ADB3CA35FD1575D763A64E94593184B0DB7C5EAE9BC0F52E47A6E7F4219CDE78A4F043A88ECC09F69BE8B2D01362201D685AF0E7A0167961EF3942348AA0268C18CA6E7C3C934C989FE834617F21BD9BBB87CB7CDCCE0DFC7FC46F1E10427DE002D4B7B8A64BDD56A22EB57C5AC0622D5F86F00C4236AA9C5C8FF87EF966E501FF9D4243D49DBE52B8C41425DCD1633E87D7547A273B13AF9B4F07D919C0709C99AB0E7E63FBA737B063F4852E20B1A9637C8B9FDFAD631FE14DD57B6CF998EAEC7C847965235CEE1D189B39E704E3367369772FF8B9C6EE4388FDD8010F479615D228446327A3AA53A77CE1538096A2374809531C6E23078DA06258F68F511C07B4197F427CFF07EBFA00B187CD1D6D0AA940B84551969C7F0ED647430EA0F21ECAC1B18283FB0E243B8FA57B28A418BB2E537DA845634B1C6BE3149CB7606F5D5ACC6B1FB1AA840D22ED524A31C677852BC68D241C7F6F5CAC806B53576730C565C8345DFECE1DF94682A29C21B0981C51CD8FD1E1ADDC709BE7E153D78C693781803D4E2EAAF3E84A8645671A304EE06ACFCC629E6FE1C3ED57BCB1A8CB7D7E9C88FF73AC5A8DF216B2A307A80E7B108F32A0967956256ABC15A541B0BAADAC7F2EBE3B3F986644DBB36E673AC7548D9CADADF7878ABDA71DC684DFE24A1C293B06C473CAC1638397CC3BD1B7E55850735EBC97550DC1CDE193EBB178A5190ECF0310C6FB09FFA1576E3304F01B3E622AC897F6657FA2361F3573AF5130932BEAE4B62D42FC883E431C32A63DC5E785DBDE147BB769E08CFC29B96600C406657264CE3AD56DBDAEE0EB707E228300123AA41D6B6D5B3D2AC4A0B555644600B83BCE38B6F08623FE0080040A44E8F05168CBAA9E626125DF902AA4070C55DDED87CC95236078DF5E7A1D67A2358C336B92DF7979228FA0267751CA3F4E2C0A012110BF56FE3B08300C64E74D3C59D39C0DABBA33274ABA6B4A96B52E0AAF351B25DD253717FF350A8E82F9A6D3072E112B151FB9821BFE021358B3068EF88A5FF31D97248C1B62B158E").unwrap();
        let e_k =
            Bytes32::from_hex("C72D033632F2272E8DAC83A1E2494E9129695893D28BD39131D44F60A4380AFF")
                .unwrap();
        let ml_kem_1024 = ParamSet::ml_kem_1024();
        let k = ml_kem_1024.decaps(&dk.as_bytes(), &c.as_bytes()).unwrap();
        assert_eq!(e_k, k);
    }

    // One of the official test vectors, to make sure we're all implementing the same FIPS203 ..
    #[test]
    fn test_ml_encap_512() {
        let ek = Bytes::from_hex("A5F799D57B310740345CF77783B5013D540F557143443A5402B1255A5B0437727113E26B516C2BB899BF1178BE7531636E810B84938DF0B95197540A39289DC3C91CA3E8201A37101221922D5A2E59719F97375D30339196F10F7E986FDD4BC27E192FEE7654F85CAB2B01AF2E52AA5420295D6429CF5B93981AACEF634DD3B055F479B72FA45B012433A16939438641245C7113951E42A78399DB1B3451AC317552440322B93577D1C0A03BC02875F0B3E93A9A24E503DED4BF9B095F0023867122BFFB16785E25BAB9D19670797A5EA812CE22B7E1DB2BFED18F513625DA434E4D1A07827277386448EEAB09A7395C0CB780EA152989D429C1AC4187C21B901CB535298AE1753BC42C33BC009839E254C1D61C1CFB5ED4C34446CBCFAC33935FAC22019498E1F8610BE012BFD637D4D330EA688D384A2AAEBC58B6E389B1B8263DC773A3D989D4A768EFEB74A3643B523947881CABFF7A7A22839572CC45841147F4C25AD590B834C3B65A9B8F4F35988E36551FD7C9ABF9C0AAE225760744BF937B4DC7B701E86C99D0B87F024B03F9651A58A075C321840B51ADEA90502C9757272FF45C8EA9302013A0BF23864ADB55762332C59C73023D8699FD15CBD81185292996B67155ADFA5A6FB51904B626320417239317CFCA687F585F9CFABC78BC08C7128405486D3AA88181C7116B114FDD3004E337B5D25B6AF2B69D0AE1BDA696A688B7C5189CC7E5F941C397474CA09D8178B99E1BB2CA590A11E08BADF18F31A14E7E24CEDA79BB94D7062FB95555283E4C4232C0E5B9E1A970DBA176F9206F1BBBCB082ABB9D6457FEB87B9E977BA1FCC846B0317F084CC0890C7B3BB1A70541BEC77183D62440139EECCC9F543218A2757006133E0C27877BA06A5D691485010938A96E29B249D19C164EE07CDF67129A1136A82C9E8DD05CB8069300DBC78DD5192C385C2A005ED4D5CE22928A6DB6A68044511F4193A45796E04CC03E832EFD8AB15C4C5F16CB84BA848A72E09D88777F3B969972D21FD4F60E1271579D32AB9A1012C9DB0D3204AAF0347925AC89B256B6AA5073CDE02584E6026ADED3696366D43E5543362749864CDB22E69A18B0124A609BE9D1A0F93C3603").unwrap();
        let m =
            Bytes32::from_hex("AF9B6CAE187C407256FC9D3F3BE37010FFAF55D0E687A128F17C7F62EB6884D3")
                .unwrap();
        let e_c = Bytes::from_hex("5A9B5EB1CBA3968976C3CAB10353472E56570BC880BE7D100337D06CF19EEA7A6EF8A5E06EADAD86FFA22E39F368824879533477751B1AC862A0D8FF50F25063B49D108868A69044CED9B3FD3144808891D5EFB95B1B965283CB4D18E2537E74A45468E9FFD7A2B69B249D963C2DB415B26858E34CADF2D8372A73443386089E34F87E2C1A985F4157A4064285B23F6F4C121449E69BF77465B2E739F19A8A43EF029CAFDE52227AE05C20CAD36DAC9BD985C748646681A7C673234768BB812DBA72C2576316D642242FA01AD3A63279FD804F9337DA93DE5638CCCC7FCDADB608D24FCECA7627FF94E5A87EFB1A7B219244E565EB0158DFB3EA0E3F8B48007709A8753BAA1486647AE7E6CCC876E20625850AA87921B5861484EBC82A7BD8A881853FDAACBF6C3ABFF50C2FACE490F9917BCA178B7860A0EB674B611931FA4F76BE5A87A1A58B7119E152BA05DDCCA66C81FF8094AD6F37B50F9909A785B421A7D010760C1EF260EADF83CEB373CE95865FB401F506F7060F96965AA5CC62BA88141AEFF9AD65836653279CB5936E2960A7AB7867B2044205BC102975E99C24ADF22D073BE5615DCFCCB7E76948664200D8F859C59903CF4DE5FA14926EE35494F19A1C5F4A0DDBB16CFADD3F9014D4BA3C258BC5F70FC0270F567CDA5D288904CE119F4DAF5FF5B8BC483347F443482B403D683E5A5754C5E42D5A2E92887EECF60F242CF1C230CF9A0E736B9FDA160B96F12A760305F21197AFA634B06B03F8328EC43452AF91793E22E3D3A5369CDEA647D5EBE2FACC41476384C8FFA598512A18F4142FA1A432DA2D326F4FD748862850377970AF2143285A9D3B81B7641055848DDD509739785751C727C398992B84BA2A0FAD06218B9456689E85DF5E1EB7FEBA7E7C410F17403E9B4A78686E22BA6A9511535A2C924F25E8BBE4885DE0BD7421E13A26BBE397253126AB5714D6AF55A5579F9341D352DF5F3A315854935D7F9D2FF26009B297BBD8FCF5794CBB877AC1AC9F5F03A599E2DF0C1CB6733A2902CAA514AAF5DA73FA8EBCDEB54B83CC81D4689869B4912FF8856C1F0C9B").unwrap();
        let e_k =
            Bytes32::from_hex("F91B8C7477A6005992EE947BB365EBF1CFE15688BD25DEEFAD54F90922B4B84C")
                .unwrap();
        let ml_kem_512 = ParamSet::ml_kem_512();
        let (k, c) = ml_kem_512.encaps(&ek.as_bytes(), &m.as_bytes()).unwrap();
        assert_eq!(e_k, k);
        assert_eq!(e_c, c);
    }

    // One of the official test vectors, to make sure we're all implementing the same FIPS203 ..
    #[test]
    fn test_ml_keygen_512() {
        let z =
            Bytes32::from_hex("1A394111163803FE2E8519C335A6867556338EADAFA22B5FC557430560CCD693")
                .unwrap();
        let d =
            Bytes32::from_hex("1EB4400A01629D517974E2CD85B9DEF59082DE508E6F9C2B0E341E12965955CA")
                .unwrap();
        let ek = Bytes::from_hex("5B318622F73E6FC6CBA5571D0537894AA890426B835640489AA218972180BB2534BCC477C62CC839135934F3B14CD0808A11557D331103B30F9A8C0CB0FA8F0A2A152E802E48E408087510D5114D4D2399A51530616C7E310528308176D0042710BC8344EC3D4CA810A92978BFABB516D81CAB0753CDF325AC2377A1F96EFC73C15F5AA367A1582A13651B0337C7943C1D54637669686BEBBD392511FFFC9E3A68CBEEC0CE2CF59A8D51C4DE288EB4641DF6610C82D09CDDA418ACD83F0DCA2859B27117E87981AAA8EBA47515812DA2C27ADF9C682E373D5AF294BE3104474B8D14173788965ECCD80322B6CA04240E7D150F2CD4B04066C1924039B9E4A9E06C2B55DBA2FDDABB4065CFE7EBC5AE01CD45C76374683CB1820C34A841836391B9D8C2AA22B29E7436CFCAB789B3CE8AE2700351C1165B7B4F72CC53E913E5668AE75170352A0DE68A5E3819443DB4113161A2019C4930C97011F31540B833E9A890503A7EC3F38C0D94BE3C7501C6161F39099E3CAC0139ACC7271B70D1664A36A89FA4D22857C6C15AD4C52D5C26E23B81DCDA9FD7A49980C5818888AB2538AD91F54E691B7558C63FAE433A7FAB51485989F4335E6187B65041401238AA0A5A932356207796AF2C70363034546F4615499245E1228BFF2C76674634A60C9A04E00FB276C6C00A114BF1B2C8961E740A082940CEEAAB464370BBBB3919C7421BC81C732415A711AA935A4C2C02CB5D0BCBB99CE830EDDBAE4C228E4F095E29FBC27EA2B881697A1D309D28C480C3E9691FB63480BC5C6239B6CCAA41CD52A6209038C2C887BC71C1BD514A0FAA21721A2A5B30ACB168227833A8260422C1F4815EC2ADB207389FB1B817D78FC96063434B6728E18469475DB5D712BC403D8231CF9C8926D0A94B6830881FA5678AD04499F40D5CA83479BA85A70B1196C32A68A6B7FFB40EA6FC3FF020768B91B27F653746546C5E256B14069E827C1616FC7647F8B70F8A32DB551CF715BBB315B7B9BC20FF76847CFC4AEAC23DDC1302EC928CFE40447C761143194DA1415D3D8389F61BAB41EB605729123A320BB54B3B3FBCBC787C46F354C7D7D60F8DFE3729375AEF1891C08A79DE237E39E860061D").unwrap();
        let dk = Bytes::from_hex("EA35075D429C8E81ADA6C4BB97D78624C602FB173DFFC78E5C744FF2FAC345B55E04A3B7325A63F58B43EEF168E259910C35A0952A7ADCF43634889C98918A1CE7B62671C1968A02688160C09B7DE9978B77A557D265A40F7C79F3124247FA0C14F8A9F2C2B939470CAD5ACA5E664ADD7B5444643B559C4BF84C524A063DA7E552C9107BB0887BC22C73F76B63ABA60955A06D1CF34491275859D46FEDEA738FEC14D9920458A7C31E657549A6160480784D5C229AD88932B384A0B0A16DD6C8FCEA56B61A50E67442FB9928B4676D92A6717A564D8E939BB7957D750BB70D2B20D8445A8912BC7940489AE01584DC7D7952A686F245D09C547EF40B74C6748B318059F66F2B76C72ECBB3A64167CAE08C27368902B76DFF684A14482AFAE5837CB6838E685987FA4E3FF7CC4A36631CBA1F77915E7C580853E3C6C84859ADC8C2A15CB131E78305F4BCB4A8100AAB206EC97D14862CF5DA4D3AE2066D4C41BBA9187BECBE0809CE6AAC1FE20BCAE87714BE1542BA9053F1802B65C82909594984A186841B2BFCB3AF38100C5E685E7B9B85515C469CA50B1F799229E024B68A4A412A185677444491689957A576F5029C742DB64C3F63614B43AA23C433A1B37811CDC1184E11BB7D9C20E587A862B364EF59651259B26F8375E4510CBB12A475816A364BA72C07566D50A2B4E4503188B7B465080DB88EE663928C894367492E1617CCBF36CF844B9D17072A3178809629B4B9729CF82C9935AA9B1205AEA6631C8EE23CEE832C46E0583FAC43C5BC5131B934527740AB12877D295C72089073E6A2567B655CCB2965FAA3DECA8315815E7B6514F05BACE8A7000B548993734DE3964F2709542496122BE58A7E536CC827839693492AE88E75EA13154AB109E6BD16F4486EE1C3EA61B8FA259283B3BC2BFB589D3206A3C77521AA08C253B4E4CC8B5F87467EEA18EBD48D92604382DB0C0087A3B501066CB55EC9602D2696380A6A5DF33C05C9B01700B61D0FC169DEBC28B3108B096B24D93B8FFE67066286B0E1461265896E33A8089D8B0B81A9781700A983B28022125A4934575220325B318622F73E6FC6CBA5571D0537894AA890426B835640489AA218972180BB2534BCC477C62CC839135934F3B14CD0808A11557D331103B30F9A8C0CB0FA8F0A2A152E802E48E408087510D5114D4D2399A51530616C7E310528308176D0042710BC8344EC3D4CA810A92978BFABB516D81CAB0753CDF325AC2377A1F96EFC73C15F5AA367A1582A13651B0337C7943C1D54637669686BEBBD392511FFFC9E3A68CBEEC0CE2CF59A8D51C4DE288EB4641DF6610C82D09CDDA418ACD83F0DCA2859B27117E87981AAA8EBA47515812DA2C27ADF9C682E373D5AF294BE3104474B8D14173788965ECCD80322B6CA04240E7D150F2CD4B04066C1924039B9E4A9E06C2B55DBA2FDDABB4065CFE7EBC5AE01CD45C76374683CB1820C34A841836391B9D8C2AA22B29E7436CFCAB789B3CE8AE2700351C1165B7B4F72CC53E913E5668AE75170352A0DE68A5E3819443DB4113161A2019C4930C97011F31540B833E9A890503A7EC3F38C0D94BE3C7501C6161F39099E3CAC0139ACC7271B70D1664A36A89FA4D22857C6C15AD4C52D5C26E23B81DCDA9FD7A49980C5818888AB2538AD91F54E691B7558C63FAE433A7FAB51485989F4335E6187B65041401238AA0A5A932356207796AF2C70363034546F4615499245E1228BFF2C76674634A60C9A04E00FB276C6C00A114BF1B2C8961E740A082940CEEAAB464370BBBB3919C7421BC81C732415A711AA935A4C2C02CB5D0BCBB99CE830EDDBAE4C228E4F095E29FBC27EA2B881697A1D309D28C480C3E9691FB63480BC5C6239B6CCAA41CD52A6209038C2C887BC71C1BD514A0FAA21721A2A5B30ACB168227833A8260422C1F4815EC2ADB207389FB1B817D78FC96063434B6728E18469475DB5D712BC403D8231CF9C8926D0A94B6830881FA5678AD04499F40D5CA83479BA85A70B1196C32A68A6B7FFB40EA6FC3FF020768B91B27F653746546C5E256B14069E827C1616FC7647F8B70F8A32DB551CF715BBB315B7B9BC20FF76847CFC4AEAC23DDC1302EC928CFE40447C761143194DA1415D3D8389F61BAB41EB605729123A320BB54B3B3FBCBC787C46F354C7D7D60F8DFE3729375AEF1891C08A79DE237E39E860061D2B87926182B602639ABB65FEBAF116F6A2FCCC167A51A2E2E6F4494C58336A2E1A394111163803FE2E8519C335A6867556338EADAFA22B5FC557430560CCD693").unwrap();

        let ml_kem_512 = ParamSet::ml_kem_512();
        let kp = ml_kem_512.keygen(&d.as_bytes(), &z.as_bytes()).unwrap();
        assert_eq!(ek, kp.public_enc);
        assert_eq!(dk, kp.private_dec);
    }
}
