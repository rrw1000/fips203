use crate::{
    basics, matrix, ntt, sample,
    types::{Bytes, Bytes32, IntRange2To3},
};
use anyhow::{Result, anyhow};

/// This could be a tuple; it is only here to ensure that before you use a (private) decryption key
/// you have to type "private_dec", which will hopefully make you think about where you are exposing it.
#[derive(Default, Clone)]
pub struct KeyPair {
    pub public_enc: Bytes,
    pub private_dec: Bytes,
}

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
    fn k_pke_keygen(&self, d: &Bytes32) -> Result<(Bytes, Bytes)> {
        let mut random_bytes = Bytes::new();
        random_bytes.as_vec_mut().extend_from_slice(d.as_bytes());
        // We know it fits because max k == 4
        random_bytes.as_vec_mut().push(self.k as u8);
        let (p, sigma) = basics::g(&random_bytes)?;
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
                &sample::sample_poly_cbd(&basics::prf(self.n1, &sigma, big_n), self.n1)?,
            );
            big_n += 1;
        }
        let mut e = matrix::Vector::new(self.k);
        for i in 0..self.k {
            e.set(
                i,
                &sample::sample_poly_cbd(&basics::prf(self.n1, &sigma, big_n), self.n1)?,
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

    fn k_pke_encrypt(&self, ek: &Bytes, m: &Bytes32, r: &Bytes32) -> Result<Bytes> {
        let mut big_n = 0;
        let mut ek_pos = 0;
        let ek_slice = ek.as_bytes();
        let mut t_hat = matrix::Vector::new(self.k);
        println!("R= {:?}\n", hex::encode(r.as_bytes()));
        for i in 0..self.k {
            // @todo clean this up - needless copying.
            let b_val = Bytes::from_bytes(&ek_slice[ek_pos..ek_pos + 384]);
            t_hat.set(i, &basics::byte_decode(&b_val, 12)?);
            ek_pos += 384;
        }
        let p = Bytes::from_bytes(&ek_slice[ek_pos..ek_pos + 32]);
        println!(
            "unpacked_key = {t_hat:?} , seed = {:?}",
            hex::encode(p.as_bytes())
        );

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
        println!("A = {a_hat_matrix:?}");
        let mut y = matrix::Vector::new(self.k);
        for i in 0..self.k {
            y.set(
                i,
                &sample::sample_poly_cbd(&basics::prf(self.n1, r, big_n), self.n1)?,
            );
            big_n += 1;
        }
        let mut e1 = matrix::Vector::new(self.k);
        for i in 0..self.k {
            e1.set(
                i,
                &sample::sample_poly_cbd(&basics::prf(self.n2, r, big_n), self.n2)?,
            );
            big_n += 1;
        }
        println!("e1 = {e1:?}");
        let e2 = sample::sample_poly_cbd(&basics::prf(self.n2, r, big_n), self.n2)?;
        println!("y= {y:?}\n e1 = {e1:?}\n e2 = {e2:?}\n");
        let y_hat = y.ntt()?;
        println!("y_hat = {y_hat:?}");
        let mut u = a_hat_matrix.compose_transpose_hat(&y_hat)?;
        println!("noninvb = {u:?}");
        u = u.inv_ntt()?;
        println!("invb = {u:?}");
        u.accumulate(&e1)?;
        let mu = basics::decompress_poly(basics::byte_decode(&Bytes::from(m), 1)?, 1);
        let mut v = t_hat.compose_transpose(&y_hat)?;
        v = ntt::inv_ntt(&v)?;
        basics::accumulate_vec(&mut v, &e2);
        basics::accumulate_vec(&mut v, &mu);
        println!("v_acc = {v:?}");
        println!("final_u = {u:?}");
        let mut c1 = Bytes::new();
        for i in 0..self.k {
            let val = basics::byte_encode(&basics::compress_poly(u.at(i), self.du), self.du)?;
            c1.accumulate(val);
        }
        let c2 = basics::byte_encode(&basics::compress_poly(v, self.dv), self.dv)?;
        c1.accumulate(c2);
        Ok(c1)
    }

    /// Alg 16 ; returns (ek,dk)
    fn keygen_internal(&self, d: &Bytes32, v: &Bytes32) -> Result<(Bytes, Bytes)> {
        let (ek_pke, dk_pke) = self.k_pke_keygen(d)?;
        let mut dk = Bytes::new();
        dk.accumulate(dk_pke);
        dk.accumulate(ek_pke.clone());
        dk.accumulate_32(basics::h(&ek_pke)?);
        dk.accumulate_32(*v);
        Ok((ek_pke, dk))
    }

    /// Alg 17: encaps_internal
    /// Returns (key, ciphertext)
    fn encaps_internal(&self, ek: &Bytes, m: &Bytes32) -> Result<(Bytes32, Bytes)> {
        let mut mek = Bytes::from(m);
        mek.accumulate_32(basics::h(ek)?);
        let (k, r) = basics::g(&mek)?;
        println!("k,r = {k:?}, {r:?}");
        let c = self.k_pke_encrypt(ek, m, &r)?;
        Ok((k, c))
    }

    /// Returns: (key, ciphertext)
    pub fn encaps(&self, ek: &Bytes, m: &Bytes32) -> Result<(Bytes32, Bytes)> {
        // Type check
        let key_len = (384 * self.k) as usize;
        if ek.len() != key_len + 32 {
            return Err(anyhow!(
                "Expected {} EK bytes, got {}",
                key_len + 32,
                ek.len()
            ));
        }
        let ek_slice = ek.as_bytes();
        let mut ptr: usize = 0;
        for i in 0..self.k {
            let the_slice = &ek_slice[ptr..ptr + 384];
            let dec = basics::byte_decode(&Bytes::from_bytes(the_slice), 12)?;
            let enc = basics::byte_encode(&dec, 12)?;
            if the_slice != enc.as_bytes() {
                return Err(anyhow!("Subkey {i} mismatch {the_slice:?} != {enc:?}"));
            }
            ptr += 384;
        }
        self.encaps_internal(ek, m)
    }

    // Somewhat pointless, but the standard does it, so we should too :-)
    pub fn keygen(&self, d: &Bytes32, z: &Bytes32) -> Result<KeyPair> {
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

    // Uses the intermediate value ke
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
        let r = ml_kem_512.keygen(&t.d, &t.z).unwrap();
        assert_eq!(t.ek, r.public_enc);
        assert_eq!(t.dk, r.private_dec);
        let (key, ct) = ml_kem_512.encaps(&t.ek, &t.m).unwrap();
        assert_eq!(t.secret, key);
        assert_eq!(t.ct, ct);
    }

    // #[test]
    // fn test_keygen_512() {
    //     // At last, we can use the official test vectors :-)
    //     let test_vectors = vec![TestSet {
    //         z: Bytes32::from_hex(
    //             "1A394111163803FE2E8519C335A6867556338EADAFA22B5FC557430560CCD693"
    //         ).unwrap(),
    //         d: Bytes32::from_hex(
    //             "1EB4400A01629D517974E2CD85B9DEF59082DE508E6F9C2B0E341E12965955CA",
    //         )
    //         .unwrap(),
    //         ek: Bytes::from_hex(
    //             "5B318622F73E6FC6CBA5571D0537894AA890426B835640489AA218972180BB2534BCC477C62CC839135934F3B14CD0808A11557D331103B30F9A8C0CB0FA8F0A2A152E802E48E408087510D5114D4D2399A51530616C7E310528308176D0042710BC8344EC3D4CA810A92978BFABB516D81CAB0753CDF325AC2377A1F96EFC73C15F5AA367A1582A13651B0337C7943C1D54637669686BEBBD392511FFFC9E3A68CBEEC0CE2CF59A8D51C4DE288EB4641DF6610C82D09CDDA418ACD83F0DCA2859B27117E87981AAA8EBA47515812DA2C27ADF9C682E373D5AF294BE3104474B8D14173788965ECCD80322B6CA04240E7D150F2CD4B04066C1924039B9E4A9E06C2B55DBA2FDDABB4065CFE7EBC5AE01CD45C76374683CB1820C34A841836391B9D8C2AA22B29E7436CFCAB789B3CE8AE2700351C1165B7B4F72CC53E913E5668AE75170352A0DE68A5E3819443DB4113161A2019C4930C97011F31540B833E9A890503A7EC3F38C0D94BE3C7501C6161F39099E3CAC0139ACC7271B70D1664A36A89FA4D22857C6C15AD4C52D5C26E23B81DCDA9FD7A49980C5818888AB2538AD91F54E691B7558C63FAE433A7FAB51485989F4335E6187B65041401238AA0A5A932356207796AF2C70363034546F4615499245E1228BFF2C76674634A60C9A04E00FB276C6C00A114BF1B2C8961E740A082940CEEAAB464370BBBB3919C7421BC81C732415A711AA935A4C2C02CB5D0BCBB99CE830EDDBAE4C228E4F095E29FBC27EA2B881697A1D309D28C480C3E9691FB63480BC5C6239B6CCAA41CD52A6209038C2C887BC71C1BD514A0FAA21721A2A5B30ACB168227833A8260422C1F4815EC2ADB207389FB1B817D78FC96063434B6728E18469475DB5D712BC403D8231CF9C8926D0A94B6830881FA5678AD04499F40D5CA83479BA85A70B1196C32A68A6B7FFB40EA6FC3FF020768B91B27F653746546C5E256B14069E827C1616FC7647F8B70F8A32DB551CF715BBB315B7B9BC20FF76847CFC4AEAC23DDC1302EC928CFE40447C761143194DA1415D3D8389F61BAB41EB605729123A320BB54B3B3FBCBC787C46F354C7D7D60F8DFE3729375AEF1891C08A79DE237E39E860061D",
    //         ).unwrap(),
    //         dk: Bytes::from_hex(
    //             "EA35075D429C8E81ADA6C4BB97D78624C602FB173DFFC78E5C744FF2FAC345B55E04A3B7325A63F58B43EEF168E259910C35A0952A7ADCF43634889C98918A1CE7B62671C1968A02688160C09B7DE9978B77A557D265A40F7C79F3124247FA0C14F8A9F2C2B939470CAD5ACA5E664ADD7B5444643B559C4BF84C524A063DA7E552C9107BB0887BC22C73F76B63ABA60955A06D1CF34491275859D46FEDEA738FEC14D9920458A7C31E657549A6160480784D5C229AD88932B384A0B0A16DD6C8FCEA56B61A50E67442FB9928B4676D92A6717A564D8E939BB7957D750BB70D2B20D8445A8912BC7940489AE01584DC7D7952A686F245D09C547EF40B74C6748B318059F66F2B76C72ECBB3A64167CAE08C27368902B76DFF684A14482AFAE5837CB6838E685987FA4E3FF7CC4A36631CBA1F77915E7C580853E3C6C84859ADC8C2A15CB131E78305F4BCB4A8100AAB206EC97D14862CF5DA4D3AE2066D4C41BBA9187BECBE0809CE6AAC1FE20BCAE87714BE1542BA9053F1802B65C82909594984A186841B2BFCB3AF38100C5E685E7B9B85515C469CA50B1F799229E024B68A4A412A185677444491689957A576F5029C742DB64C3F63614B43AA23C433A1B37811CDC1184E11BB7D9C20E587A862B364EF59651259B26F8375E4510CBB12A475816A364BA72C07566D50A2B4E4503188B7B465080DB88EE663928C894367492E1617CCBF36CF844B9D17072A3178809629B4B9729CF82C9935AA9B1205AEA6631C8EE23CEE832C46E0583FAC43C5BC5131B934527740AB12877D295C72089073E6A2567B655CCB2965FAA3DECA8315815E7B6514F05BACE8A7000B548993734DE3964F2709542496122BE58A7E536CC827839693492AE88E75EA13154AB109E6BD16F4486EE1C3EA61B8FA259283B3BC2BFB589D3206A3C77521AA08C253B4E4CC8B5F87467EEA18EBD48D92604382DB0C0087A3B501066CB55EC9602D2696380A6A5DF33C05C9B01700B61D0FC169DEBC28B3108B096B24D93B8FFE67066286B0E1461265896E33A8089D8B0B81A9781700A983B28022125A4934575220325B318622F73E6FC6CBA5571D0537894AA890426B835640489AA218972180BB2534BCC477C62CC839135934F3B14CD0808A11557D331103B30F9A8C0CB0FA8F0A2A152E802E48E408087510D5114D4D2399A51530616C7E310528308176D0042710BC8344EC3D4CA810A92978BFABB516D81CAB0753CDF325AC2377A1F96EFC73C15F5AA367A1582A13651B0337C7943C1D54637669686BEBBD392511FFFC9E3A68CBEEC0CE2CF59A8D51C4DE288EB4641DF6610C82D09CDDA418ACD83F0DCA2859B27117E87981AAA8EBA47515812DA2C27ADF9C682E373D5AF294BE3104474B8D14173788965ECCD80322B6CA04240E7D150F2CD4B04066C1924039B9E4A9E06C2B55DBA2FDDABB4065CFE7EBC5AE01CD45C76374683CB1820C34A841836391B9D8C2AA22B29E7436CFCAB789B3CE8AE2700351C1165B7B4F72CC53E913E5668AE75170352A0DE68A5E3819443DB4113161A2019C4930C97011F31540B833E9A890503A7EC3F38C0D94BE3C7501C6161F39099E3CAC0139ACC7271B70D1664A36A89FA4D22857C6C15AD4C52D5C26E23B81DCDA9FD7A49980C5818888AB2538AD91F54E691B7558C63FAE433A7FAB51485989F4335E6187B65041401238AA0A5A932356207796AF2C70363034546F4615499245E1228BFF2C76674634A60C9A04E00FB276C6C00A114BF1B2C8961E740A082940CEEAAB464370BBBB3919C7421BC81C732415A711AA935A4C2C02CB5D0BCBB99CE830EDDBAE4C228E4F095E29FBC27EA2B881697A1D309D28C480C3E9691FB63480BC5C6239B6CCAA41CD52A6209038C2C887BC71C1BD514A0FAA21721A2A5B30ACB168227833A8260422C1F4815EC2ADB207389FB1B817D78FC96063434B6728E18469475DB5D712BC403D8231CF9C8926D0A94B6830881FA5678AD04499F40D5CA83479BA85A70B1196C32A68A6B7FFB40EA6FC3FF020768B91B27F653746546C5E256B14069E827C1616FC7647F8B70F8A32DB551CF715BBB315B7B9BC20FF76847CFC4AEAC23DDC1302EC928CFE40447C761143194DA1415D3D8389F61BAB41EB605729123A320BB54B3B3FBCBC787C46F354C7D7D60F8DFE3729375AEF1891C08A79DE237E39E860061D2B87926182B602639ABB65FEBAF116F6A2FCCC167A51A2E2E6F4494C58336A2E1A394111163803FE2E8519C335A6867556338EADAFA22B5FC557430560CCD693",
    //         ).unwrap()
    //     }];
    //     let ml_kem_512 = ParamSet::ml_kem_512();
    //     for t in test_vectors {
    //         let kp = ml_kem_512.keygen(&t.d, &t.z).unwrap();
    //         println!(
    //             "kp e = {:?}\n   d = {:?}\n",
    //             hex::encode(kp.public_enc.as_bytes()),
    //             hex::encode(kp.private_dec.as_bytes())
    //         );
    //         assert_eq!(t.ek, kp.public_enc);
    //         assert_eq!(t.dk, kp.private_dec);
    //     }
    // }
}
