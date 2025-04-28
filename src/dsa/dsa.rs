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
}

#[cfg(test)]
mod tests {
    use super::*;
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
