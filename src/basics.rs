// Basic functions required by FIPS-203

use crate::types::{Bytes, Bytes32, IntRange2To3};
use sha3::{
    Shake256,
    digest::{ExtendableOutput, Update},
};

/// S4.1 PRF
pub fn prf(n: IntRange2To3, s: &Bytes32, b: &Bytes) -> Bytes {
    let mut hasher = Shake256::default();
    hasher.update(s.as_bytes());
    hasher.update(b.as_bytes());
    // 8 * 64 * n bits = 64*2 or 64*3
    match n {
        IntRange2To3::Two => {
            let mut output: [u8; 128] = [0; 128];
            hasher.finalize_xof_into(&mut output);
            Bytes::from_bytes(&output)
        }
        IntRange2To3::Three => {
            let mut output: [u8; 192] = [0; 192];
            hasher.finalize_xof_into(&mut output);
            Bytes::from_bytes(&output)
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prf() {
        let test32 =
            Bytes32::from_hex("b721ec76c0a4fc36969438de2908446f4e189f47f286411b109e8cda786d07d1")
                .unwrap();
        let testrest = Bytes::from_hex("227bf570a667ba06").unwrap();
        let result_2 = prf(IntRange2To3::Two, &test32, &testrest);
        let cmp_2 = Bytes::from_hex(
            "fc4a955618bf042cbfd736ba057053bdd1e5ec9f4394f92df5a38d30b4739ffd595c234a32f523fee38f77eb4\
                                     7f905ca308f921c567346a65a4642076c89444813b4b97aee2e44206a9da37be3431b1a076c6da2127bbab08c\
                                     0ab40171cde70f4e5352d2653a08bc94cce2e34860d156f8411be7431c0c81fb9d2fb3b55ccf1a",
        ).unwrap();
        assert_eq!(cmp_2, result_2);

        let result_3 = prf(
            IntRange2To3::Three,
            &test32,
            &Bytes::from_hex("b5fa31cdaebbf967bf").unwrap(),
        );
        let cmp_3 = Bytes::from_hex("fdbf1b1c2e74aed711cbb3e4f1c4728b846250a6e5c5b2dfb935526d5991e602f9d0b3447fbecd51c7f2199592\
                                     464dad15319885e157d6446366ba4d4ff5fc2ab3aa9e3a126ffc5b256df35edcc16e228b9a718bf6ec3fa7317f4\
                                     153db9911dee230158cf50119e9d7f262f3572fdb6e0913c45fb63bd2176a12f849b977759bc2d09768199126db\
                                     a5f5f29f8b9bd6df2b76ea3d2474eb2f4fe7806091301eb93622cd32493a2d450409262010bf35799306b959695\
                                     d59b83d6d01022c568db1").unwrap();
        assert_eq!(cmp_3, result_3);
    }
}
