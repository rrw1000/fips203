// Basic functions required by FIPS-203

use crate::types::{Bytes, IntRange2To3, XOF};
use crate::{basics, kem};
use anyhow::{Result, anyhow};

// S4.2.2
pub fn sample_ntt(b: &Bytes) -> Result<[u32; 256]> {
    let mut rv: [u32; 256] = [0; 256];
    if b.len() != 34 {
        return Err(anyhow!(
            "Expected input array length to be 34, but it is {}",
            b.len()
        ));
    }
    let mut ctx = XOF::init();
    ctx.absorb_bytes(b);
    let mut reader = ctx.finalize_xof();
    let mut j = 0;
    while j < 256 {
        let some_bytes = XOF::squeeze(&mut reader, 3)?;
        let c = some_bytes.as_bytes();
        let c0: u32 = c[0].into();
        let c1: u32 = c[1].into();
        let c2: u32 = c[2].into();
        let d1: u32 = c0 + (256 * (c1 % 16));
        let d2: u32 = (c1 / 16) + (16 * c2);
        if d1 < kem::Q {
            rv[j] = d1;
            j += 1;
        }
        if d2 < kem::Q && j < 256 {
            rv[j] = d2;
            j += 1;
        }
    }
    Ok(rv)
}

// S4.2.2
pub fn sample_poly_cbd(b: &Bytes, n: IntRange2To3) -> Result<[u32; 256]> {
    let mut f: [u32; 256] = [0; 256];
    let n_value: usize = n.value().try_into()?;
    let expected_len = n_value * 64;
    if b.len() != expected_len {
        return Err(anyhow!(
            "Expected input to be of length {}, it is {}",
            expected_len,
            b.len()
        ));
    }
    let bit_value = basics::bytes_to_bits(b)?;
    let bits = bit_value.as_slice();
    // Quite deliberately written out in full since there are limited (read: no) test vectors
    for (i, f_r) in f.iter_mut().enumerate() {
        let idx: usize = 2 * i * n_value;
        let mut x: u32 = 0;
        let mut y: u32 = 0;
        for j in 0..n_value {
            x += u32::from(bits[idx + j]);
            y += u32::from(bits[idx + n_value + j]);
        }
        *f_r = basics::sub_mod(x, y, kem::Q);
    }
    Ok(f)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sample_ntt() {
        let test_bytes = Bytes::from_hex(
            "f7192044175f4164755d9ee23f23a9b2\
             b925b9932d4c95758a62150a06b0503e\
             0102",
        )
        .unwrap();
        // This is a pure regression test, since I can't find test vectors for these.
        let expect_sampled_regress: [u32; 256] = [
            25, 2995, 1181, 1691, 3170, 2588, 223, 1599, 1325, 2713, 1823, 2715, 602, 169, 491,
            2202, 2019, 478, 844, 264, 1014, 772, 2302, 2933, 2883, 3279, 252, 1597, 3247, 1646,
            2865, 1235, 830, 3130, 592, 2192, 2778, 2516, 1744, 1333, 441, 3145, 809, 453, 2467,
            1382, 2795, 633, 1005, 841, 125, 2282, 143, 495, 1857, 645, 990, 2442, 404, 1820, 3169,
            3066, 2748, 2891, 2363, 147, 654, 2898, 3178, 2180, 1372, 372, 3260, 965, 1417, 2207,
            893, 2818, 642, 1232, 2248, 478, 961, 1437, 1040, 2282, 2199, 1133, 1068, 1943, 404,
            38, 1252, 2871, 667, 2113, 2300, 1246, 1979, 2276, 1909, 2138, 1042, 1080, 723, 1179,
            1497, 2209, 3240, 101, 1684, 1974, 2666, 1512, 1373, 2699, 407, 2429, 428, 2505, 3163,
            3119, 3243, 2246, 3064, 3164, 1402, 633, 617, 341, 398, 1111, 1993, 2076, 276, 1737,
            781, 2083, 830, 2418, 1633, 1441, 1693, 1007, 843, 2257, 726, 92, 1664, 3169, 2364,
            2703, 2395, 89, 382, 982, 2135, 2817, 2089, 57, 2806, 2605, 1224, 1813, 221, 783, 2138,
            3116, 2264, 1593, 2556, 2778, 551, 2092, 1911, 1200, 2476, 1193, 802, 2843, 1404, 2512,
            3180, 156, 728, 3081, 1133, 870, 2283, 1494, 942, 1593, 3295, 1157, 2887, 2146, 1464,
            2080, 1346, 1420, 1486, 2701, 786, 3267, 1920, 893, 262, 350, 3129, 2117, 1569, 418,
            3097, 1376, 1100, 548, 1044, 2868, 1937, 1261, 2013, 1185, 712, 2220, 247, 839, 956,
            1483, 805, 3224, 1387, 1002, 1456, 275, 3021, 1059, 2160, 1048, 3109, 2221, 983, 727,
            1356, 652, 2853, 86, 1995, 961, 2439, 2003, 1663, 2167, 1842, 2775, 1657, 1632,
        ];
        let sampled = sample_ntt(&test_bytes).unwrap();
        assert_eq!(expect_sampled_regress, sampled);
    }

    #[test]
    fn test_sample_poly_cbd() {
        let test_bytes = Bytes::from_hex(
            "f7192044175f4164755d9ee23f23a9b2\
             b925b9932d4c95758a62150a06b0503e\
             051bdd9706cc3b32214bf3eaf9415216\
             1b8dc1fa12a95084193b9f8bd6ae5e98\
             a2ca9e7bcdc716d6ee2eabc75704e462\
             87a17650d44f0e521b33317148aa073a\
             96feda2128edf6d8c943ee7fd95a36c5\
             994344b8a8d9d31a2ac0184d8a5531ff",
        )
        .unwrap();
        // This is a pure regression test, since I can't find test vectors for these.

        let expect_sampled_regress = [
            1, 0, 0, 1, 0, 1, 3328, 3328, 1, 1, 0, 0, 1, 3328, 3328, 0, 0, 1, 3328, 0, 3328, 0, 1,
            3328, 0, 2, 2, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 2, 0, 3328, 1, 3327, 3328, 0, 0, 0, 1,
            0, 3328, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 3328, 2, 0, 0, 1, 1, 3328, 3328, 1, 0, 0,
            0, 3327, 3327, 1, 2, 1, 2, 1, 1, 1, 3328, 2, 0, 0, 3328, 0, 0, 1, 3328, 1, 0, 0, 1, 1,
            1, 3328, 3328, 1, 3327, 0, 0, 1, 1, 0, 0, 0, 0, 3328, 3328, 0, 1, 1, 2, 0, 0, 1, 3328,
            0, 3328, 3328, 0, 3328, 0, 3328, 0, 1, 0, 0, 3327, 3328, 0, 1, 1, 3328, 3327, 1, 3327,
            0, 1, 0, 3328, 3328, 3328, 3328, 1, 1, 0, 1, 3327, 1, 0, 3328, 0, 3328, 3328, 1, 0, 1,
            3328, 1, 0, 0, 1, 0, 0, 3328, 3328, 0, 3328, 3328, 0, 1, 0, 1, 1, 2, 2, 1, 2, 1, 1,
            3328, 3328, 0, 0, 1, 0, 0, 2, 0, 0, 3328, 0, 0, 3328, 1, 1, 3328, 1, 3328, 3328, 0, 0,
            3328, 3328, 0, 3327, 2, 3328, 3328, 3328, 0, 1, 0, 3328, 0, 0, 0, 2, 0, 3327, 0, 0, 2,
            3328, 3328, 3328, 3328, 1, 3328, 0, 0, 3328, 2, 3328, 0, 1, 0, 1, 0, 3327, 3328, 1,
            3328, 3328, 0, 3328, 0, 0, 1, 2, 0, 0,
        ];

        let sampled = sample_poly_cbd(&test_bytes, IntRange2To3::Two).unwrap();
        assert_eq!(expect_sampled_regress, sampled);
    }
}
