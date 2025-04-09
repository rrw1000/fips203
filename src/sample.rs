// Basic functions required by FIPS-203

use crate::basics;
use crate::types::{Bytes, IntRange2To3, XOF};
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
        let d1: u32 = c0 + 256 * (c1 % 16);
        let d2: u32 = (c1 / 16) + 16 * c2;
        if d1 < basics::Q {
            rv[j] = d1;
            j += 1;
        }
        if d2 < basics::Q && j == 256 {
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
        *f_r = basics::sub_mod(x, y, basics::Q);
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
        let expect_sampled_regress = [
            25, 1181, 3170, 223, 1325, 1823, 2715, 169, 844, 1014, 772, 2933, 3279, 1597, 2865,
            830, 3130, 2192, 2516, 441, 809, 1382, 1005, 125, 143, 1857, 990, 404, 3169, 2891, 654,
            3178, 1372, 372, 965, 2207, 2818, 1232, 478, 961, 1437, 2282, 2199, 1068, 38, 2871,
            2113, 1246, 2276, 1080, 1179, 1497, 3240, 1684, 1974, 1512, 428, 3119, 2246, 3164, 633,
            341, 1993, 1737, 781, 2083, 2418, 1633, 1441, 1007, 2257, 92, 2364, 2395, 382, 2135,
            2089, 2806, 1224, 221, 3116, 1593, 2778, 2092, 2476, 802, 1404, 3180, 728, 1133, 870,
            1494, 942, 3295, 2887, 2146, 1346, 2701, 3267, 262, 2117, 418, 1376, 1044, 2868, 1261,
            1185, 247, 956, 3224, 1002, 275, 3021, 2160, 3109, 2221, 1356, 2853, 1995, 2439, 1663,
            1842, 1657, 500, 2764, 2304, 1034, 73, 97, 2391, 2717, 1901, 2813, 1268, 1831, 2738,
            906, 2626, 1430, 2826, 536, 1937, 614, 3211, 1899, 555, 3122, 190, 466, 787, 2233,
            2433, 2688, 1361, 1160, 588, 643, 1900, 929, 929, 2096, 2907, 3039, 1036, 1957, 2980,
            2080, 301, 2969, 2206, 2364, 883, 1545, 2560, 99, 2133, 1762, 1087, 2425, 675, 2020,
            639, 1290, 738, 986, 873, 2888, 880, 2343, 1861, 2113, 965, 1774, 1123, 2195, 1629,
            2025, 849, 2325, 3189, 43, 1300, 1106, 1756, 1657, 3219, 819, 1911, 1691, 1679, 2417,
            12, 2881, 1434, 2836, 221, 1831, 2373, 1917, 90, 1705, 343, 1759, 1372, 1210, 1059,
            987, 1335, 106, 3081, 118, 1172, 1482, 830, 1997, 2712, 2075, 1322, 1395, 2052, 2627,
            3234, 2909, 2852, 2722, 1223, 644, 1550, 2213, 205, 3042, 3283, 793, 137, 1361, 919,
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
