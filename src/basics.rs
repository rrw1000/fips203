// Basic functions required by FIPS-203

use crate::types::{Bits, Bytes, Bytes32, IntRange2To3};
use anyhow::{Result, anyhow};
use sha3::{
    Digest, Sha3_256, Sha3_512, Shake256,
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

/// S4.1 H
pub fn h(s: &Bytes) -> Result<Bytes32> {
    let mut hasher = Sha3_256::new();
    Update::update(&mut hasher, s.as_bytes());
    let rv = hasher.finalize();
    Bytes32::try_from(rv.as_slice())
}

/// S4.1 J
pub fn j(s: &Bytes) -> Result<Bytes32> {
    let mut hasher = Shake256::default();
    hasher.update(s.as_bytes());
    let mut output = [0u8; 32];
    hasher.finalize_xof_into(&mut output);
    Ok(Bytes32::new(output))
}

/// S4 1.G
pub fn g(s: &Bytes) -> Result<(Bytes32, Bytes32)> {
    let mut hasher = Sha3_512::new();
    Update::update(&mut hasher, s.as_bytes());
    let rv: [u8; 64] = hasher.finalize().as_slice().try_into()?;
    Ok((
        Bytes32::try_from(&rv[0..32])?,
        Bytes32::try_from(&rv[32..64])?,
    ))
}

/// S4 4.2.1
pub fn bits_to_bytes(b: &Bits) -> Result<Bytes> {
    if b.len() % 8 != 0 {
        return Err(anyhow!(
            "Attempt to give bits_to_bytes() a bitstring not an integer multiple of 8 long"
        ));
    }
    let mut result = Vec::new();
    let bit_vec = b.as_vec().as_slice();
    let mut bitstring_idx = 0;
    // Using for() loops here to give the compiler a bit more info.
    for _ in 0..(b.len() / 8) {
        let mut val: u8 = 0;
        for bit_idx in 0..8 {
            val |= bit_vec[bitstring_idx] << bit_idx;
            bitstring_idx += 1;
        }
        result.push(val)
    }
    Ok(Bytes::from(result))
}

/// S4 4.2.1
pub fn bytes_to_bits(b: &Bytes) -> Result<Bits> {
    let mut result: Vec<u8> = Vec::new();
    b.as_vec().iter().for_each(|x| {
        let mut val = *x;
        for _ in 0..8 {
            result.push(val & 1);
            val >>= 1;
        }
    });
    Ok(Bits::from(result))
}

const Q: u32 = 3329;

/// S4 4.2.1 - recall that we're told q = 3329 and bit len(q) == 12
/// This means we can get away with this - there are probably more
/// elegant ways, but..
pub fn compress(x: u32, d: u8) -> u32 {
    // This is fine, because log_2 3329 < 13, d is at most 12, and 13+12 < 31
    // Do an extra shift left, then add 1 to get the rounding, then shift right.
    let num = (((x << (d + 1)) / Q) + 1) >> 1;
    // We hope the compiler is bright enough to turn this into a mask op.
    num % (1 << d)
}

/// S4 4.2.1
pub fn decompress(y: u32, d: u8) -> u32 {
    // Same trick again for rounding
    ((((Q << 1) * y) / (1 << d)) + 1) >> 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compress() {
        // Test data worked out via speedcrunch from the spec.
        assert_eq!(compress(144, 10), 44);
        assert_eq!(compress(2367, 6), 46);
    }

    #[test]
    fn test_decompress() {
        assert_eq!(decompress(1410, 6), 73342);
        assert_eq!(decompress(44, 10), 143);
        assert_eq!(decompress(46, 6), 2393);
    }

    #[test]
    fn test_bits_to_bytes() {
        // 0x7081 - little endian bitstrings
        let test_input = Bits::from_bitstring("1110000000011000").unwrap();
        let expected_output = Bytes::from_hex("0718").unwrap();
        assert_eq!(expected_output, bits_to_bytes(&test_input).unwrap());
    }

    #[test]
    fn test_bytes_to_bits() {
        let expected_output = Bits::from_bitstring("0100111101101001").unwrap();
        let test_input = Bytes::from_hex("f296").unwrap();
        assert_eq!(expected_output, bytes_to_bits(&test_input).unwrap());
    }

    #[test]
    fn test_g() {
        let test_input = Bytes::from_hex("a453b5addf93e22deecd1c263041aa6f").unwrap();
        let expected_output = Bytes::from_hex(
            "1e967a12a311e816dfcfcae978aa908e2d83d03409047e9423440de6491d5048e8692661ed379ab44fd766cfd5185fd835a8d98e3e096e40a801d0d544025195"
        ).unwrap();
        let some_bytes = expected_output.as_bytes();
        let expect_left = Bytes32::try_from(&some_bytes[0..32]).unwrap();
        let expect_right = Bytes32::try_from(&some_bytes[32..64]).unwrap();
        let (left, right) = g(&test_input).unwrap();
        assert_eq!(expect_left, left);
        assert_eq!(expect_right, right);
    }

    #[test]
    fn test_j() {
        let test_input = Bytes::from_hex("3664bd9582a9044ae82b10d5bd368f7a").unwrap();
        let expected_output =
            Bytes32::from_hex("b4c117f00bb0ebfadcdcf9e7aff1cbf0bf44e5c6f0eec5686e29d0c25e72cd85")
                .unwrap();
        assert_eq!(expected_output, j(&test_input).unwrap());
    }

    #[test]
    fn test_h() {
        let test_input = Bytes::from_hex("f20e7d42bbe29148c4d3e9e0556d14d2").unwrap();
        let expected_output =
            Bytes32::from_hex("cb86ad64b294f4ace98d08c736752476e6c033b5e54bee859f8ec168bb3d53d5")
                .unwrap();
        assert_eq!(expected_output, h(&test_input).unwrap());
    }

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
