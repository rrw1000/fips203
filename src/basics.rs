// Basic functions required by FIPS-203

use crate::types::{Bits, Bytes, Bytes32, IntRange2To3};
use anyhow::{Result, anyhow};
use sha3::{
    Digest, Sha3_256, Sha3_512, Shake256,
    digest::{ExtendableOutput, Update},
};

pub const Q: u32 = 3329;
pub const Z: u32 = 17;

pub fn accumulate_vec(a: &mut [u32; 256], b: &[u32; 256]) {
    for idx in 0..256 {
        a[idx] = (a[idx] + b[idx]) % Q
    }
}

pub fn subtract_vec(a: &mut [u32; 256], b: &[u32; 256]) {
    for idx in 0..256 {
        a[idx] = sub_mod(a[idx], b[idx], Q)
    }
}

/// Z^(2*BitRev7(i) + 1) mod Q
/// == Z^BitRev7(i) * Z^BitRev7(i) * Z
/// log_2 Q = 12, so size of the result is max 24 + 5 = 29 bits and we can afford to reduce only once.
pub fn zeta_2(i: u32) -> u32 {
    let z = zeta_mod(i);
    (z * z * Z) % Q
}

/// The very specific operation x = Z^Bitrev(i) mod Q.
/// Observe that Z^Bitrev(i) mod Q = (Z^(2^0) mod Q + Z^(2^1) mod Q ..) mod Q and that Q^2 < u32::max()
/// We could type the table in the spec in, but that would involve a lot of tedious typing..
/// @todo Compute zeta_mod in build.rs
pub fn zeta_mod(i: u32) -> u32 {
    // A table of Z^(2^i) - you could get away with a much smaller table by
    // recomputing, but this is just as easy.
    // We index backwards to achieve the bitrev() part of our spec.
    const Z_2_I_MOD_Q: [u32; 7] = [
        1729, // 17^64 mod Q
        2580, // 17^32 mod Q
        2642, // 17^16 mod Q
        1062, // 17^8 mod Q
        296,  // 17^4 mod Q
        289,  // 17^2 mod Q
        17,   // 17^1 mod Q
    ];
    // This could be quite efficiently vectorised using dot product and bitrev(), but since
    // we don't have these here ..
    let mut sum: u32 = 1;
    for (bit, val) in Z_2_I_MOD_Q.iter().enumerate() {
        if (i >> bit) & 1 != 0 {
            sum = (sum * val) % Q;
        }
    }
    sum % Q
}

/// Subtract and modulus
/// The comment in alg 8 suggests that mod is defined wrapped here - ie. 2-3 mod 5 = 4, not 1, because 2-3+4 = 5.
/// Luckily the two definitions are identical for +ve numbers - it's only when you have -ves turning up that
/// the difference is perceptible and you need this function to get the answer that (I think) the spec is looking for.
/// Computes a-b mod n
pub fn sub_mod(a: u32, b: u32, n: u32) -> u32 {
    if a < b {
        n - (a.abs_diff(b) % n)
    } else {
        a.abs_diff(b) % n
    }
}

/// S4.1 PRF
pub fn prf(n: IntRange2To3, s: &Bytes32, b: u8) -> Bytes {
    let mut hasher = Shake256::default();
    hasher.update(s.as_bytes());
    hasher.update(&[b]);
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

/// @todo: make this more efficient - probably need a custom poly type.
pub fn compress_poly(x: [u32; 256], d: u8) -> [u32; 256] {
    let mut result: [u32; 256] = [0; 256];
    for i in 0..256 {
        result[i] = compress(x[i], d);
    }
    result
}

/// @todo: make this more efficient - probably need a custom poly type.
pub fn decompress_poly(x: [u32; 256], d: u8) -> [u32; 256] {
    let mut result: [u32; 256] = [0; 256];
    for i in 0..256 {
        result[i] = decompress(x[i], d);
    }
    result
}

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

// S4 4.2.1
// @todo this could be speeded up substantially by not copying the data twice.
pub fn byte_encode(values: &[u32], d: u8) -> Result<Bytes> {
    let mut some_bits = Vec::<u8>::new();
    let m = if d < 12 { 1 << d } else { Q };
    if values.len() != 256 {
        return Err(anyhow!("Require 256 values for byte_encode"));
    }
    values.iter().for_each(|x| {
        let mut x_tmp = *x % m;
        for _ in 0..d {
            some_bits.push((x_tmp & 1).try_into().unwrap());
            x_tmp >>= 1;
        }
    });
    bits_to_bytes(&Bits::from(some_bits))
}

// S4 4.2.1
pub fn byte_decode(bits: &Bytes, d: u8) -> Result<[u32; 256]> {
    let mut result: [u32; 256] = [0; 256];
    let bits = bytes_to_bits(bits)?;
    if bits.len() != usize::from(d) * 256 {
        return Err(anyhow!(
            "Expecting {} bits, got {}",
            usize::from(d) * 256,
            bits.len()
        ));
    }
    let bit_slice = bits.as_slice();
    let mut bit_idx = 0;
    for elem in result.iter_mut() {
        let mut bit = 1;
        for _ in 0..d {
            if bit_slice[bit_idx] != 0 {
                *elem |= bit;
            }
            bit_idx += 1;
            bit <<= 1;
        }
    }
    // We don't protect against values higher than Q, because it's an extra
    // op and we're "supposed" never to get them.
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_byte_decode() {
        let mut test_bytes: [u32; 256] = [14; 256];
        test_bytes[0] = 596;
        test_bytes[1] = 18345;
        let mut test_bytes_mod_6 = test_bytes;
        test_bytes_mod_6[0] = test_bytes[0] % (1 << 6);
        test_bytes_mod_6[1] = test_bytes[1] % (1 << 6);

        let mut test_bytes_mod_12 = test_bytes;
        test_bytes_mod_12[0] = test_bytes[0] % Q;
        test_bytes_mod_12[1] = test_bytes[1] % Q;
        let encoded_6 = byte_encode(&test_bytes, 6).unwrap();
        let encoded_12 = byte_encode(&test_bytes, 12).unwrap();
        let decoded_6 = byte_decode(&encoded_6, 6).unwrap();
        let decoded_12 = byte_decode(&encoded_12, 12).unwrap();
        assert_eq!(decoded_6, test_bytes_mod_6);
        assert_eq!(decoded_12, test_bytes_mod_12);
    }

    #[test]
    fn test_byte_encode() {
        {
            let test_bytes: [u32; 256] = [0; 256];
            let mut expected_result = Vec::<u8>::new();
            for _ in 0..352 {
                expected_result.push(0);
            }
            let enc_1 = byte_encode(&test_bytes, 11).unwrap();
            assert_eq!(Bytes::from(expected_result), enc_1)
        }
        {
            let mut test_bytes: [u32; 256] = [0; 256];
            test_bytes[0] = 5963;
            // 5963 mod 2048 = 1867. In hex, 11 bits 0x74B, backwards 0x4B7
            let mut expected_result = Vec::<u8>::new();
            expected_result.push(0x4B);
            expected_result.push(0x07);
            for _ in 0..350 {
                expected_result.push(0);
            }
            let enc_1 = byte_encode(&test_bytes, 11).unwrap();
            assert_eq!(Bytes::from(expected_result), enc_1)
        }
    }

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
        let result_2 = prf(IntRange2To3::Two, &test32, 0x78);
        let cmp_2 = Bytes::from_hex(
            "291f9389b72dffdab2ee42d4d28f960d3532b4ae78d3a7324fa428ae5896ec46f1c1f1967c7a14d67d06501b7167fa81fc8a9\
             ccaaaa2c93a2a973648e46dc363125ab579d15a9dc0b86f9f1eb47da73b4298b5e06785db53ae607b4a66f9c68464832e0db5\
             cca8ed2146f7b6694f0f6390a6429b2468549c5ee1036e119eb0e3"
        ).unwrap();
        assert_eq!(cmp_2, result_2);

        let result_3 = prf(IntRange2To3::Three, &test32, 0x40);
        let cmp_3 = Bytes::from_hex("dd2b2ffec4eb0e21ab994c3c6a4f45e10975688927fdbeaac224ce3e94b4719c18451e29ecabd25\
                                     6a6fdc7d2731f595f67105d35bfdc9cf4ebc9a21ac5d4054812256b3d73587a9df1dabd14ceeef7\
                                     cc45a3e052e89839e3e6fc352c01852da6972e3de3bd871db12121bac870bb0219c82a6828c9dfc\
                                     577a274771d795e2a134f4b2c12253a035496cd06070eee9ae927360a53b228312b5914febd10206\
                                     4fce6d279cf18db2fb4cb63d82138bb0ea69d8e9da513cbc98689a877a195b323b1").unwrap();
        assert_eq!(cmp_3, result_3);
    }

    #[test]
    fn test_sub_mod() {
        assert_eq!(sub_mod(0, 0, 5), 0);
        assert_eq!(sub_mod(3, 2, 5), 1);
        assert_eq!(sub_mod(2, 3, 5), 4);
        assert_eq!(sub_mod(297, 1486, 5), 1);
    }

    #[test]
    fn test_zeta_mod() {
        // These values computed by python..
        // Z^1 mod Q
        assert_eq!(zeta_mod(2), 2580);
        // Z^219 mod Q (0xdb, which is a binary palindrome)
        assert_eq!(zeta_mod(3), 3289);
        // Z^0x7a , bitrev = 0x5e
        assert_eq!(zeta_mod(6), 1897);
    }

    #[test]
    fn test_zeta_2() {
        // 17^3 % Q
        assert_eq!(zeta_2(1), Q - 17);
        assert_eq!(zeta_2(8), 1637);
        assert_eq!(zeta_2(26), 2156);
    }

    #[test]
    fn test_subtract_vec() {
        let mut a: [u32; 256] = [12; 256];
        let mut b: [u32; 256] = [1; 256];
        b[3] = 47;
        b[255] = 2700;
        subtract_vec(&mut a, &b);
        subtract_vec(&mut a, &b);
        let expected_a: [u32; 256] = [
            10, 10, 10, 3247, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 1270,
        ];
        assert_eq!(expected_a, a);
    }

    #[test]
    fn test_accumulate_vec() {
        let mut a: [u32; 256] = [12; 256];
        let mut b: [u32; 256] = [1; 256];
        b[3] = 47;
        b[255] = 2700;
        accumulate_vec(&mut a, &b);
        accumulate_vec(&mut a, &b);
        let expected_a: [u32; 256] = [
            14, 14, 14, 106, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
            14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
            14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
            14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
            14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
            14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
            14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
            14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
            14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
            14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
            14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
            14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 2083,
        ];
        assert_eq!(expected_a, a);
    }
}
