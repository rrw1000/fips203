// Basic functions common across the FIPS PQ standards

use crate::types::{Bits, Bytes};
use anyhow::{Result, anyhow};

pub fn accumulate_vec(a: &mut [u32; 256], b: &[u32; 256], n: u32) {
    for idx in 0..256 {
        a[idx] = (a[idx] + b[idx]) % n
    }
}

pub fn subtract_vec(a: &mut [u32; 256], b: &[u32; 256], n: u32) {
    for idx in 0..256 {
        a[idx] = sub_mod(a[idx], b[idx], n)
    }
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
pub fn bytes_to_bits(b: &[u8]) -> Result<Bits> {
    let mut result: Vec<u8> = Vec::new();
    b.iter().for_each(|x| {
        let mut val = *x;
        for _ in 0..8 {
            result.push(val & 1);
            val >>= 1;
        }
    });
    Ok(Bits::from(result))
}

#[cfg(test)]
mod tests {
    use super::*;

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
        assert_eq!(
            expected_output,
            bytes_to_bits(&test_input.as_bytes()).unwrap()
        );
    }

    #[test]
    fn test_sub_mod() {
        assert_eq!(sub_mod(0, 0, 5), 0);
        assert_eq!(sub_mod(3, 2, 5), 1);
        assert_eq!(sub_mod(2, 3, 5), 4);
        assert_eq!(sub_mod(297, 1486, 5), 1);
    }

    #[test]
    fn test_subtract_vec() {
        let mut a: [u32; 256] = [12; 256];
        let mut b: [u32; 256] = [1; 256];
        b[3] = 47;
        b[255] = 2700;
        subtract_vec(&mut a, &b, 3329);
        subtract_vec(&mut a, &b, 3329);
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
        accumulate_vec(&mut a, &b, 3329);
        accumulate_vec(&mut a, &b, 3329);
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
