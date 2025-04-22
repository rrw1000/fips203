// Conversions between data types

use crate::{
    dsa::dsa,
    format,
    types::{Bits, Bytes, IntRange2Or4},
};
use anyhow::Result;

// @todo this could be optimised quite heavily.
pub fn integer_to_bits(x: u32, alpha: u32) -> Bits {
    let mut rv = Bits::default();
    let vec = rv.as_vec_mut();
    let mut rest = x;
    for _ in 0..alpha {
        vec.push((rest % 2) as u8);
        rest >>= 1;
    }
    rv
}

// @todo this could be optimised quite heavily
pub fn bits_to_integer(y: &Bits, alpha: u32) -> u32 {
    let mut result: u32 = 0;
    let y_slice = y.as_slice();
    for (idx, bit) in y_slice.iter().take(alpha as usize).enumerate() {
        result |= (*bit as u32) << idx;
    }
    result
}

pub fn integer_to_bytes(x: u32, alpha: u32) -> Bytes {
    let mut result = Bytes::default();
    let writer = result.as_vec_mut();
    let mut rest = x;
    for _ in 0..alpha {
        writer.push((rest & 0xff) as u8);
        rest >>= 8;
    }
    result
}

pub fn coeff_from_three_bytes(b0: u8, b1: u8, b2: u8) -> Option<u32> {
    let result: u32 = ((b2 & 0x7f) as u32) << 16 | ((b1 as u32) << 8) | (b0 as u32);
    if result < dsa::Q { Some(result) } else { None }
}

pub fn coeff_from_half_byte(b: u8, n: IntRange2Or4) -> Option<i32> {
    match n {
        IntRange2Or4::Two => {
            if b == 15 {
                None
            } else {
                Some(2 - ((b as i32) % 5))
            }
        }
        IntRange2Or4::Four => {
            if b < 9 {
                Some(4 - (b as i32))
            } else {
                None
            }
        }
    }
}

pub fn bitlen(b: u32) -> u32 {
    if b == 0 { 0 } else { b.ilog2() + 1 }
}

pub fn simple_bit_pack(b: u32, r: &[u32; 256]) -> Result<Bytes> {
    let mut z = Bits::default();
    let nr_bits = bitlen(b);
    for v in r.iter() {
        let mut this_coeff = integer_to_bits(*v, nr_bits);
        z.as_vec_mut().append(this_coeff.as_vec_mut());
    }
    format::bits_to_bytes(&z)
}

// Here, since the modulus isn't specified, we have vectors over i32.
pub fn bit_pack(w: &[i32; 256], a: u32, b: u32) -> Result<Bytes> {
    let mut z = Bits::default();
    let nr_bits = bitlen(a + b);
    for v in w.iter() {
        let mut this_coeff = integer_to_bits(((b as i32) - v) as u32, nr_bits);
        z.as_vec_mut().append(this_coeff.as_vec_mut());
    }
    format::bits_to_bytes(&z)
}

pub fn simple_bit_unpack(b: u32, v: &[u8]) -> Result<[u32; 256]> {
    let mut result: [u32; 256] = [0; 256];
    let nr_bits = bitlen(b);
    let z = format::bytes_to_bits(v)?;
    for (idx, item) in result.iter_mut().enumerate() {
        let bit_idx = idx * (nr_bits as usize);
        *item = bits_to_integer(&z.interval(bit_idx..bit_idx + (nr_bits as usize)), nr_bits);
    }
    Ok(result)
}

pub fn bit_unpack(v: &[u8], a: u32, b: u32) -> Result<[i32; 256]> {
    let mut result: [i32; 256] = [0; 256];
    let c = bitlen(a + b);
    let z = format::bytes_to_bits(v)?;
    for (idx, elem) in result.iter_mut().enumerate() {
        let bit_idx = idx * (c as usize);
        *elem =
            (b as i32) - (bits_to_integer(&z.interval(bit_idx..bit_idx + (c as usize)), c) as i32);
    }
    Ok(result)
}

pub fn hint_bit_pack(v: &Vec<[i32; 256]>, w: u32) -> Bytes {
    let mut result = Bytes::default();
    let vec = result.as_vec_mut();
    let k = v.len();
    vec.resize((w as usize) + k, 0);
    let mut index = 0;
    for (i, h) in v.iter().enumerate() {
        for (j, elem) in h.iter().enumerate() {
            if *elem != 0 {
                vec[index] = j as u8;
                index += 1;
            }
        }
        vec[(w as usize) + i] = index as u8;
    }
    result
}

pub fn hint_bit_unpack(y: &[u8], w: u32) -> Option<Vec<[i32; 256]>> {
    // First, recover k.
    let k = y.len() - (w as usize);
    let mut result: Vec<[i32; 256]> = vec![[0; 256]; k];
    let mut index = 0;
    for (i, v) in result.iter_mut().enumerate() {
        if y[(w as usize) + i] < index || y[(w as usize) + i] > (w as u8) {
            return None;
        }
        let first = index;
        while index < y[(w as usize) + i] {
            if index > first && y[(index as usize) - 1] >= y[index as usize] {
                return None;
            }
            v[y[index as usize] as usize] = 1;
            index += 1;
        }
    }
    for i in (index as usize)..((w - 1) as usize) {
        if y[i] != 0 {
            return None;
        }
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integer_to_bits() {
        // LSB first, remember.
        let exp_bit_rep = Bits::from_bitstring("01001100").unwrap();
        let bit_rep = integer_to_bits(0x132, 8);
        assert_eq!(exp_bit_rep, bit_rep);
    }

    #[test]
    fn test_bits_to_integer() {
        let bit_rep = Bits::from_bitstring("1010110011101").unwrap();
        // Expecting 100110101 = 0x135
        assert_eq!(0x135, bits_to_integer(&bit_rep, 9));
        assert_eq!(0x35, bits_to_integer(&bit_rep, 7));

        assert_eq!(0x49, bits_to_integer(&integer_to_bits(0x49, 8), 8));
    }

    #[test]
    fn test_integer_to_bytes() {
        let exp_byte_rep = Bytes::from_hex("bc48").unwrap();
        assert_eq!(exp_byte_rep, integer_to_bytes(0xa148bc, 2));
        assert_eq!(
            Bytes::from_hex("bc480000").unwrap(),
            integer_to_bytes(0x48bc, 4)
        );
    }

    #[test]
    fn test_simple_bit_pack_unpack() {
        let mut r: [u32; 256] = [12; 256];
        r[0] = 45;
        let result = simple_bit_pack(46, &r).unwrap();
        println!("result = {}", hex::encode(result.as_vec()));
        // regression test.
        let expected = Bytes::from_hex("2dc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc330").unwrap();
        assert_eq!(expected, result);
        let unpacked = simple_bit_unpack(46, &result.as_bytes()).unwrap();
        assert_eq!(r, unpacked);
    }

    #[test]
    fn test_bit_pack() {
        let mut r: [i32; 256] = [4; 256];
        r[0] = -1;
        r[1] = 56;
        r[2] = -19;
        let packed = bit_pack(&r, 20, 58).unwrap();
        let unpacked = bit_unpack(&packed.as_bytes(), 20, 58).unwrap();
        assert_eq!(r, unpacked)
    }

    #[test]
    fn test_coeff_from_three_bytes() {
        assert_eq!(Some(0x12345), coeff_from_three_bytes(0x45, 0x23, 0x01));
        // Because the top bit gets sliced off.
        assert_eq!(Some(0x12345), coeff_from_three_bytes(0x45, 0x23, 0x81));
        // > Q.
        assert_eq!(None, coeff_from_three_bytes(0xff, 0xff, 0x7f));
    }

    #[test]
    fn test_hint_bits() {
        let mut v: Vec<[i32; 256]> = vec![[0; 256]; 12];
        v[0][1] = 1;
        v[2][23] = 1;
        v[7][14] = 1;
        v[8][2] = 1;
        let packed = hint_bit_pack(&v, 16);
        let unpacked = hint_bit_unpack(&packed.as_bytes(), 16).unwrap();
        assert_eq!(v, unpacked);
    }
}
