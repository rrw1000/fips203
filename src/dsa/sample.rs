// Sampling

use crate::{
    dsa::dsa,
    format,
    types::{Bits, Bytes, IntRange2Or4},
};
use anyhow::Result;
use sha3::{
    Shake256, Shake256Reader,
    digest::{ExtendableOutput, Update, XofReader},
};

// H is SHAKE256

pub fn sample_in_ball(p: &[u8], tau: u32) -> Result<[i32; 256]> {
    let mut c: [i32; 256] = [0; 256];
    let mut xof = Shake256::default();
    xof.update(p);
    let mut reader = xof.finalize_xof();
    let mut s: [u8; 8] = [0; 8];
    reader.read(&mut s);
    let h_val = format::bytes_to_bits(&s[..])?;
    let h = h_val.as_slice();
    for i in 256 - (tau as usize)..256 {
        let mut j: [u8; 1] = [0];
        reader.read(&mut j);
        while (j[0] as usize) > i {
            reader.read(&mut j);
        }
        c[i] = c[j[0] as usize];
        // -1^1 == -1 , -1 ^^ 2 = 1.
        c[j[0] as usize] = if h[i + (tau as usize) - 256] % 2 == 1 {
            -1
        } else {
            1
        };
    }
    Ok(c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sample_in_ball() {
        // Pure regression ..
        let p = Bytes::from_hex("1234567890").unwrap();
        let in_ball = sample_in_ball(p.as_bytes(), 20).unwrap();
        println!("result {:?}", in_ball);
        let expected = [
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, -1, 0, 1, 0, 0, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0,
        ];
        assert_eq!(expected, in_ball);
    }
}
