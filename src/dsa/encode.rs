use crate::{
    dsa::convert,
    types::{Bytes, Bytes32},
};
use anyhow::{Result, anyhow};

// bitlen is bitlen(q-1)-d
pub fn pk_encode(p: &Bytes32, t1: &Vec<[u32; 256]>, bitlen: u32) -> Result<Bytes> {
    let mut pk = Bytes::new();
    let max_coeff = (1 << bitlen) - 1;
    pk.accumulate_32(*p);
    for t in t1.iter() {
        pk.accumulate(convert::simple_bit_pack(max_coeff, t)?);
    }
    Ok(pk)
}

pub fn pk_decode(pk: &Bytes, bitlen: u32) -> Result<(Bytes32, Vec<[u32; 256]>)> {
    // @todo optimise.
    let p = pk.interval(0..32).to_bytes32()?;
    let vec_len = (32 * bitlen) as usize;
    let max_coeff = (1 << bitlen) - 1;
    let mut t1: Vec<[u32; 256]> = Vec::new();
    for start_idx in (32..pk.len()).step_by(vec_len) {
        t1.push(convert::simple_bit_unpack(
            max_coeff,
            &pk.as_bytes()[start_idx..start_idx + vec_len],
        )?);
    }
    Ok((p, t1))
}

#[derive(Clone, Default, Debug)]
pub struct SK {
    p: Bytes32,
    big_k: Bytes32,
    tr: Bytes,
    s1: Vec<[i32; 256]>,
    s2: Vec<[i32; 256]>,
    t0: Vec<[i32; 256]>,
    n: u32,
    d: u32,
}

#[allow(clippy::too_many_arguments)]
// @todo we might want to validate the lengths of these vectors, but
// this is other side of optimisation, when the data structures will
// likely change.
pub fn sk_encode(sk: &SK) -> Result<Bytes> {
    let mut encoded = Bytes::new();
    encoded.append_32(&sk.p);
    encoded.append_32(&sk.big_k);
    encoded.append(&sk.tr);
    for v in sk.s1.iter() {
        encoded.accumulate(convert::bit_pack(v, sk.n, sk.n)?);
    }
    for v in sk.s2.iter() {
        encoded.accumulate(convert::bit_pack(v, sk.n, sk.n)?);
    }
    let d_exp = 1 << (sk.d - 1);
    for v in sk.t0.iter() {
        encoded.accumulate(convert::bit_pack(v, d_exp - 1, d_exp)?)
    }
    Ok(encoded)
}

pub fn sk_decode(encoded: &Bytes, n: u32, d: u32, k: u32, l: u32) -> Result<SK> {
    let mut result = SK::default();
    let byte_slice = encoded.as_bytes();
    // For efficiency, we do this by slice.
    let bl = (32 * convert::bitlen(2 * n)) as usize;
    result.p = Bytes32::try_from(&byte_slice[0..32])?;
    result.big_k = Bytes32::try_from(&byte_slice[32..64])?;
    result.tr = Bytes::from_bytes(&byte_slice[64..128]);
    // Hopefully this will be constant-propagated by the compiler, so it's really
    // just handy readability sugar.
    let mut offset = 128;
    fn push_if_valid(val: &mut Vec<[i32; 256]>, cand: [i32; 256], n: u32) -> Result<()> {
        let nx = n as i32;
        for c in cand {
            if c < -nx || c > nx {
                return Err(anyhow!(
                    "Invalid coefficient in unpacked value - {c} with n = {n}"
                ));
            }
        }
        val.push(cand);
        Ok(())
    }
    for _ in 0..l {
        let cand = convert::bit_unpack(&byte_slice[offset..offset + bl], n, n)?;
        push_if_valid(&mut result.s1, cand, n)?;
        offset += bl;
    }
    for _ in 0..k {
        let cand = convert::bit_unpack(&byte_slice[offset..offset + bl], n, n)?;
        push_if_valid(&mut result.s2, cand, n)?;
        offset += bl;
    }
    let bytes_per_t0 = (d << 5) as usize;
    let d_exp = 1 << (d - 1);
    for _ in 0..k {
        result.t0.push(convert::bit_unpack(
            &byte_slice[offset..offset + bytes_per_t0],
            d_exp - 1,
            d_exp,
        )?);
        offset += bytes_per_t0;
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pk_encoding() {
        let p =
            Bytes32::from_hex("160903758735d4ea86db4e6f431ef68185bdef10c8da5cc64ead6d90f9822137")
                .unwrap();
        let t1: Vec<[u32; 256]> = vec![[2; 256], [8; 256], [6; 256], [31; 256]];
        let coded = pk_encode(&p, &t1, 5).unwrap();
        let (dec_p, dec_t1) = pk_decode(&coded, 5).unwrap();
        assert_eq!(p, dec_p);
        assert_eq!(t1, dec_t1);
    }
}
