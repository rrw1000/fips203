use crate::{
    dsa::convert,
    types::{Bytes, Bytes32},
};
use anyhow::Result;

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
            &pk.interval(start_idx..start_idx + vec_len),
        )?);
    }
    Ok((p, t1))
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
