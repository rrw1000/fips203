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

#[derive(Clone, Default, Debug, PartialEq, Eq)]
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
    result.n = n;
    result.d = d;
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

#[derive(Clone, Default, Debug, PartialEq, Eq)]
pub struct Sig {
    c: Bytes,
    z: Vec<[i32; 256]>,
    h: Vec<[i32; 256]>,
    gamma_1: u32,
    lambda: u32,
    l: u32,
    w: u32,
}

pub fn sig_encode(sig: &Sig) -> Result<Bytes> {
    let mut rv = Bytes::default();
    if sig.c.len() != ((sig.lambda >> 2) as usize) {
        return Err(anyhow!(
            "c is of length {}, but lambda is {}",
            sig.c.len(),
            sig.lambda
        ));
    }
    rv.append(&sig.c);
    for v in sig.z.iter() {
        rv.accumulate(convert::bit_pack(v, sig.gamma_1 - 1, sig.gamma_1)?);
    }
    rv.accumulate(convert::hint_bit_pack(&sig.h, sig.w));
    Ok(rv)
}

pub fn sig_decode(sig: &[u8], gamma_1: u32, lambda: u32, l: u32, w: u32) -> Result<Sig> {
    let c_len = (lambda >> 2) as usize;
    let c = Bytes::from(&sig[0..c_len]);
    let mut z = Vec::new();
    let mut offset = c_len;
    let z_bits = 32 * (1 + convert::bitlen(gamma_1 - 1)) as usize;
    for _ in 0..l {
        z.push(convert::bit_unpack(
            &sig[offset..offset + z_bits],
            gamma_1 - 1,
            gamma_1,
        )?);
        offset += z_bits;
    }
    let h =
        convert::hint_bit_unpack(&sig[offset..], w).ok_or(anyhow!("Cannot unpack hint bits"))?;
    Ok(Sig {
        c,
        z,
        h,
        gamma_1,
        lambda,
        l,
        w,
    })
}

// @todo test this - though there are no test vectors and no inverse, so let's wait until
// we can do it in the context of ML-DSA.Sign
pub fn w1_encode(w1: &[[u32; 256]], q: u32, gamma_2: u32) -> Result<Bytes> {
    let max_coeff = (q - 1) / ((gamma_2 * 2) - 1);
    let mut result = Bytes::default();
    for v in w1.iter() {
        result.accumulate(convert::simple_bit_pack(max_coeff, v)?);
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

    #[test]
    fn test_sk_encoding() {
        let p =
            Bytes32::from_hex("a58866d2974b947eaee704c6b4b623c0847255df5738603d36f6ab4772bb9ccb")
                .unwrap();
        let big_k =
            Bytes32::from_hex("29cc48f65c9980b7be57b3fe8e54f3fb90d76c2d8df1281ce52a9ed617d94525")
                .unwrap();
        let tr = Bytes::from_hex(
            "29cc48f65c9980b7be57b3fe8e54f3fb90d76c2d8df1281ce52a9ed617d94525f1abf1b151cce457765bd958d09d698783422ef375fce67874dddab649e05202",
        ).unwrap();
        // ML-DSA-44
        let n = 2;
        let d = 13;
        let k = 4;
        let l = 4;
        // l vectors with coeffs in -n .. n
        let s1 = vec![[1; 256], [2; 256], [-1; 256], [0; 256]];
        // k vectors with coeffs in -n .. n
        let s2 = vec![[0; 256], [2; 256], [-2; 256], [0; 256]];
        // k vectors with coeffs in -2^13+1 .. 2^13
        // 2^13 = 8192
        let t0 = vec![[4095; 256], [-3128; 256], [0; 256], [456; 256]];
        let sk_in = SK {
            p,
            big_k,
            tr,
            s1,
            s2,
            t0,
            n,
            d,
        };
        let encoded = sk_encode(&sk_in).unwrap();
        let decoded = sk_decode(&encoded, n, d, k, l).unwrap();
        assert_eq!(sk_in, decoded);
    }

    #[test]
    fn test_sig_encoding() {
        // lambda>>2 is 32
        let c = Bytes::from_hex("1faef05caa3a39f4589225fa27573c2643f400d13de1047fd09d37bc1ce5d068")
            .unwrap();
        // l = 4 k = 4
        let z = vec![[0; 256], [1; 256], [2; 256], [3; 256]];
        let h = vec![[0; 256], [0; 256], [0; 256], [0; 256]];
        let sig = Sig {
            c,
            z,
            h,
            gamma_1: (1 << 17),
            lambda: 128,
            l: 4,
            w: 80,
        };
        let encoded = sig_encode(&sig).unwrap();
        let restored = sig_decode(encoded.as_bytes(), 1 << 17, 128, 4, 80).unwrap();
        assert_eq!(sig, restored);
    }

    #[test]
    fn test_w1_encode() {
        // Pure regression test
        let input_data: [[u32; 256]; 4] = [[0; 256], [1; 256], [4; 256], [3; 256]];
        let encoded = w1_encode(&input_data, 8380417, (8380417 - 1) / 88).unwrap();
        let expected = Bytes::from_hex("000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004411004044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110044110c3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300cc3300c").unwrap();
        assert_eq!(expected, encoded);
    }
}
