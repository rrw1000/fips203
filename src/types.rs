use anyhow::{Result, anyhow};
use hex;
use sha3::{
    Shake128, Shake128Reader,
    digest::{ExtendableOutput, Update},
};
use std::iter::Iterator;
use std::ops::Index;
use std::{fmt, io::Read};

pub struct XOF(Shake128);

impl XOF {
    // Named for the function in the standard.
    pub fn init() -> Self {
        Self(Shake128::default())
    }

    pub fn absorb_bytes(&mut self, b: &Bytes) {
        self.0.update(b.as_bytes());
    }

    pub fn absorb_slice(&mut self, b: &[u8]) {
        self.0.update(b);
    }

    pub fn finalize_xof(&mut self) -> Shake128Reader {
        self.0.clone().finalize_xof()
    }

    pub fn squeeze(rdr: &mut Shake128Reader, l: usize) -> Result<Bytes> {
        let mut rv = vec![0; l];
        rdr.read_exact(rv.as_mut_slice())?;
        Ok(Bytes::from(rv))
    }
}

/// Represented like this because we want variable-sized bitfields to avoid having to specify the length
/// in the type. We may change later..
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Bits(Vec<u8>);

impl Bits {
    pub fn as_slice(&self) -> &[u8] {
        self.0.as_slice()
    }

    pub fn as_vec(&self) -> &Vec<u8> {
        &self.0
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn as_vec_mut(&mut self) -> &mut Vec<u8> {
        &mut self.0
    }

    pub fn from_bitstring(bit_string: &str) -> Result<Self> {
        let mut result = Vec::new();
        // Hard without try_map()
        for c in bit_string.chars() {
            match c {
                '0' => result.push(0),
                '1' => result.push(1),
                _ => {
                    return Err(anyhow!(
                        "Invalid bit string element '{:#02x}'",
                        u32::from(c)
                    ));
                }
            }
        }
        Ok(Self(result))
    }
}

impl From<Vec<u8>> for Bits {
    fn from(bits: Vec<u8>) -> Self {
        Self(bits)
    }
}

#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct Bytes(Vec<u8>);

impl Bytes {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn accumulate_32(&mut self, other: Bytes32) {
        self.0.extend(other.as_bytes())
    }

    pub fn accumulate(&mut self, other: Bytes) {
        self.0.extend(other.as_bytes())
    }

    pub fn as_vec_mut(&mut self) -> &mut Vec<u8> {
        &mut self.0
    }

    pub fn from_bytes(bytes: &[u8]) -> Self {
        Self(Vec::from(bytes))
    }

    pub fn as_bytes(&self) -> &[u8] {
        self.0.as_slice()
    }

    pub fn as_vec(&self) -> &Vec<u8> {
        &self.0
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn from_hex(hex_string: &str) -> Result<Self> {
        let hex_string = hex_string.trim_start_matches("0x");
        let bytes = hex::decode(hex_string)?;
        Ok(Self(bytes))
    }
}

impl From<&Bytes32> for Bytes {
    fn from(bytes: &Bytes32) -> Self {
        Self(bytes.0.to_vec())
    }
}

impl From<Vec<u8>> for Bytes {
    fn from(bytes: Vec<u8>) -> Self {
        Self(bytes)
    }
}

#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
pub struct Bytes32([u8; 32]);

impl Bytes32 {
    pub fn new(bytes: [u8; 32]) -> Self {
        Self(bytes)
    }

    pub fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }

    pub fn into_bytes(self) -> [u8; 32] {
        self.0
    }

    pub fn from_hex(hex_string: &str) -> Result<Self> {
        let hex_string = hex_string.trim_start_matches("0x");

        if hex_string.len() != 64 {
            return Err(anyhow::anyhow!(
                "Invalid hex string length for Bytes32: expected 64 characters"
            ));
        }
        let decoded = hex::decode(hex_string)?;
        let mut bytes = [0u8; 32];
        bytes.copy_from_slice(&decoded);
        Ok(Self(bytes))
    }

    pub fn try_from(bytes: &[u8]) -> Result<Self> {
        if bytes.len() != 32 {
            return Err(anyhow::anyhow!(
                "Invalid byte slice length: expected 32 bytes, got {}",
                bytes.len()
            ));
        }

        let mut array = [0u8; 32];
        array.copy_from_slice(bytes);

        Ok(Self(array))
    }
}

impl fmt::Display for Bytes32 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "0x")?;
        for byte in self.0.iter() {
            write!(f, "{:02x}", byte)?;
        }
        Ok(())
    }
}

impl From<[u8; 32]> for Bytes32 {
    fn from(bytes: [u8; 32]) -> Self {
        Self(bytes)
    }
}

impl From<Bytes32> for [u8; 32] {
    fn from(bytes: Bytes32) -> Self {
        bytes.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntRange2To3 {
    Two = 2,
    Three = 3,
}

impl IntRange2To3 {
    pub fn value(&self) -> u32 {
        match self {
            IntRange2To3::Two => 2,
            IntRange2To3::Three => 3,
        }
    }

    pub fn try_from(value: u32) -> Option<Self> {
        match value {
            2 => Some(IntRange2To3::Two),
            3 => Some(IntRange2To3::Three),
            _ => None,
        }
    }
}

impl fmt::Display for IntRange2To3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.value())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BitVector {
    // Underlying storage as bytes
    data: Vec<u8>,
}

impl BitVector {
    /// Create a new empty BitVector
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    /// Create a BitVector with the specified number of bytes, initialized to 0
    pub fn with_capacity_bytes(bytes: usize) -> Self {
        Self {
            data: vec![0; bytes],
        }
    }

    /// Create a BitVector from an existing byte vector
    pub fn from_bytes(bytes: Vec<u8>) -> Self {
        Self { data: bytes }
    }

    /// Returns the number of bits in the vector (always a multiple of 8)
    pub fn len(&self) -> usize {
        self.data.len() * 8
    }

    /// Returns the number of bytes in the vector
    pub fn byte_len(&self) -> usize {
        self.data.len()
    }

    /// Returns true if the bit vector is empty
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Get the bit at the specified index
    pub fn get(&self, index: usize) -> Option<bool> {
        if index >= self.len() {
            return None;
        }

        let byte_index = index / 8;
        let bit_index = index % 8;
        let mask = 1 << bit_index;

        Some((self.data[byte_index] & mask) != 0)
    }

    /// Set the bit at the specified index
    pub fn set(&mut self, index: usize, value: bool) -> Result<()> {
        if index >= self.len() {
            return Err(anyhow::anyhow!(
                "Index out of bounds: {} >= {}",
                index,
                self.len()
            ));
        }

        let byte_index = index / 8;
        let bit_index = index % 8;
        let mask = 1 << bit_index;

        if value {
            self.data[byte_index] |= mask;
        } else {
            self.data[byte_index] &= !mask;
        }

        Ok(())
    }

    /// Get the underlying bytes
    pub fn as_bytes(&self) -> &[u8] {
        &self.data
    }

    /// Get a mutable reference to the underlying bytes
    pub fn as_bytes_mut(&mut self) -> &mut Vec<u8> {
        &mut self.data
    }

    /// Count the number of bits set to 1
    pub fn count_ones(&self) -> usize {
        self.data
            .iter()
            .map(|&byte| byte.count_ones() as usize)
            .sum()
    }

    /// Creates an iterator over the bits
    pub fn iter(&self) -> BitVectorIterator {
        BitVectorIterator {
            bitvector: self,
            index: 0,
        }
    }
}

impl Default for BitVector {
    fn default() -> Self {
        Self::new()
    }
}

impl Index<usize> for BitVector {
    type Output = bool;

    fn index(&self, index: usize) -> &Self::Output {
        if self.get(index).unwrap_or(false) {
            &true
        } else {
            &false
        }
    }
}

/// Iterator over the bits in a BitVector
pub struct BitVectorIterator<'a> {
    bitvector: &'a BitVector,
    index: usize,
}

impl Iterator for BitVectorIterator<'_> {
    type Item = bool;

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.bitvector.get(self.index);
        if result.is_some() {
            self.index += 1;
        }
        result
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.bitvector.len().saturating_sub(self.index);
        (remaining, Some(remaining))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bytes_from_hex() {
        // Test valid hex
        let result = Bytes::from_hex("deadbeef").unwrap();
        assert_eq!(result.as_bytes(), &[0xde, 0xad, 0xbe, 0xef]);

        // Test with 0x prefix
        let result = Bytes::from_hex("0xdeadbeef").unwrap();
        assert_eq!(result.as_bytes(), &[0xde, 0xad, 0xbe, 0xef]);

        // Test empty string
        let result = Bytes::from_hex("").unwrap();
        assert_eq!(result.as_bytes().len(), 0);
    }

    #[test]
    fn test_bytes32_from_hex() {
        // Test valid 32-byte hex
        let hex_str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
        let result = Bytes32::from_hex(hex_str).unwrap();

        let expected: [u8; 32] = [
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            24, 25, 26, 27, 28, 29, 30, 31,
        ];
        assert_eq!(result.as_bytes(), &expected);

        // Test with 0x prefix
        let with_prefix = format!("0x{}", hex_str);
        let result = Bytes32::from_hex(&with_prefix).unwrap();
        assert_eq!(result.as_bytes(), &expected);

        // Test invalid length
        let result = Bytes32::from_hex("deadbeef");
        assert!(result.is_err());
    }

    #[test]
    fn test_bytes32_try_from() {
        // Test valid 32-byte array
        let bytes: Vec<u8> = (0..32).collect();
        let result = Bytes32::try_from(&bytes).unwrap();

        let expected: [u8; 32] = [
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            24, 25, 26, 27, 28, 29, 30, 31,
        ];
        assert_eq!(result.as_bytes(), &expected);

        // Test invalid length
        let short_bytes = vec![0; 16];
        let result = Bytes32::try_from(&short_bytes);
        assert!(result.is_err());

        let long_bytes = vec![0; 64];
        let result = Bytes32::try_from(&long_bytes);
        assert!(result.is_err());
    }

    #[test]
    fn test_defaults() {
        // Test Bytes default
        let default_bytes = Bytes::default();
        assert_eq!(default_bytes.as_bytes().len(), 0);

        // Test Bytes32 default
        let default_bytes32 = Bytes32::default();
        let all_zeros = [0u8; 32];
        assert_eq!(default_bytes32.as_bytes(), &all_zeros);
    }

    #[test]
    fn test_bit_vector() {
        // Test creation
        let mut bv = BitVector::with_capacity_bytes(3); // 3 bytes = 24 bits
        assert_eq!(bv.len(), 24);
        assert_eq!(bv.byte_len(), 3);
        assert_eq!(bv.as_bytes().len(), 3);

        // Test setting and getting bits
        for i in 0..24 {
            assert_eq!(bv.get(i), Some(false));
            if i % 2 == 0 {
                bv.set(i, true).unwrap();
            }
        }

        for i in 0..24 {
            assert_eq!(bv.get(i), Some(i % 2 == 0));
        }

        // Test count_ones
        assert_eq!(bv.count_ones(), 12); // 12 bits set to 1

        // Test out of bounds
        assert_eq!(bv.get(24), None);
        assert!(bv.set(24, true).is_err());

        // Test from_bytes
        let data = vec![0b10101010, 0b01010101];
        let bv2 = BitVector::from_bytes(data);
        assert_eq!(bv2.len(), 16);
        assert_eq!(bv2.byte_len(), 2);

        // First byte (10101010) should have bits 0,2,4,6 set
        assert_eq!(bv2.get(0), Some(false));
        assert_eq!(bv2.get(1), Some(true));
        assert_eq!(bv2.get(2), Some(false));
        assert_eq!(bv2.get(3), Some(true));
        assert_eq!(bv2.get(4), Some(false));
        assert_eq!(bv2.get(5), Some(true));
        assert_eq!(bv2.get(6), Some(false));
        assert_eq!(bv2.get(7), Some(true));

        // Second byte (01010101) should have bits 8,10,12,14 set
        assert_eq!(bv2.get(8), Some(true));
        assert_eq!(bv2.get(9), Some(false));
        assert_eq!(bv2.get(10), Some(true));
        assert_eq!(bv2.get(11), Some(false));
        assert_eq!(bv2.get(12), Some(true));
        assert_eq!(bv2.get(13), Some(false));
        assert_eq!(bv2.get(14), Some(true));
        assert_eq!(bv2.get(15), Some(false));

        // Test iteration
        let iter_values: Vec<bool> = bv2.iter().collect();
        assert_eq!(iter_values.len(), 16);
        assert_eq!(iter_values[0], false);
        assert_eq!(iter_values[1], true);

        // Test IndexMut
        assert_eq!(bv2[0], false);
        assert_eq!(bv2[1], true);
    }
}
