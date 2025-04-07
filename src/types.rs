use anyhow::Result;
use hex;
use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bytes(Vec<u8>);

impl Bytes {
    pub fn new() -> Self {
        Self(Vec::new())
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

    pub fn from_hex(hex_string: &str) -> Result<Self> {
        let hex_string = hex_string.trim_start_matches("0x");
        let bytes = hex::decode(hex_string)?;
        Ok(Self(bytes))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
}
