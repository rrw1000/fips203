#![cfg(test)]
use crate::types::{Bytes, Bytes32};
use anyhow::{Result, anyhow};
use std::fs;

#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct KEMTestVector {
    /// Randomness.
    pub d: Bytes32,
    pub z: Bytes32,
    /// Randomness which determines the encapsulation key
    pub m: Bytes32,
    /// (public) encryption key.
    pub ek: Bytes,
    /// Decryption key.
    pub dk: Bytes,
    /// The secret key itself (derived from m and ek)
    pub secret: Bytes32,
    /// The ciphertext which conveys the secret key to the holder of dk.
    pub ct: Bytes,
    /// Implicit reject ct
    pub implicit_reject_ct: Bytes,
    /// Implict reject value
    pub implicit_reject: Bytes32,
}

pub fn read_test_vectors(file_name: &str) -> Result<Vec<KEMTestVector>> {
    let contents = fs::read_to_string(file_name)?;
    let mut result: Vec<KEMTestVector> = Vec::new();
    // Annoyingly, the reference code test vectors aren't very well structured so we need to do this in a
    // somewhat ad-hoc manner.
    let mut current = KEMTestVector::default();
    let mut first: bool = true;
    for line in contents.split('\n') {
        if line.contains("***") {
            if !first {
                result.push(current);
                current = KEMTestVector::default();
            } else {
                first = false;
            }
        } else {
            let fields: Vec<&str> = line.trim().split(':').collect();
            if line.is_empty() {
                continue;
            }
            let (name, value) = (fields[0], fields[1].trim());
            match name {
                "d" => current.d = Bytes32::from_hex(value)?,
                "z" => current.z = Bytes32::from_hex(value)?,
                "m" => current.m = Bytes32::from_hex(value)?,
                "Public Key" => current.ek = Bytes::from_hex(value)?,
                "Secret Key" => current.dk = Bytes::from_hex(value)?,
                "Shared Secret A" => (),
                "Shared Secret B" => current.secret = Bytes32::from_hex(value)?,
                "Ciphertext" => current.ct = Bytes::from_hex(value)?,
                "BadCT" => current.implicit_reject_ct = Bytes::from_hex(value)?,
                "Pseudorandom shared Secret A" => {
                    current.implicit_reject = Bytes32::from_hex(value)?
                }
                _ => {
                    return Err(anyhow!("Invalid field in test vector file: {name}"));
                }
            }
        }
    }
    Ok(result)
}
