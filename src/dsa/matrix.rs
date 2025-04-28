// DSA matrix ops
// @todo should turn these into generic matrix ops at some point.

use crate::dsa::{basics, ntt};
use anyhow::{Result, anyhow};
use std::{fmt, iter::zip};

#[derive(Default, Clone, Eq, PartialEq)]
pub struct Vector {
    values: Vec<[i32; 256]>,
}

impl Vector {
    pub fn zero(l: usize) -> Self {
        Self {
            values: vec![[0_i32; 256]; l],
        }
    }

    pub fn as_slice(&self) -> &[[i32; 256]] {
        self.values.as_slice()
    }

    pub fn as_vec(&self) -> &Vec<[i32; 256]> {
        &self.values
    }

    pub fn as_vec_mut(&mut self) -> &mut Vec<[i32; 256]> {
        &mut self.values
    }

    pub fn at(&self, idx: usize) -> &[i32; 256] {
        &self.values[idx]
    }

    pub fn at_mut(&mut self, idx: usize) -> &mut [i32; 256] {
        &mut self.values[idx]
    }

    pub fn len(&self) -> usize {
        self.values.len()
    }
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    pub fn add(&self, other: &Vector) -> Result<Vector> {
        if self.len() != other.len() {
            return Err(anyhow!(
                "Attempt to add vectors of length {} and {}",
                self.len(),
                other.len()
            ));
        }
        let values = zip(self.values.iter(), other.values.iter())
            .map(|(a, b)| basics::poly_add(a, b))
            .collect();
        Ok(Self { values })
    }

    pub fn sub(&self, other: &Vector) -> Result<Vector> {
        if self.len() != other.len() {
            return Err(anyhow!(
                "Attempt to subtract vectors of unequal length {} vs {}",
                self.len(),
                other.len()
            ));
        }
        let values = zip(self.values.iter(), other.values.iter())
            .map(|(a, b)| basics::poly_sub(a, b))
            .collect();
        Ok(Self { values })
    }

    pub fn mod_pm(&self) -> Vector {
        let values = self.values.iter().map(basics::poly_pm).collect();
        Self { values }
    }

    pub fn minus(&self) -> Vector {
        let values = self.values.iter().map(basics::minus).collect();
        Self { values }
    }

    pub fn add_ntt(&self, other: &Vector) -> Vector {
        let values = zip(self.values.iter(), other.values.iter())
            .map(|(a, b)| basics::poly_add(a, b))
            .collect();
        Self { values }
    }

    pub fn mul_vec(&self, c: &[i32; 256]) -> Vector {
        let values = self
            .values
            .iter()
            .map(|vi| basics::poly_mulntt(c, vi))
            .collect();
        Self { values }
    }

    pub fn ntt(&self) -> Vector {
        let values = self.values.iter().map(|vi| ntt::ntt(vi)).collect();
        Self { values }
    }

    pub fn inv_ntt(&self) -> Vector {
        let values = self.values.iter().map(|vi| ntt::inv_ntt(vi)).collect();
        Self { values }
    }

    pub fn power2_round(&self, d: i32) -> (Vector, Vector) {
        let (v1, v2) = self
            .values
            .iter()
            .map(|vi| basics::poly_power2_round(vi, d))
            .unzip();
        (Self { values: v1 }, Self { values: v2 })
    }

    pub fn high_bits(&self, gamma_2: i32) -> Vector {
        let values = self
            .values
            .iter()
            .map(|vi| basics::high_bits_vec(vi, gamma_2))
            .collect();
        Self { values }
    }

    pub fn make_hint(&self, r: &Vector, gamma_2: i32) -> Result<Vector> {
        if self.len() != r.len() {
            return Err(anyhow!(
                "Attempt to make_hint on vectors of unequal length {}, {}",
                self.len(),
                r.len()
            ));
        }
        let values = zip(self.values.iter(), r.values.iter())
            .map(|(a, b)| basics::poly_make_hint(a, b, gamma_2))
            .collect();
        Ok(Self { values })
    }

    pub fn low_bits(&self, gamma_2: i32) -> Vector {
        let values = self
            .values
            .iter()
            .map(|vi| basics::low_bits_vec(vi, gamma_2))
            .collect();
        Self { values }
    }

    pub fn metric_inf(&self) -> i32 {
        self.values
            .iter()
            .map(|v1| v1.iter().map(|x| i32::abs(*x)).max().unwrap_or(0))
            .max()
            .unwrap_or(0)
    }

    pub fn count_ones(&self) -> usize {
        self.values.iter().map(basics::poly_count_ones).sum()
    }

    pub fn reduce(&self) -> Vector {
        let values = self
            .values
            .iter()
            .map(|vi| basics::mod_vec(vi, basics::Q))
            .collect();
        Self { values }
    }
}

impl From<Vec<[i32; 256]>> for Vector {
    fn from(values: Vec<[i32; 256]>) -> Vector {
        Self { values }
    }
}

impl fmt::Debug for Vector {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for k in 0..self.values.len() {
            write!(f, "K({}) x N: [ {} ", k, self.values[k][0])?;
            for v in self.values[k].iter().skip(1) {
                write!(f, ", {v}")?;
            }
            writeln!(f, "]")?;
        }
        writeln!(f)
    }
}

#[derive(Clone, Eq, PartialEq)]
pub struct Matrix {
    // matrix entries in row-major format.
    values: Vec<[i32; 256]>,
    cols: usize,
}

// Corresponds with leancrypto's output.
impl fmt::Debug for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for k in 0..self.rows() {
            for l in 0..self.cols {
                write!(f, "K({}) x L({}) x N:", k, l)?;
                let coeff = self.at(k, l);
                for c in coeff.iter() {
                    write!(f, " 0x{:08x} ", c)?;
                }
                writeln!(f)?;
            }
        }
        writeln!(f)
    }
}

impl Matrix {
    /// Create a new k by l matrix.
    pub fn new(k: usize, l: usize) -> Self {
        // The matrix is stored with all the ls together and the array is addressed
        // (cols * k_idx) + l_offset .
        Matrix {
            values: vec![[0; 256]; k * l],
            cols: l,
        }
    }

    pub fn rows(&self) -> usize {
        self.values.len() / self.cols
    }

    pub fn from_vec(values: Vec<[i32; 256]>, cols: usize) -> Self {
        Matrix { values, cols }
    }

    pub fn values(&mut self) -> &mut Vec<[i32; 256]> {
        &mut self.values
    }

    pub fn value(&mut self, k: usize, l: usize) -> &mut [i32; 256] {
        let idx = k * self.cols + l;
        &mut self.values[idx]
    }
    pub fn at(&self, k: usize, l: usize) -> &[i32; 256] {
        &self.values[self.cols * k + l]
    }

    pub fn mul_ntt(&self, v: &Vector) -> Vector {
        let k = self.rows();
        let l = self.cols;
        let mut w_hat = Vector::zero(k);
        for i in 0..k {
            for j in 0..l {
                *w_hat.at_mut(i) =
                    basics::poly_add(w_hat.at(i), &basics::poly_mulntt(self.at(i, j), v.at(j)));
            }
        }
        w_hat
    }
}
