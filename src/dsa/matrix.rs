// DSA matrix ops
// @todo should turn these into generic matrix ops at some point.

use crate::dsa::{basics, ntt};
use anyhow::{Result, anyhow};
use std::iter::zip;

#[derive(Default, Debug, Clone, Eq, PartialEq)]
pub struct Vector {
    values: Vec<[i32; 256]>,
}

impl Vector {
    pub fn zero(l: usize) -> Self {
        Self {
            values: vec![[0_i32; 256]; l],
        }
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
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Matrix {
    // matrix entries in row-major format.
    values: Vec<[i32; 256]>,
    cols: usize,
}

impl Matrix {
    pub fn new(x: usize, y: usize) -> Self {
        Matrix {
            values: vec![[0; 256]; x * y],
            cols: x,
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

    pub fn value(&mut self, x: usize, y: usize) -> &mut [i32; 256] {
        let idx = x + self.cols * y;
        &mut self.values[idx]
    }
    pub fn at(&self, x: usize, y: usize) -> &[i32; 256] {
        &self.values[x + self.cols * y]
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
