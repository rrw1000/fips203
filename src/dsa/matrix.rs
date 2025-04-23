// DSA matrix ops
// @todo should turn these into generic matrix ops at some point.

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

    pub fn as_vec_mut(&mut self) -> &mut Vec<[i32; 256]> {
        &mut self.values
    }

    pub fn at(&self, idx: usize) -> &[i32; 256] {
        &self.values[idx]
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
}
