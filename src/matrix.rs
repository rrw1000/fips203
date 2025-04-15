use crate::{basics, mul, ntt};
use anyhow::{Result, anyhow};

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Vector {
    values: [[u32; 256]; 4],
    dimension: u32,
}

impl Vector {
    pub fn new(dimension: u32) -> Self {
        Self {
            values: [[0_u32; 256]; 4],
            dimension,
        }
    }

    pub fn set(&mut self, idx: u32, v: &[u32; 256]) {
        self.values[idx as usize] = *v;
    }

    pub fn ntt(&self) -> Result<Vector> {
        let mut new_vec = Vector::new(self.dimension);
        for i in 0..self.dimension {
            new_vec.values[i as usize] = ntt::ntt(&self.values[i as usize])?;
        }
        Ok(new_vec)
    }

    pub fn inv_ntt(&self) -> Result<Vector> {
        let mut new_vec = Vector::new(self.dimension);
        for i in 0..self.dimension {
            new_vec.values[i as usize] = ntt::inv_ntt(&self.values[i as usize])?;
        }
        Ok(new_vec)
    }

    pub fn is_empty(&self) -> bool {
        self.dimension == 0
    }

    pub fn len(&self) -> u32 {
        self.dimension
    }

    pub fn at(&self, idx: u32) -> [u32; 256] {
        self.values[idx as usize]
    }

    pub fn accumulate(&mut self, v: &Vector) -> Result<()> {
        if v.dimension != self.dimension {
            return Err(anyhow!(
                "Dimension mismatch {} != {}",
                v.dimension,
                self.dimension
            ));
        }
        for j in 0..self.dimension {
            basics::accumulate_vec(&mut self.values[j as usize], &v.values[j as usize]);
        }
        Ok(())
    }

    pub fn add(&self, v_hat: &Vector) -> Result<Vector> {
        if v_hat.dimension != self.dimension {
            return Err(anyhow!(
                "Expected v_hat dimension {}, got {}",
                self.dimension,
                v_hat.dimension
            ));
        }
        let mut result = Vector::new(self.dimension);
        let k = self.dimension;
        for j in 0..k {
            result.values[j as usize] = self.values[j as usize];
            basics::accumulate_vec(&mut result.values[j as usize], &v_hat.values[j as usize]);
        }
        Ok(result)
    }

    pub fn compose_transpose(&self, v_hat: &Vector) -> Result<[u32; 256]> {
        if v_hat.dimension != self.dimension {
            return Err(anyhow!(
                "Expected v_hat dimension {}, got {}",
                self.dimension,
                v_hat.dimension
            ));
        }
        let k = self.dimension;
        let mut z = [0; 256];
        for j in 0..k {
            basics::accumulate_vec(
                &mut z,
                &mul::multiply_ntts(self.values[j as usize], v_hat.values[j as usize])?,
            );
        }
        Ok(z)
    }
}

// Because vectors are slow(er), we always store the maximum size matrix (4x4), and
// just use the upper elements for smaller values of k.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct SquareMatrix {
    values: [[u32; 256]; 16],
    dimension: u32,
}

impl SquareMatrix {
    pub fn new(dimension: u32) -> Self {
        Self {
            values: [[0; 256]; 16],
            dimension,
        }
    }

    pub fn set(&mut self, i: u32, j: u32, v: &[u32; 256]) {
        self.values[self.idx(i, j) as usize] = *v
    }

    fn idx(&self, i: u32, j: u32) -> u32 {
        (j * self.dimension) + i
    }

    pub fn assign_elem(&mut self, x: u32, y: u32, v: [u32; 256]) -> Result<()> {
        let elem = self.idx(x, y);
        self.values[elem as usize] = v;
        Ok(())
    }

    // Compose this matrix with a 256-entry hatted vector via S2.12
    pub fn compose_hat(&self, u: &Vector) -> Result<Vector> {
        let k = self.dimension;
        if u.len() != k {
            return Err(anyhow!(
                "Attempt to multiply a matrix of dimension {} by a vector of dimension {}",
                self.dimension,
                u.len()
            ));
        }
        let mut result = Vector::new(self.dimension);
        for i in 0..k {
            for j in 0..k {
                basics::accumulate_vec(
                    &mut result.values[i as usize],
                    &mul::multiply_ntts(
                        self.values[self.idx(i, j) as usize],
                        u.values[j as usize],
                    )?,
                );
            }
        }
        Ok(result)
    }

    pub fn compose_transpose_hat(&self, u: &Vector) -> Result<Vector> {
        let k = self.dimension;
        if u.len() != k {
            return Err(anyhow!(
                "Attempt to multiply a matrix of dimension {} by a vector of dimension {}",
                self.dimension,
                u.len()
            ));
        }
        let mut result = Vector::new(self.dimension);
        for i in 0..k {
            for j in 0..k {
                basics::accumulate_vec(
                    &mut result.values[i as usize],
                    &mul::multiply_ntts(
                        self.values[self.idx(j, i) as usize],
                        u.values[j as usize],
                    )?,
                );
            }
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    // use super::*;
    // At some point, we'll implement regression tests here. For now, it's easier to debug with the test vectors.
}
