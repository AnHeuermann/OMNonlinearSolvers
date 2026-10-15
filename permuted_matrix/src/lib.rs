/*
 * This file belongs to the OpenModelica Run-Time System
 *
 * Copyright (c) 1998-2026, Open Source Modelica Consortium (OSMC), c/o Linköpings
 * universitet, Department of Computer and Information Science, SE-58183 Linköping, Sweden. All rights
 * reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF THE BSD NEW LICENSE OR THE
 * AGPL VERSION 3 LICENSE OR THE OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.8. ANY
 * USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES RECIPIENT'S
 * ACCEPTANCE OF THE BSD NEW LICENSE OR THE OSMC PUBLIC LICENSE OR THE AGPL
 * VERSION 3, ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the OSMC (Open Source Modelica Consortium) Public License
 * (OSMC-PL) are obtained from OSMC, either from the above address, from the URLs:
 * http://www.openmodelica.org or https://github.com/OpenModelica/ or
 * http://www.ida.liu.se/projects/OpenModelica, and in the OpenModelica distribution. GNU
 * AGPL version 3 is obtained from: https://www.gnu.org/licenses/licenses.html#GPL. The BSD NEW
 * License is obtained from: http://www.opensource.org/licenses/BSD-3-Clause.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY
 * SET FORTH IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF
 * OSMC-PL.
 *
 */

//! A view of a [`faer::Mat`] with permuted rows and columns.
//!
//! [`PermutedMatrix`] swaps rows and columns by swapping entries of two index
//! permutations, the data of the matrix never moves. Elements are accessed by
//! their logical (permuted) position with `m[(i, j)]`, which reads or writes
//! the physical position `(row_permutation[i], column_permutation[j])`.
//! [`PermutedMatrix::to_mat`] creates a new matrix with the permutations
//! applied.
//!
//! This is useful for algorithms like Gaussian elimination with complete
//! pivoting, which swap rows and columns in every step.
//!
//! # Example
//!
//! ```
//! use faer::mat;
//! use permuted_matrix::PermutedMatrix;
//!
//! let mut m = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
//! let mut p = PermutedMatrix::new(&mut m);
//!
//! p.swap_rows(0, 1);
//! p.swap_columns(0, 2);
//! assert_eq!(p[(0, 0)], 6.0);
//! assert_eq!(p.to_mat(), mat![[6.0, 5.0, 4.0], [3.0, 2.0, 1.0]]);
//!
//! // Writes go to the physical position
//! p[(0, 0)] = 10.0;
//! drop(p);
//! assert_eq!(m, mat![[1.0, 2.0, 3.0], [4.0, 5.0, 10.0]]);
//! ```

use faer::Mat;
use std::ops::{Index, IndexMut};

/// View of a borrowed matrix with permuted rows and columns
pub struct PermutedMatrix<'a, T> {
    matrix: &'a mut Mat<T>,
    row_permutation: Box<[usize]>,
    column_permutation: Box<[usize]>,
}

impl<'a, T> PermutedMatrix<'a, T> {
    /// New permutation from matrix without any permutation
    pub fn new(matrix: &'a mut Mat<T>) -> Self {
        let nrows = matrix.nrows();
        let ncols = matrix.ncols();
        PermutedMatrix {
            matrix,
            row_permutation: (0..nrows).collect(),
            column_permutation: (0..ncols).collect(),
        }
    }

    /// Return the number of rows of the matrix
    pub fn nrows(&self) -> usize {
        self.matrix.nrows()
    }

    /// Return the number of columns of the matrix
    pub fn ncols(&self) -> usize {
        self.matrix.ncols()
    }

    /// Get row permutation
    pub fn row_permutation(&self) -> &[usize] {
        &self.row_permutation
    }

    /// Get column permutation
    pub fn column_permutation(&self) -> &[usize] {
        &self.column_permutation
    }

    /// Apply permutation and return permuted matrix
    pub fn to_mat(&self) -> Mat<T>
    where
        T: Clone,
    {
        Mat::from_fn(self.matrix.nrows(), self.matrix.ncols(), |i, j| {
            self[(i, j)].clone()
        })
    }

    /// Switch row a with b
    pub fn swap_rows(&mut self, a: usize, b: usize) {
        self.row_permutation.swap(a, b);
    }

    /// Switch column a with b
    pub fn swap_columns(&mut self, a: usize, b: usize) {
        self.column_permutation.swap(a, b);
    }
}

impl<'a, T> Index<(usize, usize)> for PermutedMatrix<'a, T> {
    type Output = T;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        &self.matrix[(self.row_permutation[i], self.column_permutation[j])]
    }
}

impl<'a, T> IndexMut<(usize, usize)> for PermutedMatrix<'a, T> {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut Self::Output {
        &mut self.matrix[(self.row_permutation[i], self.column_permutation[j])]
    }
}

#[cfg(test)]
mod unit_tests;
