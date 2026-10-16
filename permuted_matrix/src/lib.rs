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
//! applied, [`PermutedMatrix::row`], [`PermutedMatrix::row_iter`] and
//! [`PermutedMatrix::diagonal`] iterate over elements in logical order, and
//! the `Display` implementation prints the logical matrix.
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

use std::fmt::{self};
use std::iter::Iterator;

use faer::Mat;
use std::ops::{Index, IndexMut};

/// View of a borrowed matrix with permuted rows and columns.
///
/// All positions in the methods of this type are logical positions, see the
/// [crate documentation](crate). The `Debug` output shows the internal state,
/// the physical matrix and both permutations, while `Display` shows the
/// logical matrix.
#[derive(Debug)]
pub struct PermutedMatrix<'a, T> {
    /// Borrowed matrix, its memory layout is never changed by the view
    matrix: &'a mut Mat<T>,
    /// Physical row of each logical row
    row_permutation: Box<[usize]>,
    /// Physical column of each logical column
    column_permutation: Box<[usize]>,
}

impl<'a, T> PermutedMatrix<'a, T> {
    /// Creates a view of `matrix` without any permutation, logical and
    /// physical positions are the same.
    ///
    /// The matrix stays borrowed while the view exists. Writes through the
    /// view change `matrix`, the data is never copied.
    pub fn new(matrix: &'a mut Mat<T>) -> Self {
        let nrows = matrix.nrows();
        let ncols = matrix.ncols();
        PermutedMatrix {
            matrix,
            row_permutation: (0..nrows).collect(),
            column_permutation: (0..ncols).collect(),
        }
    }

    /// Returns the number of rows of the matrix.
    pub fn nrows(&self) -> usize {
        self.matrix.nrows()
    }

    /// Returns the number of columns of the matrix.
    pub fn ncols(&self) -> usize {
        self.matrix.ncols()
    }

    /// Returns the row permutation.
    ///
    /// Entry `i` is the physical row of logical row `i`.
    pub fn row_permutation(&self) -> &[usize] {
        &self.row_permutation
    }

    /// Returns the column permutation.
    ///
    /// Entry `j` is the physical column of logical column `j`. Use it to map
    /// results computed in logical order back to the original columns, for
    /// example the entries of a solution vector.
    pub fn column_permutation(&self) -> &[usize] {
        &self.column_permutation
    }

    /// Returns an iterator over the logical diagonal, the elements at
    /// `(k, k)` for `k` in `0..min(nrows, ncols)`.
    ///
    /// This is not the physical diagonal of the borrowed matrix, unless no
    /// rows or columns were swapped.
    pub fn diagonal(&self) -> impl DoubleEndedIterator<Item = &T> + ExactSizeIterator<Item = &T> {
        let ndiag = self.matrix.ncols().min(self.matrix.nrows());
        (0..ndiag).map(move |i| &self[(i, i)])
    }

    /// Returns an iterator over the elements of logical row `i`, from the
    /// first to the last logical column.
    ///
    /// # Panics
    ///
    /// Panics while iterating if `i` is not smaller than [`Self::nrows`].
    pub fn row(
        &self,
        i: usize,
    ) -> impl DoubleEndedIterator<Item = &T> + ExactSizeIterator<Item = &T> {
        (0..self.matrix.ncols()).map(move |j| &self[(i, j)])
    }

    /// Returns an iterator over all logical rows, each one an iterator over
    /// its elements as returned by [`Self::row`].
    pub fn row_iter(
        &self,
    ) -> impl DoubleEndedIterator<Item = impl DoubleEndedIterator<Item = &T>>
    + ExactSizeIterator<Item = impl ExactSizeIterator<Item = &T>> {
        (0..self.matrix.nrows()).map(move |i| self.row(i))
    }

    /// Returns a new matrix with the permutations applied, its element
    /// `(i, j)` is the logical element `(i, j)` of the view.
    ///
    /// The borrowed matrix is not changed.
    pub fn to_mat(&self) -> Mat<T>
    where
        T: Clone,
    {
        Mat::from_fn(self.matrix.nrows(), self.matrix.ncols(), |i, j| {
            self[(i, j)].clone()
        })
    }

    /// Swaps logical rows `a` and `b`.
    ///
    /// Only the row permutation changes, the data of the matrix doesn't move.
    ///
    /// # Panics
    ///
    /// Panics if `a` or `b` is not smaller than [`Self::nrows`].
    pub fn swap_rows(&mut self, a: usize, b: usize) {
        self.row_permutation.swap(a, b);
    }

    /// Swaps logical columns `a` and `b`.
    ///
    /// Only the column permutation changes, the data of the matrix doesn't
    /// move.
    ///
    /// # Panics
    ///
    /// Panics if `a` or `b` is not smaller than [`Self::ncols`].
    pub fn swap_columns(&mut self, a: usize, b: usize) {
        self.column_permutation.swap(a, b);
    }

    /// Write the logical (permuted) matrix to `f`, formatting each element
    /// with `fmt_elem`.
    ///
    /// Each row is written on its own line as `[a, b, c]`, without a trailing
    /// newline after the last row. Used by the `Display` implementation, which
    /// passes `fmt::Display::fmt` as `fmt_elem`.
    ///
    /// `fmt_elem` receives the same formatter `f`, so formatting options of
    /// the caller, like width and precision in `{:8.3}`, apply to every
    /// element.
    fn fmt_with(
        &self,
        f: &mut fmt::Formatter<'_>,
        fmt_elem: impl Fn(&T, &mut fmt::Formatter<'_>) -> fmt::Result,
    ) -> fmt::Result {
        for (i, row) in self.row_iter().enumerate() {
            if i > 0 {
                writeln!(f)?;
            }
            write!(f, "[")?;
            for (j, val) in row.enumerate() {
                if j > 0 {
                    write!(f, ", ")?;
                }
                fmt_elem(val, f)?;
            }
            write!(f, "]")?;
        }

        Ok(())
    }
}

/// Read access to the logical element `(i, j)`.
///
/// # Panics
///
/// Panics if `i` or `j` is out of bounds.
impl<'a, T> Index<(usize, usize)> for PermutedMatrix<'a, T> {
    type Output = T;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        &self.matrix[(self.row_permutation[i], self.column_permutation[j])]
    }
}

/// Write access to the logical element `(i, j)`, the write goes to the
/// physical element `(row_permutation[i], column_permutation[j])`.
///
/// # Panics
///
/// Panics if `i` or `j` is out of bounds.
impl<'a, T> IndexMut<(usize, usize)> for PermutedMatrix<'a, T> {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut Self::Output {
        &mut self.matrix[(self.row_permutation[i], self.column_permutation[j])]
    }
}

/// Prints the logical matrix, one row per line, like
///
/// ```text
/// [6, 5, 4]
/// [3, 2, 1]
/// ```
///
/// Formatting options apply to every element, `{:.2}` prints all elements
/// with two decimals.
impl<'a, T: fmt::Display> fmt::Display for PermutedMatrix<'a, T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.fmt_with(f, fmt::Display::fmt)
    }
}

#[cfg(test)]
mod unit_tests;
