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

use faer::Mat;
use std::ops::{Index, IndexMut};

struct PermutedMatrix<'a, T> {
    matrix: &'a mut Mat<T>,
    row_permutation: Box<[usize]>,
    column_permutation: Box<[usize]>,
}

impl<'a, T> PermutedMatrix<'a, T> {
    /// New permutation from matrix without any permutation
    fn new(matrix: &'a mut Mat<T>) -> Self {
        let nrows = matrix.nrows();
        let ncols = matrix.ncols();
        PermutedMatrix {
            matrix,
            row_permutation: (0..nrows).collect(),
            column_permutation: (0..ncols).collect(),
        }
    }

    /// Return the number of rows of the matrix
    fn nrows(&self) -> usize {
        self.matrix.nrows()
    }

    /// Return the number of rows of the matrix
    fn ncols(&self) -> usize {
        self.matrix.ncols()
    }

    /// Get row permutation
    fn row_permutation(&self) -> &[usize] {
        &self.row_permutation
    }

    /// Get column permutation
    fn column_permutation(&self) -> &[usize] {
        &self.column_permutation
    }

    /// Apply permutation and return permuted matrix
    fn to_mat(&self) -> Mat<T>
    where
        T: Clone,
    {
        Mat::from_fn(self.matrix.nrows(), self.matrix.ncols(), |i, j| {
            self[(i, j)].clone()
        })
    }

    /// Switch row a with b
    fn swap_rows(&mut self, a: usize, b: usize) {
        self.row_permutation.swap(a, b);
    }

    /// Switch column a with b
    fn swap_columns(&mut self, a: usize, b: usize) {
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
mod permuted_matrix_tests {
    use super::*;
    use faer::mat;

    #[test]
    fn identity_permutation_test() {
        let mut m = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];

        let m_permutation = PermutedMatrix::new(&mut m);

        let m = m_permutation.to_mat();
        assert_eq!(m, mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
    }

    #[test]
    fn swap_rows_test() {
        let mut m = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];

        let mut m_permutation = PermutedMatrix::new(&mut m);

        m_permutation.swap_rows(0, 1);
        let m = m_permutation.to_mat();
        assert_eq!(m, mat![[4.0, 5.0, 6.0], [1.0, 2.0, 3.0]]);
    }

    #[test]
    fn swap_columns_test() {
        let mut m = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];

        let mut m_permutation = PermutedMatrix::new(&mut m);

        m_permutation.swap_columns(0, 1);
        let m = m_permutation.to_mat();
        assert_eq!(
            m,
            mat![[2.0, 1.0, 3.0], [5.0, 4.0, 6.0]],
            "Permutation 0 <-> 1 failed!"
        );

        m_permutation.swap_columns(2, 1);
        let m = m_permutation.to_mat();
        assert_eq!(
            m,
            mat![[2.0, 3.0, 1.0], [5.0, 6.0, 4.0]],
            "Permutation 2 <-> 1 failed!"
        );
    }

    #[test]
    fn mat_outliving_perm_test() {
        let mut m = mat![[1.0, 2.0], [3.0, 4.0]];

        let m_permutation = PermutedMatrix::new(&mut m);
        drop(m_permutation);

        assert_eq!(m, mat![[1.0, 2.0], [3.0, 4.0]])
    }

    #[test]
    fn new_test() {
        let mut m = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];

        let m_permutation = PermutedMatrix::new(&mut m);

        assert_eq!(m_permutation.nrows(), 2);
        assert_eq!(m_permutation.ncols(), 3);
        assert_eq!(m_permutation.row_permutation(), &[0, 1]);
        assert_eq!(m_permutation.column_permutation(), &[0, 1, 2]);
    }

    #[test]
    fn swap_rows_repeated_test() {
        let mut m = mat![[1.0], [2.0], [3.0]];

        let mut m_permutation = PermutedMatrix::new(&mut m);

        // Second swap moves a row that was already moved by the first one
        m_permutation.swap_rows(0, 1);
        m_permutation.swap_rows(2, 1);
        assert_eq!(m_permutation.row_permutation(), &[1, 2, 0]);
        assert_eq!(m_permutation.to_mat(), mat![[2.0], [3.0], [1.0]]);
    }

    #[test]
    fn swap_same_index_test() {
        let mut m = mat![[1.0, 2.0], [3.0, 4.0]];

        let mut m_permutation = PermutedMatrix::new(&mut m);

        m_permutation.swap_rows(1, 1);
        m_permutation.swap_columns(0, 0);
        assert_eq!(m_permutation.row_permutation(), &[0, 1]);
        assert_eq!(m_permutation.column_permutation(), &[0, 1]);
    }

    #[test]
    fn swap_rows_and_columns_test() {
        let mut m = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];

        let mut m_permutation = PermutedMatrix::new(&mut m);

        m_permutation.swap_rows(0, 1);
        m_permutation.swap_columns(0, 2);
        assert_eq!(
            m_permutation.to_mat(),
            mat![[6.0, 5.0, 4.0], [3.0, 2.0, 1.0]]
        );
    }

    #[test]
    fn index_test() {
        let mut m = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];

        let mut m_permutation = PermutedMatrix::new(&mut m);

        m_permutation.swap_rows(0, 1);
        m_permutation.swap_columns(1, 2);

        // Logical (0, 0) is physical (1, 0), logical (1, 1) is physical (0, 2)
        assert_eq!(m_permutation[(0, 0)], 4.0);
        assert_eq!(m_permutation[(0, 1)], 6.0);
        assert_eq!(m_permutation[(1, 1)], 3.0);
        assert_eq!(m_permutation[(1, 2)], 2.0);
    }

    #[test]
    fn index_mut_test() {
        let mut m = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];

        let mut m_permutation = PermutedMatrix::new(&mut m);

        m_permutation.swap_rows(0, 1);
        m_permutation.swap_columns(0, 2);
        m_permutation[(0, 0)] = 10.0;
        assert_eq!(m_permutation[(0, 0)], 10.0);
        drop(m_permutation);

        // Writes go to the physical position, the memory layout is unchanged
        assert_eq!(m, mat![[1.0, 2.0, 3.0], [4.0, 5.0, 10.0]]);
    }

    #[test]
    fn to_mat_keeps_original_test() {
        let mut m = mat![[1.0, 2.0], [3.0, 4.0]];

        let mut m_permutation = PermutedMatrix::new(&mut m);

        m_permutation.swap_rows(0, 1);
        let mut copy = m_permutation.to_mat();
        copy[(0, 0)] = 10.0;
        assert_eq!(m_permutation[(0, 0)], 3.0);
        drop(m_permutation);

        assert_eq!(m, mat![[1.0, 2.0], [3.0, 4.0]]);
    }

    #[test]
    #[should_panic]
    fn index_row_out_of_bounds_test() {
        let mut m = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];

        let m_permutation = PermutedMatrix::new(&mut m);

        let _ = m_permutation[(2, 0)];
    }

    #[test]
    #[should_panic]
    fn index_column_out_of_bounds_test() {
        let mut m = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];

        let m_permutation = PermutedMatrix::new(&mut m);

        let _ = m_permutation[(0, 3)];
    }

    #[test]
    #[should_panic]
    fn swap_rows_out_of_bounds_test() {
        let mut m = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];

        let mut m_permutation = PermutedMatrix::new(&mut m);

        m_permutation.swap_rows(0, 2);
    }

    #[test]
    #[should_panic]
    fn swap_columns_out_of_bounds_test() {
        let mut m = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];

        let mut m_permutation = PermutedMatrix::new(&mut m);

        m_permutation.swap_columns(0, 3);
    }
}
