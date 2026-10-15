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
use permuted_matrix::PermutedMatrix;

/// Gaussian elimination with complete ("total") pivoting
///
/// Solve a linear equation system
/// ```tex
/// A · y = 0,   with one entry of y fixed to 1
/// ```
/// where `A` is a `n×(n+1)` column-major matrix, `y` a vector of size `n+1`.
fn total_pivot_search(matrix: &Mat<f64>, mut vector_y: &[f64], pos: Option<usize>) {
    let nrows = matrix.nrows();
    let ncols = matrix.ncols();
    assert!(
        nrows + 1 == ncols,
        "Wrong matrix shape {0}x{1}",
        nrows,
        ncols
    );

    assert!(
        pos.unwrap_or(0) <= ncols,
        "Position argument `pos` out of range!"
    );

    // Assume full rank of matrix [dim_1 × dim_2] = [n × (n+1)]
    let rank = nrows;

    let mut row_permutation: Box<[usize]> = (0..nrows).collect();
    let mut column_permutation: Box<[usize]> = (0..ncols).collect();

    let mut pivot_element = nrows;
    match pos {
        None => {
            pivot_element = ncols;
        }
        Some(pos) => {
            column_permutation[nrows] = pos;
            column_permutation[pos] = nrows;
        }
    };

    for i in 0..nrows {}

    // return position_larges_value;
}

#[cfg(test)]
mod total_pivot_tests {
    use super::*;
    use faer::mat;

    #[test]
    #[ignore = "total_pivot_search is not implemented yet"]
    /// Solve:
    /// ```text
    /// 2·x₁ + 1·x₂ = 5
    /// 1·x₁ + 3·x₂ = 10
    /// ```
    ///represented as
    /// ```text
    ///        x₁   x₂   last column
    /// A = [   2    1    -5  ]
    ///     [   1    3   -10  ]
    /// ```
    fn total_pivot_test_1() {
        let a = mat![[2.0, 1.0, -5.0], [1.0, 3.0, -10.0],];

        let vector_y: [f64; 3] = [0.0, 0.0, 0.0];

        total_pivot_search(&a, &vector_y, Some(2));

        assert_eq!(vector_y, [1., 3., 1.], "Test failed");
    }
}

/// Pivot element found by [`pivot_search`]
struct Pivot {
    /// Logical row of the pivot element
    row: usize,
    /// Logical column of the pivot element
    column: usize,
    /// Absolute value of the pivot element
    abs_value: f64,
}

/// Find the pivot element for elimination step `stage`.
///
/// Searches the not yet eliminated part of `mat`, logical rows and columns
/// `stage..`, for the entry with the largest absolute value. Rows are searched
/// one after another, and an entry only replaces the current candidate if it
/// is strictly larger, so on ties the first entry found wins, starting with
/// the diagonal entry `(stage, stage)`. This matches
/// `getIndicesOfPivotElement` of the C runtime.
///
/// Returns the [`Pivot`] with the logical position of the entry in `mat` and
/// its absolute value. An absolute value below `f64::EPSILON` means the
/// remaining part is (numerically) zero and the matrix is rank-deficient.
///
/// # Panics
///
/// Panics if `stage` is not a valid row and column index of `mat`.
fn pivot_search(mat: &PermutedMatrix<f64>, stage: usize) -> Pivot {
    let mut max_abs_value = mat[(stage, stage)].abs();
    let mut col_index = stage;
    let mut row_index = stage;

    for row in stage..mat.nrows() {
        for col in stage..mat.ncols() {
            let current_abs_value = mat[(row, col)].abs();
            if (current_abs_value > max_abs_value) {
                max_abs_value = current_abs_value;
                col_index = col;
                row_index = row;
            }
        }
    }

    return Pivot {
        row: row_index,
        column: col_index,
        abs_value: max_abs_value,
    };
}

#[cfg(test)]
mod pivot_search_tests {
    use super::*;
    use faer::mat;

    /// Run `pivot_search` and return `(row, column, abs_max)`
    fn search(m: &PermutedMatrix<f64>, stage: usize) -> (usize, usize, f64) {
        let pivot = pivot_search(m, stage);
        (pivot.row, pivot.column, pivot.abs_value)
    }

    #[test]
    fn largest_absolute_value_test() {
        let mut a = mat![[1.0, -7.0], [3.0, 2.0]];
        let m = PermutedMatrix::new(&mut a);

        // The largest absolute value is negative
        assert_eq!(search(&m, 0), (0, 1, 7.0));
    }

    #[test]
    fn tie_keeps_first_test() {
        let mut a = mat![[2.0, 5.0], [5.0, 1.0]];
        let m = PermutedMatrix::new(&mut a);

        // Row by row, (0, 1) is found before (1, 0). Strict `>` keeps the first.
        assert_eq!(search(&m, 0), (0, 1, 5.0));
    }

    #[test]
    fn diagonal_wins_tie_test() {
        let mut a = mat![[4.0, -4.0], [4.0, 4.0]];
        let m = PermutedMatrix::new(&mut a);

        assert_eq!(search(&m, 0), (0, 0, 4.0));
    }

    #[test]
    fn zero_matrix_test() {
        let mut a = mat![[0.0, 0.0], [0.0, 0.0]];
        let m = PermutedMatrix::new(&mut a);

        assert_eq!(search(&m, 0), (0, 0, 0.0));
    }

    #[test]
    fn stage_ignores_eliminated_part_test() {
        let mut a = mat![[9.0, 8.0, 7.0], [6.0, 1.0, 2.0], [5.0, 3.0, 1.0]];
        let m = PermutedMatrix::new(&mut a);

        // At stage 1 only rows 1..3 and columns 1..3 are searched
        assert_eq!(search(&m, 1), (2, 1, 3.0));
    }

    #[test]
    fn permuted_matrix_test() {
        let mut a = mat![[1.0, 2.0], [3.0, 4.0]];
        let mut m = PermutedMatrix::new(&mut a);
        m.swap_rows(0, 1);

        // Logical matrix is [[3, 4], [1, 2]], the position is logical too
        assert_eq!(search(&m, 0), (0, 1, 4.0));
    }
}
