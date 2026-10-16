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

use faer::{Mat, diag};
use log::{trace, warn};
use permuted_matrix::PermutedMatrix;

use crate::pivot_search::pivot_search;

/// Gaussian elimination with complete ("total") pivoting
///
/// Solve a linear equation system
/// ```tex
/// A · x = 0,   with one entry of x fixed to 1
/// ```
/// where `A` is a `n×(n+1)` column-major matrix, `x` a vector of size `n+1`.
///
/// From C function `solveSystemWithTotalPivotSearch`.
fn solve(matrix: &mut Mat<f64>, pos: Option<usize>) {
    let n_rows = matrix.nrows(); // n
    let n_cols = matrix.ncols(); // n+1

    assert!(n_rows + 1 == n_cols, "Wrong matrix shape {n_rows}x{n_cols}");
    assert!(
        pos.unwrap_or(0) <= n_cols,
        "Position argument `pos` out of range!"
    );
    assert!(matrix.is_all_finite(), "Matrix has infinite or NaN values!");

    let mut matrix = PermutedMatrix::new(matrix);

    let mut pivot_element = if pos.is_some() { n_rows } else { n_cols };

    if let Some(pos) = pos {
        matrix.swap_columns(n_rows, pos);
    }

    let rank = forward_elimination(&mut matrix);
    if rank < n_rows {
        warn!(target: "NLS", "Matrix (nearly) singular. Continuing anyway.");
        trace!(target: "NLS", "rank = {rank}");
        trace!(target: "NLS", "position = {:?}", pos);
    }
    trace!(target: "NLS", "Linear System Matrix after forward elimination:\n{0}", matrix);

    let determinant_abs = matrix.diagonal().product::<f64>().abs();

    // return position_larges_value;
}

#[cfg(test)]
mod total_pivot_tests {
    use super::*;
    use faer::mat;

    #[test]
    // #[ignore = "total_pivot_search is not implemented yet"]
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
    /// Solution: `x = [1, 3]`, with the fixed last entry `[1, 3, 1]`.
    fn total_pivot_test_1() {
        let mut matrix = mat![[2.0, 1.0, -5.0], [1.0, 3.0, -10.0],];

        solve(&mut matrix, Some(2));
    }

    #[test]
    #[ignore = "total_pivot_search is not implemented yet"]
    /// Solve a system with zeros, including on the diagonal:
    /// ```text
    /// 0·x₁ + 2·x₂ + 0·x₃ = 4
    /// 3·x₁ + 0·x₂ + 1·x₃ = 6
    /// 0·x₁ + 1·x₂ + 4·x₃ = 14
    /// ```
    /// represented as
    /// ```text
    ///        x₁   x₂   x₃   last column
    /// A = [   0    2    0    -4  ]
    ///     [   3    0    1    -6  ]
    ///     [   0    1    4   -14  ]
    /// ```
    /// Solution: `x = [1, 2, 3]`, with the fixed last entry `[1, 2, 3, 1]`.
    ///
    /// Without pivoting the first step would divide by `A[0, 0] = 0`. Total
    /// pivoting starts with the largest entry, `4` at `(2, 2)`.
    fn total_pivot_test_zeros() {
        let mut matrix = mat![
            [0.0, 2.0, 0.0, -4.0],
            [3.0, 0.0, 1.0, -6.0],
            [0.0, 1.0, 4.0, -14.0],
        ];

        solve(&mut matrix, Some(3));
    }

    #[test]
    #[should_panic = "Matrix has infinite or NaN values!"]
    /// A matrix with NaN should fail with a clear error message
    fn nan_test() {
        let mut m = mat![
            [0.0, 1.0, 0.0, -1.0],
            [1.0, 0.0, 0.0, f64::NAN],
            [0.0, 1.0, 4.0, -3.0]
        ];

        solve(&mut m, Some(3));

        print!("{:#?}", m);
    }
}

/// Forward elimination with complete pivoting
///
/// Returning rank of `matrix`.
fn forward_elimination(matrix: &mut PermutedMatrix<f64>) -> usize {
    let n_rows = matrix.nrows(); // n
    let mut rank: usize = n_rows;

    for stage in 0..n_rows {
        let pivot = pivot_search(&matrix, stage);

        if pivot.abs_value < f64::EPSILON {
            rank = stage;
            break;
        }

        // Swap
        matrix.swap_rows(stage, pivot.row);
        matrix.swap_columns(stage, pivot.column);

        // Gauss elimination
        eliminate_below_pivot(&mut *matrix, stage);
    }

    rank
}

/// Single forward elimination step
///
/// Eliminates the entries below the pivot `(i,i)` by subtracting multiples of
/// row `i`.
fn eliminate_below_pivot(matrix: &mut PermutedMatrix<f64>, i: usize) {
    for k in i + 1..matrix.nrows() {
        let multiplier = -matrix[(k, i)] / matrix[(i, i)];

        for j in i + 1..matrix.ncols() {
            matrix[(k, j)] = matrix[(k, j)] + multiplier * matrix[(i, j)];
        }

        matrix[(k, i)] = 0.0;
    }
}

#[cfg(test)]
mod forward_elimination_tests {
    use super::*;
    use equator::assert;
    use faer::mat;
    use faer::utils::approx::{ApproxEq, CwiseMat};

    use permuted_matrix::PermutedMatrix;

    #[test]
    /// Example from Wiki: https://en.wikipedia.org/wiki/Gaussian_elimination
    fn gauss_elimination_test() {
        let mut m = mat![
            [1.0, 3.0, 1.0, -9.0],
            [1.0, 1.0, -1.0, -1.0],
            [3.0, 11.0, 5.0, -35.0]
        ];
        let mut p = PermutedMatrix::new(&mut m);

        // Step 0: subtract 1·row 0 from row 1 and 3·row 0 from row 2
        let expected_result = mat![
            [1.0, 3.0, 1.0, -9.0],
            [0.0, -2.0, -2.0, 8.0],
            [0.0, 2.0, 2.0, -8.0]
        ];
        eliminate_below_pivot(&mut p, 0);
        assert_eq!(p.to_mat(), expected_result, "Wrong result after step 0");

        // Step 1: add row 1 to row 2, which becomes zero (the matrix has rank 2)
        let expected_result = mat![
            [1.0, 3.0, 1.0, -9.0],
            [0.0, -2.0, -2.0, 8.0],
            [0.0, 0.0, 0.0, 0.0]
        ];
        eliminate_below_pivot(&mut p, 1);
        assert_eq!(p.to_mat(), expected_result, "Wrong result after step 1");

        // Step 2: no rows below row 2, nothing changes
        eliminate_below_pivot(&mut p, 2);
        assert_eq!(p.to_mat(), expected_result, "Wrong result after step 2");
    }

    #[test]
    /// Forward elimination of
    /// ```text
    /// 0·x₁ + 2·x₂ + 0·x₃ = 4
    /// 3·x₁ + 0·x₂ + 1·x₃ = 6
    /// 0·x₁ + 1·x₂ + 4·x₃ = 14
    /// ```
    /// with the last column `-b`. All four columns are pivot candidates, as in
    /// the C implementation with `pos = -1`.
    ///
    /// | Stage | Pivot                   | Swaps                     |
    /// |-------|-------------------------|---------------------------|
    /// | 0     | `-14` at `(2, 3)`       | rows 0 ↔ 2, columns 0 ↔ 3 |
    /// | 1     | `3` at logical `(1, 1)` | none                      |
    /// | 2     | `12/7` at `(2, 2)`      | none                      |
    ///
    /// The result is upper triangular in logical order:
    /// ```text
    /// [ -14   0     1     4  ]
    /// [   0   3   -3/7  -5/7 ]
    /// [   0   0   12/7  -8/7 ]
    /// ```
    fn forward_elimination_test() {
        let mut m = mat![
            [0.0, 2.0, 0.0, -4.0],
            [3.0, 0.0, 1.0, -6.0],
            [0.0, 1.0, 4.0, -14.0]
        ];
        let mut p = PermutedMatrix::new(&mut m);

        let rank = forward_elimination(&mut p);

        assert_eq!(rank, 3, "Wrong rank");
        assert_eq!(p.row_permutation(), &[2, 1, 0], "Wrong row permutation");
        assert_eq!(
            p.column_permutation(),
            &[3, 0, 1, 2],
            "Wrong column permutation"
        );

        let expected_result = mat![
            [-14.0, 0.0, 1.0, 4.0],
            [0.0, 3.0, -3.0 / 7.0, -5.0 / 7.0],
            [0.0, 0.0, 12.0 / 7.0, -8.0 / 7.0]
        ];
        let approx_eq = CwiseMat(ApproxEq::eps());
        assert!(p.to_mat() ~ expected_result);
    }
}
