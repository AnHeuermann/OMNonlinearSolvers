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

use crate::solve::SolveError;
use permuted_matrix::PermutedMatrix;

/// Back substitution on a matrix after [`forward_elimination`].
///
/// Solves the upper triangular system in logical order from the last row up
/// and returns `x` with `n+1` entries in the original (physical) column order.
/// The entry of the fixed column (logical column `n`) is set to 1.
///
/// Rows from `rank` on have no pivot. If their right-hand side is (almost)
/// zero, `|rhs| <= 1e-6`, the system is consistent and the unknown of that row
/// is set to 0, otherwise the system has no solution.
///
/// # Errors
///
/// Returns [`SolveError::UnderDeterminde`] if a row without pivot has a
/// right-hand side larger than `1e-6`.
///
/// [`forward_elimination`]: crate::forward_elimination::forward_elimination
pub(crate) fn backward_elimination(
    matrix: &PermutedMatrix<f64>,
    rank: usize,
) -> Result<Box<[f64]>, SolveError> {
    let n_rows = matrix.nrows(); // n
    let n_cols = matrix.ncols(); // n+1

    let mut x: Box<[f64]> = vec![0.0; n_cols].into_boxed_slice();
    for i in (0..n_rows).rev() {
        if i >= rank {
            if matrix[(i, n_rows)].abs() > 1e-6 {
                return Err(SolveError::UnderDeterminde);
            } else {
                let col_index_i = matrix.column_permutation()[i];
                x[col_index_i] = 0.0;
            }
        } else {
            let col_index_i = matrix.column_permutation()[i];
            x[col_index_i] = -matrix[(i, n_rows)];

            for j in (i + 1..n_rows).rev() {
                let col_index_j = matrix.column_permutation()[j];
                x[col_index_i] -= matrix[(i, j)] * x[col_index_j];
            }
            x[col_index_i] /= matrix[(i, i)];
        }
    }
    let col_index_n = matrix.column_permutation()[n_rows];
    x[col_index_n] = 1.0;

    Ok(x)
}

#[cfg(test)]
mod backward_elimination_tests {
    use super::*;
    use equator::assert;
    use faer::mat;
    use faer::utils::approx::{ApproxEq, CwiseMat};
    use faer::{Col, ColRef};

    use permuted_matrix::PermutedMatrix;

    #[test]
    fn backward_elimination_test() {
        let mut m = mat![
            [1.0, 3.0, 1.0, -9.0],
            [0.0, -2.0, -2.0, 8.0],
            [0.0, 0.0, 1.0, 2.0]
        ];
        let p = PermutedMatrix::new(&mut m);

        let x = backward_elimination(&p, 3).unwrap();

        let residual = &m * ColRef::from_slice(&x);
        let approx_eq = CwiseMat(ApproxEq::eps());
        assert!(residual ~ Col::<f64>::zeros(3));
    }

    #[test]
    /// Rank-deficient but consistent: the Wikipedia example
    /// <https://en.wikipedia.org/wiki/Gaussian_elimination> after the forward
    /// elimination, with rank 2. The last row is zero, including the right-hand
    /// side.
    ///
    /// The unknown without a pivot, `x₃`, is set to 0. That gives
    /// `x = [-3, 4, 0]`, the solution `[-3 + 2·x₃, 4 - x₃, x₃]` of the
    /// original system for `x₃ = 0`.
    fn backward_elimination_rank_deficient_test() {
        let mut m = mat![
            [1.0, 3.0, 1.0, -9.0],
            [0.0, -2.0, -2.0, 8.0],
            [0.0, 0.0, 0.0, 0.0]
        ];
        let p = PermutedMatrix::new(&mut m);

        let x = backward_elimination(&p, 2).unwrap();

        assert!(*x == [-3.0, 4.0, 0.0, 1.0]);

        // x also solves the original system before the elimination
        let original = mat![
            [1.0, 3.0, 1.0, -9.0],
            [1.0, 1.0, -1.0, -1.0],
            [3.0, 11.0, 5.0, -35.0]
        ];
        let residual = &original * ColRef::from_slice(&x);
        let approx_eq = CwiseMat(ApproxEq::eps());
        assert!(residual ~ Col::<f64>::zeros(3));
    }

    #[test]
    /// Rank-deficient and inconsistent: the same matrix, but the zero row has
    /// a right-hand side of 5. Its equation reads `0 = 5`, so there is no
    /// solution.
    fn backward_elimination_under_determined_test() {
        let mut m = mat![
            [1.0, 3.0, 1.0, -9.0],
            [0.0, -2.0, -2.0, 8.0],
            [0.0, 0.0, 0.0, 5.0]
        ];
        let p = PermutedMatrix::new(&mut m);

        let result = backward_elimination(&p, 2);

        assert!(matches!(result, Err(SolveError::UnderDeterminde)));
    }
}
