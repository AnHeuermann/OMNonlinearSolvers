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

struct Position {
    row: usize,
    column: usize,
}

/// Calculate pivot element on permuted matrix.
///
/// At elimination step `step` find the entry with the larges absolute value in the
/// not yet eliminated part.
/// Return `position: Position` of pivot element and its absolute value
fn pivot_search(a: &Mat<f64>, step: usize) -> (Position, f64) {
    return (Position { row: 0, column: 0 }, 0.0);
}
