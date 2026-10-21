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

//! Scaling of variables, residuals and linear systems, shared by the Newton
//! and homotopy solvers.
//!
//! Ports of the scaling in `nonlinearSolverHomotopy.c`: the variable scaling
//! `xScaling`, the residual scaling `resScaling` and `scaleMatrixRows`.

use faer::reborrow::{Reborrow, ReborrowMut};
use faer::{MatMut, MatRef};

/// Compute the variable scaling `factor[i] = max(nominal[i], |x[i]|)`.
///
/// Each unknown is scaled with its nominal value, or with its current
/// magnitude if that is larger. The solvers compute it at the start point and
/// use it to scale steps and Jacobian columns, so unknowns of very different
/// magnitudes are treated alike. From `xScaling` in `solveHomotopy`.
///
/// # Panics
///
/// Panics if `nominal`, `x` and `factor` don't have the same length.
pub fn variables(nominal: &[f64], x: &[f64], factor: &mut [f64]) {
    assert_eq!(factor.len(), nominal.len());
    assert_eq!(factor.len(), x.len());

    for ((out, x), nom) in factor.iter_mut().zip(x).zip(nominal) {
        *out = nom.max(x.abs());
    }
}

/// Compute the residual scaling `factor[i] = Σⱼ |J[i,j]|`, the L1 norm of each
/// row of the Jacobian.
///
/// `matrix` is the `n×(n+1)` matrix `[J | b]` of a linear system. Only the
/// Jacobian `J`, the first `n` columns, is used, the last column is ignored.
/// The solvers divide the residuals by these factors to measure the error of
/// residuals of very different magnitudes alike. From `resScaling`, computed
/// with `matVecMultAbsBB` and a vector of ones.
///
/// # Panics
///
/// Panics if `matrix` is not `n×(n+1)` or `factor` doesn't have `n` entries.
pub fn residuals(matrix: MatRef<'_, f64>, factor: &mut [f64]) {
    assert_eq!(factor.len(), matrix.nrows());
    assert_eq!(matrix.nrows() + 1, matrix.ncols(), "Matrix has wrong shape");

    let matrix = matrix.subcols(0, matrix.nrows());
    for (out, row) in factor.iter_mut().zip(matrix.row_iter()) {
        *out = row.norm_l1();
    }
}

/// Scale each row of the linear system `[J | b]` by the largest absolute value
/// of its row in `J`.
///
/// `matrix` is the `n×(n+1)` matrix `[J | b]`. The maximum is taken over the
/// first `n` columns only, but the whole row, including the right-hand side,
/// is divided by it. So the scaled system has the same solution, and its rows
/// have entries of comparable size, which improves the pivot search. A row
/// whose `J` part is zero is left unchanged. From `scaleMatrixRows`.
///
/// # Panics
///
/// Panics if `matrix` is not `n×(n+1)`.
pub fn rows(mut matrix: MatMut<'_, f64>) {
    const DELTA: f64 = 0.0;

    let n = matrix.nrows();
    assert_eq!(matrix.ncols(), n + 1, "Matrix has wrong shape");

    for mut row in matrix.rb_mut().row_iter_mut() {
        let maximum = row.rb().subcols(0, n).norm_max();
        let scale = if maximum <= DELTA { 1.0 } else { maximum };
        row /= scale;
    }
}
