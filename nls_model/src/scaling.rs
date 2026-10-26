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
//! `xScaling`, the residual scaling `resScaling`, `scaleMatrixRows` and the
//! column scaling of the Jacobian.

use faer::reborrow::{Reborrow, ReborrowMut};
use faer::{MatMut, MatRef};
use log::trace;

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
    assert_eq!(
        factor.len(),
        matrix.nrows(),
        "Scaling vector has wrong length"
    );
    assert_eq!(matrix.nrows() + 1, matrix.ncols(), "Matrix has wrong shape");

    let jacobian = matrix.subcols(0, matrix.nrows());
    for (out, row) in factor.iter_mut().zip(jacobian.row_iter()) {
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

/// Scale each column `j` of the Jacobian in the linear system `[J | b]` by
/// `factor[j]`.
///
/// `matrix` is the `n×(n+1)` matrix `[J | b]` and `factor` the [`variables`]
/// scaling of the `n` unknowns. Scaling the unknowns as `x = factor·x̃` turns
/// `J·x = −b` into `(J·diag(factor))·x̃ = −b`, so only the columns of `J`
/// change and the right-hand side `b` stays unchanged. A solution `x̃` of the
/// scaled system must be multiplied by `factor` to get `x`. From the column
/// scaling in `getAnalyticalJacobianHomotopy` and
/// `getNumericalJacobianHomotopy`.
///
/// # Panics
///
/// Panics if `matrix` is not `n×(n+1)` or `factor` doesn't have `n` entries.
pub fn columns(matrix: MatMut<'_, f64>, factor: &[f64]) {
    let n = matrix.nrows();
    assert_eq!(factor.len(), n, "Scaling vector has wrong length");
    assert_eq!(matrix.ncols(), n + 1, "Matrix has wrong shape");

    let jacobian = matrix.subcols_mut(0, n);
    for (mut col, factor) in jacobian.col_iter_mut().zip(factor) {
        col *= factor;
    }
}

/// Replace infinite or NaN scaling factors with `1.0`, and log one warning if
/// any were replaced.
///
/// A non-finite residual scaling would make the scaled error `0` or NaN, for
/// example after a Jacobian entry overflowed. From `vecMakeFinite`.
pub fn make_finite(x: &mut [f64]) {
    let mut changed = false;
    for elem in x.iter_mut() {
        if !elem.is_finite() {
            changed = true;
            *elem = 1.0;
        }
    }
    if changed {
        trace!(target: "NLS", "Entry of the residual scaling is inf or NaN, set to 1.0");
    }
}

/// Scale a solution `x̃` of a system with [`columns`] scaling back to the
/// unknowns, `x = |factor|·x̃`. Entries with a factor of zero stay unscaled,
/// like `vecMultScaling` in the C runtime.
///
/// # Panics
///
/// Panics if `scaled`, `factor` and `x` don't have the same length.
pub fn unscale(scaled: &[f64], factor: &[f64], x: &mut [f64]) {
    assert_eq!(scaled.len(), factor.len());
    assert_eq!(scaled.len(), x.len());

    for ((x, scaled), factor) in x.iter_mut().zip(scaled).zip(factor) {
        *x = if factor.abs() > 0.0 {
            scaled * factor.abs()
        } else {
            *scaled
        };
    }
}

/// Squared error `‖f‖²` of a residual
pub fn squared_error(f: &[f64]) -> f64 {
    f.iter().map(|f| f * f).sum()
}

/// Squared scaled error `‖f / |res_scaling|‖²` of a residual. Entries with a
/// scaling of zero stay unscaled, like `vecDivScaling` in the C runtime.
pub fn scaled_squared_error(f: &[f64], res_scaling: &[f64]) -> f64 {
    debug_assert_eq!(f.len(), res_scaling.len());

    f.iter()
        .zip(res_scaling)
        .map(|(f, scaling)| {
            let scaled = if scaling.abs() > 0.0 {
                f / scaling.abs()
            } else {
                *f
            };
            scaled * scaled
        })
        .sum()
}
