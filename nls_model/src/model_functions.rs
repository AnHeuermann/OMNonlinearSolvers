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

//! The interface of a nonlinear system and the error of a failed evaluation.

use faer::MatMut;

use crate::SystemInfo;

/// Why evaluating a model failed.
#[derive(Debug)]
pub enum EvalError {
    /// The model raised an assert, e.g. a function outside its domain
    Assert,
    /// A constraint of a casual tearing set is violated
    ConstraintViolated,
}

/// A nonlinear system `F(x) = 0`.
///
/// Implemented by the type of each concrete system, its fields hold the data
/// needed to evaluate it. The facts about the system, its sizes, nominal
/// values, bounds and tearing kind, are a [`SystemInfo`] returned by
/// [`info`](Self::info). A system with an additional unknown, like a homotopy
/// `H(x, λ)` with `n` residuals and `n+1` unknowns, is non-square.
///
/// [`info`](Self::info) and [`eval`](Self::eval) are required,
/// [`derivative`](Self::derivative) defaults to a numerical Jacobian.
pub trait ModelFunctions {
    /// Facts about the system: sizes, nominal values, bounds and tearing kind.
    fn info(&self) -> &SystemInfo;

    /// Evaluate the residual `res := F(x)`.
    ///
    /// `x` has `info().num_unknowns` entries, `res` has `info().num_residuals`
    /// entries.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if `F` can't be evaluated at `x`.
    fn eval(&mut self, x: &[f64], res: &mut [f64]) -> Result<(), EvalError>;

    /// Evaluate the Jacobian `jacobian := J(x)` with `J = ∂F/∂x`.
    ///
    /// `fx` is the residual `F(x)` at the same `x`, the caller has computed it
    /// already. `jacobian` is a `num_residuals × num_unknowns` view of a
    /// column-major matrix.
    ///
    /// The default computes forward differences, one evaluation of `F` per
    /// column, with the step `h = sqrt(20·ε)·(1 + |xⱼ|)` of the C runtime.
    /// Unlike the C runtime it doesn't respect the bounds of `x` yet. Systems
    /// with an analytical Jacobian override it.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if `F` can't be evaluated at a perturbed `x`.
    ///
    /// # Panics
    ///
    /// Panics if the shape of `jacobian` doesn't match `x` and `fx`.
    fn derivative(
        &mut self,
        x: &[f64],
        fx: &[f64],
        jacobian: MatMut<'_, f64>,
    ) -> Result<(), EvalError> {
        assert_eq!(x.len(), jacobian.ncols());
        assert_eq!(fx.len(), jacobian.nrows());

        let delta = (f64::EPSILON * 20.0).sqrt();

        let mut x_tilde = x.to_vec();

        for (j, col) in jacobian.col_iter_mut().enumerate() {
            let h = delta * (1.0 + x[j].abs());
            x_tilde[j] = x[j] + h;

            // col := F(x + h·eⱼ)
            let column: &mut [f64] = col.try_as_col_major_mut().unwrap().as_slice_mut();
            self.eval(&x_tilde, column)?;

            // col := (F(x + h·eⱼ) − F(x)) / h
            for i in 0..column.len() {
                column[i] = (column[i] - fx[i]) / h;
            }

            x_tilde[j] = x[j]; // undo the perturbation
        }

        Ok(())
    }
}
