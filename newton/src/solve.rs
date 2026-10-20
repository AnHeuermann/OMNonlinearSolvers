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

//! The Newton-Raphson algorithm, its workspace and its results.

// TODO: Remove once the port is complete
#![allow(dead_code)]

use crate::settings::Settings;
use faer::Mat;
use nls_model::{EvalError, ModelFunctions};

/// Buffers of the Newton-Raphson algorithm.
///
/// Allocated once and reused by every call of [`solve`], so the iterations
/// don't allocate. Sized for a system with `num_residuals` residuals and
/// `num_unknowns` unknowns.
pub struct NewtonWorkspace {
    /// Residual of the latest evaluated point, `F(x)` or `F(x1)`
    f: Box<[f64]>,
    /// Linear system `[J | f]` of the Newton step, `num_residuals` rows and
    /// `num_unknowns + 1` columns. The total pivot solver solves it as
    /// `J·dy + f = 0`, so the last column holds `f` itself.
    jacobian: Mat<f64>,
    /// Newton step, solution of `J(x)·dy = −f`
    dy: Box<[f64]>,
    /// Trial point `x + λ·dy`, candidate for next iteration
    x1: Box<[f64]>,

    /// Variable scaling, one factor per unknown, see
    /// [`nls_model::scaling::variables`]
    x_scaling: Box<[f64]>,
    /// Residual scaling, one factor per residual, see
    /// [`nls_model::scaling::residuals`]
    res_scaling: Box<[f64]>,
}

impl NewtonWorkspace {
    /// Allocates zero-initialized buffers for a system with `num_unknowns`
    /// unknowns and `num_residuals` residuals.
    pub fn new(num_unknowns: usize, num_residuals: usize) -> NewtonWorkspace {
        NewtonWorkspace {
            f: vec![0.0; num_residuals].into_boxed_slice(),
            jacobian: Mat::zeros(num_residuals, num_unknowns + 1),
            dy: vec![0.0; num_unknowns].into_boxed_slice(),
            x1: vec![0.0; num_unknowns].into_boxed_slice(),

            x_scaling: vec![0.0; num_unknowns].into_boxed_slice(),
            res_scaling: vec![0.0; num_residuals].into_boxed_slice(),
        }
    }
}

/// How a successful [`solve`] ended.
pub enum NewtonStatus {
    /// Converged: the residual and the last step are below their tolerances
    Solved,
    /// The steps became tiny while the residual is small, but not below its
    /// tolerance. The point is accepted as a solution with less accuracy.
    SolvedLessAccurate,
}

/// Why [`solve`] found no solution.
pub enum NewtonError {
    /// The maximum number of iterations of the [`Settings`] is reached
    MaxIterations,
    /// The steps became tiny while the residual is still large, for example
    /// in a local minimum of `‖F‖`
    Stuck,
    /// The iteration moves away from any solution: the error grew strongly,
    /// cycled or stagnated for many iterations
    Diverging,
    /// The Jacobian is singular and the linear system has no solution
    Singular,
    /// The model can't be evaluated: an assert or, for a casual tearing set,
    /// a violated constraint
    ModelError(EvalError),
}

/// Solve `F(x) = 0` for the `model` with a damped Newton-Raphson method.
///
/// `x` holds the start point on entry and the solution on success. `settings`
/// sets the tolerances and the maximum number of iterations, `work` provides
/// the buffers, sized for the model.
///
/// One iteration:
///
/// ```text
/// 1. f  := F(x)                      residual at the current point
/// 2. dy := solution of J(x)·dy = -f  the full Newton step
/// 3. x1 := x + λ·dy                  trial point, λ ∈ (0, 1] from the line search
/// 4. f  := F(x1)                     residual at the trial point
/// 5. accepted?  x := x1              the trial point becomes the new iterate
/// ```
///
/// From `newtonAlgorithm` of the C runtime.
///
/// # Errors
///
/// Returns a [`NewtonError`] if the iteration doesn't converge, the Jacobian
/// is singular or the model can't be evaluated.
pub fn solve<M: ModelFunctions>(
    model: &mut M,
    x: &mut [f64],
    settings: &Settings,
    work: &mut NewtonWorkspace,
) -> Result<NewtonStatus, NewtonError> {
    #![allow(unused_variables)]
    todo!()
}

#[cfg(test)]
mod newton_workspace_tests {
    use super::*;

    #[test]
    /// `NewtonWorkspace::new` allocates zero-initialized buffers of the
    /// documented sizes. Different numbers of unknowns and residuals catch
    /// swapped dimensions.
    fn newton_workspace_new_test() {
        let (num_unknowns, num_residuals) = (3, 2);

        let work = NewtonWorkspace::new(num_unknowns, num_residuals);

        assert_eq!(*work.f, [0.0; 2], "f: one entry per residual");
        assert_eq!(*work.dy, [0.0; 3], "dy: one entry per unknown");
        assert_eq!(*work.x1, [0.0; 3], "x1: one entry per unknown");
        assert_eq!(
            *work.x_scaling, [0.0; 3],
            "x_scaling: one entry per unknown"
        );
        assert_eq!(
            *work.res_scaling, [0.0; 2],
            "res_scaling: one entry per residual"
        );

        // [J | f]: one row per residual, one column per unknown and one for f
        assert_eq!(
            (work.jacobian.nrows(), work.jacobian.ncols()),
            (num_residuals, num_unknowns + 1),
            "jacobian: shape of [J | f]"
        );
        assert_eq!(work.jacobian, Mat::<f64>::zeros(2, 4));
    }
}
