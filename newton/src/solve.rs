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

use crate::line_search::{LineSearchError, Start, line_search};
use crate::settings::Settings;
use crate::tracker::{Decision, Step, Tracker};
use faer::{ColRef, Mat};
use log::{debug, trace};
use nls_model::scaling::{self, scaled_squared_error, squared_error};
use nls_model::{EvalError, ModelFunctions};
use total_pivot_solver::{SolveError, SolveStatus};

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

    /// The residual scaling of the latest Jacobian of [`solve`].
    pub fn res_scaling(&self) -> &[f64] {
        &self.res_scaling
    }
}

/// How a successful [`solve`] ended.
#[derive(Debug, PartialEq)]
pub enum NewtonStatus {
    /// Converged: the residual and the last step are below their tolerances
    Solved,
    /// The steps became tiny while the residual is small, but not below its
    /// tolerance. The point is accepted as a solution with less accuracy.
    SolvedLessAccurate,
}

/// Why [`solve`] found no solution.
#[derive(Debug)]
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

/// A model that fails to evaluate is a [`NewtonError::ModelError`]
impl From<EvalError> for NewtonError {
    fn from(error: EvalError) -> Self {
        NewtonError::ModelError(error)
    }
}

/// A linear system the total pivot solver can't solve is singular
impl From<SolveError> for NewtonError {
    fn from(_: SolveError) -> Self {
        NewtonError::Singular
    }
}

/// A model that asserted at every step factor of the line search fails like a
/// single assert
impl From<LineSearchError> for NewtonError {
    fn from(error: LineSearchError) -> Self {
        match error {
            LineSearchError::BackedOffBelowLambdaMin => NewtonError::ModelError(EvalError::Assert),
            LineSearchError::ModelError(error) => NewtonError::ModelError(error),
        }
    }
}

/// Solve `F(x) = 0` for the `model` with a damped Newton-Raphson method.
///
/// `x` holds the start point on entry. On success it holds the solution, on
/// failure the last accepted iterate. `settings` sets the tolerances and the
/// maximum number of iterations, `work` provides the buffers, sized for the
/// model.
///
/// One iteration:
///
/// ```text
/// 1. dy := solution of J(x)·dy = −F(x)   the full Newton step
/// 2. x1 := x + λ·dy                      trial point, λ ∈ (0, 1] from the line search
/// 3. stop?                               see the tracker's decision
/// 4. x  := x1, refresh J(x)              the trial point becomes the new iterate
/// ```
///
/// The unknowns are scaled with their nominal values or magnitudes at the
/// start point, the residuals with the row sums of the Jacobian. From
/// `newtonAlgorithm` of the C runtime, together with the setup in
/// `solveHomotopy` before it.
///
/// # Errors
///
/// Returns a [`NewtonError`] if the iteration doesn't converge, the Jacobian
/// is singular or the model can't be evaluated.
///
/// # Panics
///
/// Panics if `x` or `work` don't match the size of the model.
pub fn solve<M: ModelFunctions>(
    model: &mut M,
    x: &mut [f64],
    settings: &Settings,
    work: &mut NewtonWorkspace,
) -> Result<NewtonStatus, NewtonError> {
    let n = model.info().num_unknowns;
    let casual_tearing_set = model.info().casual_tearing_set;
    assert_eq!(x.len(), n, "Start point has wrong length");
    assert_eq!(work.x1.len(), n, "Workspace has wrong size");

    scaling::variables(&model.info().nominal, x, &mut work.x_scaling);
    evaluate_start(model, x, work)?;

    let mut err = squared_error(&work.f);
    let mut err_scaled = scaled_squared_error(&work.f, &work.res_scaling);
    let mut tracker = Tracker::new(settings, err, err_scaled);
    debug!(target: "NLS_V", "Newton iteration started, error {err:e}, scaled {err_scaled:e}");

    for iteration in 1.. {
        let solution = total_pivot_solver::solve(&mut work.jacobian, Some(n), casual_tearing_set)?;
        // The C runtime ignores a vanishing determinant of the first linear
        // system, which `solveHomotopy` solves before the iteration starts
        let determinant_zero =
            iteration > 1 && matches!(solution.status, SolveStatus::SolvedDeterminantZero);
        scaling::unscale(&solution.x[..n], &work.x_scaling, &mut work.dy);

        let start = Start {
            err,
            err_scaled,
            g: -2.0 * err,
            g_scaled: -2.0 * err_scaled,
        };
        let mut eval = |x: &[f64], f: &mut [f64]| model.eval(x, f);
        let step_factor = line_search(
            &mut eval,
            x,
            &work.dy,
            &work.res_scaling,
            &start,
            &mut work.x1,
            &mut work.f,
        )?;

        err = squared_error(&work.f);
        err_scaled = scaled_squared_error(&work.f, &work.res_scaling);
        let step = Step {
            err,
            err_scaled,
            dy: ColRef::from_slice(&work.dy).squared_norm_l2(),
            dy_scaled: scaled_squared_error(&work.dy, &work.x_scaling),
            step_factor,
        };
        tracker.update(step);
        trace!(
            target: "NLS_V",
            "Newton iteration {iteration}: error {err:e}, scaled {err_scaled:e}, ‖dy‖² {:e}, λ {step_factor}",
            step.dy
        );

        match tracker.decide(determinant_zero) {
            Decision::Continue => {}
            Decision::Converged { status, take_step } => {
                if take_step {
                    x.copy_from_slice(&work.x1);
                }
                debug!(target: "NLS_V", "Newton iteration {iteration}: {status:?}");
                return Ok(status);
            }
            Decision::Failed(error) => {
                debug!(target: "NLS_V", "Newton iteration {iteration}: {error:?}");
                return Err(error);
            }
        }

        x.copy_from_slice(&work.x1);
        refresh_jacobian(model, x, work)?;
    }
    unreachable!("the tracker stops after the maximum number of iterations")
}

/// Evaluate the start point `x` and set up the first linear system.
///
/// Writes `F(x)` into `work.f`, and the scaled `[J | f]` and the residual
/// scaling as [`refresh_jacobian`] does. Non-finite factors of the residual
/// scaling are set to 1, like `solveHomotopy` does before the C Newton
/// iteration starts.
///
/// # Errors
///
/// Returns the [`EvalError`] of the model if `F` or `J` can't be evaluated at
/// `x`.
fn evaluate_start<M: ModelFunctions>(
    model: &mut M,
    x: &[f64],
    work: &mut NewtonWorkspace,
) -> Result<(), EvalError> {
    model.eval(x, &mut work.f)?;
    refresh_jacobian(model, x, work)?;
    scaling::make_finite(&mut work.res_scaling);
    Ok(())
}

/// Set up the linear system `[J | f]` of the next Newton step at `x`.
///
/// Expects `F(x)` in `work.f`, as the line search leaves it. Evaluates the
/// Jacobian `J(x)` into `work.jacobian`, with `f` as its last column, scales
/// the columns with `work.x_scaling`, computes `work.res_scaling` from the
/// column-scaled `J` and finally scales the rows. From the end of the loop in
/// `newtonAlgorithm`.
///
/// # Errors
///
/// Returns the [`EvalError`] of the model if `J` can't be evaluated at `x`.
/// `work.jacobian` is unspecified then.
fn refresh_jacobian<M: ModelFunctions>(
    model: &mut M,
    x: &[f64],
    work: &mut NewtonWorkspace,
) -> Result<(), EvalError> {
    let jacobian = work.jacobian.subcols_mut(0, work.jacobian.nrows());

    model.derivative(x, &work.f, jacobian)?;
    work.jacobian
        .col_mut(work.jacobian.nrows())
        .copy_from(ColRef::from_slice(&work.f));

    scaling::columns(work.jacobian.as_mut(), &work.x_scaling);
    scaling::residuals(work.jacobian.as_ref(), &mut work.res_scaling);
    scaling::rows(work.jacobian.as_mut());
    Ok(())
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
