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

//! The strategy of `solveHomotopy` in the C runtime: which solver runs when.

use crate::path_follower::{Direction, follow_path};
use crate::settings::Settings;
use crate::{FixpointHomotopy, NewtonHomotopy};
use faer::{ColRef, Mat};
use log::debug;
use newton::NewtonWorkspace;
use newton::settings::Settings as NewtonSettings;
use nls_model::scaling::{self, scaled_squared_error, squared_error};
use nls_model::{EvalError, ModelFunctions};
use total_pivot_solver::SolveError;

/// A start point is already a solution if its squared error is below this
/// factor times `f_tol²`
const SOLVED_AT_START_FACTOR: f64 = 1e-4;
/// The start point moves by these factors of the variable scaling on the
/// first and second retry
const RETRY_SHIFTS: [f64; 2] = [0.01, 0.1];

/// A homotopy the strategy follows after Newton's method failed, in the order
/// of the runs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HomotopyRun {
    /// The Newton homotopy, starting towards increasing `λ`
    NewtonPositive,
    /// The Newton homotopy, starting towards decreasing `λ`
    NewtonNegative,
    /// The fixpoint homotopy, starting towards increasing `λ`
    Fixpoint,
}

/// The runs in the order of the C runtime
const RUNS: [HomotopyRun; 3] = [
    HomotopyRun::NewtonPositive,
    HomotopyRun::NewtonNegative,
    HomotopyRun::Fixpoint,
];

/// How [`HomotopySolver::solve`] found the solution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Solved {
    /// The start point already solves the system
    StartPoint,
    /// Newton's method from the start point
    Newton,
    /// Newton's method from the end of the path of a homotopy
    Homotopy(HomotopyRun),
}

/// Why [`HomotopySolver`] found no solution.
#[derive(Debug)]
pub enum HomotopyError {
    /// The model can't be evaluated at the start point, also not at the
    /// moved start points, or it violates a constraint there
    ModelError(EvalError),
    /// Newton's method and all homotopy runs failed
    NotSolved,
}

/// Why the linear system at a point failed.
enum PointError {
    ModelError(EvalError),
    Singular,
}

impl From<EvalError> for PointError {
    fn from(error: EvalError) -> Self {
        PointError::ModelError(error)
    }
}

impl From<SolveError> for PointError {
    fn from(_: SolveError) -> Self {
        PointError::Singular
    }
}

/// The start point after its checks.
enum Start {
    /// The start point solves the system
    Solved,
    /// Newton's method and the homotopies can start from it
    Ready,
}

/// The damped Newton solver combined with homotopies, `solveHomotopy` of the
/// C runtime.
///
/// It keeps the buffers and the residual scaling of the latest successful
/// solve between calls, like the C runtime keeps them in its solver data.
pub struct HomotopySolver {
    newton_settings: NewtonSettings,
    settings: Settings,
    newton: NewtonWorkspace,
    /// The start point given by the caller
    x_start: Box<[f64]>,
    /// The start point of Newton's method and the homotopies, maybe moved
    x0: Box<[f64]>,
    /// `F(x₀)`
    fx0: Box<[f64]>,
    /// Variable scaling of the `n` unknowns and `λ`
    x_scaling: Box<[f64]>,
    /// The linear system `[J | f]` the start point is checked with
    jacobian: Mat<f64>,
    /// Residual scaling of the latest successful solve, zero otherwise. Only
    /// the test of an already solved start point uses it.
    res_scaling: Box<[f64]>,
}

impl HomotopySolver {
    /// Creates a solver for systems with `size` residuals and unknowns.
    ///
    /// # Panics
    ///
    /// Panics if `size` is zero.
    pub fn new(size: usize, newton_settings: NewtonSettings, settings: Settings) -> Self {
        assert!(size > 0, "Invalid size {size}");
        HomotopySolver {
            newton_settings,
            settings,
            newton: NewtonWorkspace::new(size, size),
            x_start: vec![0.0; size].into(),
            x0: vec![0.0; size].into(),
            fx0: vec![0.0; size].into(),
            x_scaling: vec![1.0; size + 1].into(),
            jacobian: Mat::zeros(size, size + 1),
            res_scaling: vec![0.0; size].into(),
        }
    }

    /// Solve `F(x) = 0` for `model` from the start point `x`.
    ///
    /// 1. A start point that already solves the system is taken. Otherwise the
    ///    model must be evaluated there and its first linear system must be
    ///    solvable, else the start point is moved by `1%` and then `10%` of
    ///    the variable scaling.
    /// 2. Newton's method from the start point.
    /// 3. If that fails, the paths of the Newton homotopy in both directions
    ///    and of the fixpoint homotopy, each from the start point, until a
    ///    path ends at a point from which Newton's method converges.
    ///
    /// A casual tearing set gives up after the first failure. `x` holds the
    /// solution on success and is unchanged otherwise. From `solveHomotopy`,
    /// without the hybrd fallback and without the handling of mixed systems,
    /// which belongs to the C interface.
    ///
    /// # Errors
    ///
    /// Returns a [`HomotopyError`] if the model can't be evaluated at the start
    /// point or no solution is found.
    ///
    /// # Panics
    ///
    /// Panics if `model` or `x` don't have the size of the solver.
    pub fn solve<M: ModelFunctions>(
        &mut self,
        model: &mut M,
        x: &mut [f64],
    ) -> Result<Solved, HomotopyError> {
        let n = self.x_start.len();
        assert_eq!(model.info().num_unknowns, n, "Model has wrong size");
        assert_eq!(model.info().num_residuals, n, "Model has wrong size");
        assert_eq!(x.len(), n, "Start point has wrong length");

        self.x_start.copy_from_slice(x);
        scaling::variables(&model.info().nominal, x, &mut self.x_scaling[..n]);
        self.x_scaling[n] = 1.0;

        let result = self.solve_from_start(model, x);
        match &result {
            Ok(solved) => debug!(target: "NLS", "Solved with {solved:?}"),
            Err(error) => {
                debug!(target: "NLS", "Not solved: {error:?}");
                x.copy_from_slice(&self.x_start);
                self.res_scaling.fill(0.0);
            }
        }
        result
    }

    /// Solve with a model that has its own homotopy, the `homotopy()`
    /// operator of Modelica, from the start point `x`.
    ///
    /// `model` has `λ` as last unknown. Its path is followed from `(x, 0)` in
    /// the `first` direction, and if that fails in the opposite one. The end
    /// of the path is the solution, there is no Newton's method before or
    /// after it. Returns the direction of the successful path. From
    /// `solveHomotopy` with `initHomotopy`.
    ///
    /// # Errors
    ///
    /// Returns [`HomotopyError::NotSolved`] if both paths fail. `x` is
    /// unchanged then.
    ///
    /// # Panics
    ///
    /// Panics if `model` doesn't have the size of the solver with `λ` as
    /// additional unknown, or `x` doesn't have the size of the solver.
    pub fn solve_own_homotopy<M: ModelFunctions>(
        &mut self,
        model: &mut M,
        x: &mut [f64],
        first: Direction,
    ) -> Result<Direction, HomotopyError> {
        let n = self.x_start.len();
        assert_eq!(model.info().num_residuals, n, "Model has wrong size");
        assert_eq!(
            model.info().num_unknowns,
            n + 1,
            "Model needs λ as last unknown"
        );
        assert_eq!(x.len(), n, "Start point has wrong length");

        self.x_start.copy_from_slice(x);
        scaling::variables(&model.info().nominal[..n], x, &mut self.x_scaling[..n]);
        self.x_scaling[n] = 1.0;

        for direction in [first, first.opposite()] {
            x.copy_from_slice(&self.x_start);
            let path = follow_path(
                model,
                &self.x_scaling,
                direction,
                self.newton_settings.f_tol(),
                &self.settings,
                x,
            );
            match path {
                Ok(_) => return Ok(direction),
                Err(error) => {
                    debug!(target: "NLS", "Homotopy in {direction:?} direction: {error:?}")
                }
            }
        }
        x.copy_from_slice(&self.x_start);
        self.res_scaling.fill(0.0);
        Err(HomotopyError::NotSolved)
    }

    /// The steps of [`HomotopySolver::solve`] after the setup.
    fn solve_from_start<M: ModelFunctions>(
        &mut self,
        model: &mut M,
        x: &mut [f64],
    ) -> Result<Solved, HomotopyError> {
        let casual_tearing_set = model.info().casual_tearing_set;
        let start = self.check_start(model)?;
        x.copy_from_slice(&self.x0);
        if let Start::Solved = start {
            return Ok(Solved::StartPoint);
        }

        if self.newton(model, x) {
            return Ok(Solved::Newton);
        }
        if casual_tearing_set {
            return Err(HomotopyError::NotSolved);
        }

        for run in RUNS {
            x.copy_from_slice(&self.x0);
            if let Err(error) = self.follow(model, run, x) {
                debug!(target: "NLS", "Homotopy {run:?}: {error:?}");
                continue;
            }
            // The C runtime gives up if the end of a path can't start Newton's
            // method
            model
                .eval(x, &mut self.fx0)
                .map_err(|_| HomotopyError::NotSolved)?;
            check_linear_system(model, x, &self.fx0, &self.x_scaling, &mut self.jacobian)
                .map_err(|_| HomotopyError::NotSolved)?;
            if self.newton(model, x) {
                return Ok(Solved::Homotopy(run));
            }
        }
        Err(HomotopyError::NotSolved)
    }

    /// Check the start point and move it while the model can't be evaluated
    /// there or its linear system is singular. Leaves the start point to use
    /// in `x0` and its residual in `fx0`.
    fn check_start<M: ModelFunctions>(&mut self, model: &mut M) -> Result<Start, HomotopyError> {
        let n = self.x_start.len();
        let casual_tearing_set = model.info().casual_tearing_set;
        let f_tol_squared = self.newton_settings.f_tol().powi(2);

        for attempt in 0..=RETRY_SHIFTS.len() {
            let last_attempt = attempt == RETRY_SHIFTS.len();
            self.x0.copy_from_slice(&self.x_start);
            if attempt > 0 {
                let shift = RETRY_SHIFTS[attempt - 1];
                for (i, (x0, scaling)) in self.x0.iter_mut().zip(&self.x_scaling).enumerate() {
                    *x0 += scaling * i as f64 / n as f64 * shift;
                }
            }

            let checked = model
                .eval(&self.x0, &mut self.fx0)
                .map_err(PointError::from);
            let checked = checked.and_then(|()| {
                let solved = squared_error(&self.fx0) < f_tol_squared * SOLVED_AT_START_FACTOR
                    || scaled_squared_error(&self.fx0, &self.res_scaling)
                        < f_tol_squared * SOLVED_AT_START_FACTOR;
                if solved {
                    return Ok(Start::Solved);
                }
                check_linear_system(
                    model,
                    &self.x0,
                    &self.fx0,
                    &self.x_scaling,
                    &mut self.jacobian,
                )?;
                Ok(Start::Ready)
            });

            match checked {
                Ok(start) => return Ok(start),
                Err(PointError::ModelError(error @ EvalError::ConstraintViolated)) => {
                    return Err(HomotopyError::ModelError(error));
                }
                Err(PointError::ModelError(error)) if casual_tearing_set || last_attempt => {
                    return Err(HomotopyError::ModelError(error));
                }
                Err(PointError::Singular) if casual_tearing_set => {
                    return Err(HomotopyError::NotSolved);
                }
                // Like the C runtime, a singular system at the last moved
                // start point is left to Newton's method and the homotopies
                Err(PointError::Singular) if last_attempt => return Ok(Start::Ready),
                Err(_) => {}
            }
        }
        unreachable!("the last attempt returns")
    }

    /// Newton's method from `x`. Keeps its residual scaling on success.
    fn newton<M: ModelFunctions>(&mut self, model: &mut M, x: &mut [f64]) -> bool {
        match newton::solve(model, x, &self.newton_settings, &mut self.newton) {
            Ok(_) => {
                self.res_scaling.copy_from_slice(self.newton.res_scaling());
                true
            }
            Err(error) => {
                debug!(target: "NLS", "Newton's method: {error:?}");
                false
            }
        }
    }

    /// Follow the path of the homotopy of `run` from `x₀`, writing its end
    /// into `x`.
    fn follow<M: ModelFunctions>(
        &mut self,
        model: &mut M,
        run: HomotopyRun,
        x: &mut [f64],
    ) -> Result<usize, crate::PathError> {
        let n = self.x0.len();
        let f_tol = self.newton_settings.f_tol();
        match run {
            HomotopyRun::NewtonPositive | HomotopyRun::NewtonNegative => {
                let direction = if run == HomotopyRun::NewtonPositive {
                    Direction::Positive
                } else {
                    Direction::Negative
                };
                let mut homotopy = NewtonHomotopy::new(model, &self.fx0);
                follow_path(
                    &mut homotopy,
                    &self.x_scaling,
                    direction,
                    f_tol,
                    &self.settings,
                    x,
                )
            }
            HomotopyRun::Fixpoint => {
                let mut homotopy = FixpointHomotopy::new(model, &self.x0, &self.x_scaling[..n]);
                follow_path(
                    &mut homotopy,
                    &self.x_scaling,
                    Direction::Positive,
                    f_tol,
                    &self.settings,
                    x,
                )
            }
        }
    }
}

/// Evaluate the Jacobian at `x` with the residual `f = F(x)` and solve its
/// scaled linear system `[J | f]`, the first Newton step from `x`.
///
/// The C runtime solves it as the first step of its Newton iteration. Here
/// [`newton::solve`] sets up its own, so it only tells whether Newton's method
/// can start from `x`.
fn check_linear_system<M: ModelFunctions>(
    model: &mut M,
    x: &[f64],
    f: &[f64],
    x_scaling: &[f64],
    jacobian: &mut Mat<f64>,
) -> Result<(), PointError> {
    let n = x.len();
    model.derivative(x, f, jacobian.as_mut().subcols_mut(0, n))?;
    jacobian.col_mut(n).copy_from(ColRef::from_slice(f));
    scaling::columns(jacobian.as_mut(), &x_scaling[..n]);
    scaling::rows(jacobian.as_mut());
    total_pivot_solver::solve(jacobian, Some(n), model.info().casual_tearing_set)?;
    Ok(())
}
