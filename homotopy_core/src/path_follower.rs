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

//! Following the path of a homotopy from `λ = 0` to `λ = 1` with predictor
//! and corrector steps, `homotopyAlgorithm` of the C runtime.

use crate::corrector::{Correction, CorrectorError, CorrectorWorkspace, corrector};
use crate::predictor::{PredictorError, limit_to_end, predictor};
use crate::settings::Settings;
use crate::step_control::{Control, Outcome, StepFailure, bend, control};
use crate::tangent::{TangentError, orient, tangent};
use faer::Mat;
use log::trace;
use nls_model::{EvalError, ModelFunctions};

/// The direction the path starts in from `λ = 0`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Direction {
    /// Towards increasing `λ`
    Positive,
    /// Towards decreasing `λ`
    Negative,
}

impl Direction {
    /// `+1` or `−1`, the `startDirection` of the C runtime
    fn sign(self) -> f64 {
        match self {
            Direction::Positive => 1.0,
            Direction::Negative => -1.0,
        }
    }
}

/// Why the path following didn't reach `λ = 1`.
#[derive(Debug)]
pub enum PathError {
    /// The homotopy or its Jacobian can't be evaluated
    ModelError(EvalError),
    /// The linear system of the tangent has no solution
    Singular,
    /// The homotopy asserts at every predictor point down to `tau_min`
    PredictorBackedOff,
    /// The corrector neither reached the path nor moved
    CorrectorStuck,
    /// A step was tried `max_tries` times
    TooManyTries,
    /// The path went below `λ = −1`
    LambdaBelowMinusOne,
    /// The path didn't reach `λ = 1` within `max_lambda_steps` steps
    TooManySteps,
}

impl From<EvalError> for PathError {
    fn from(error: EvalError) -> Self {
        PathError::ModelError(error)
    }
}

impl From<TangentError> for PathError {
    fn from(error: TangentError) -> Self {
        match error {
            TangentError::ModelError(error) => PathError::ModelError(error),
            TangentError::Singular => PathError::Singular,
        }
    }
}

impl From<PredictorError> for PathError {
    fn from(error: PredictorError) -> Self {
        match error {
            PredictorError::BackedOffBelowTauMin => PathError::PredictorBackedOff,
            PredictorError::ModelError(error) => PathError::ModelError(error),
        }
    }
}

impl From<StepFailure> for PathError {
    fn from(failure: StepFailure) -> Self {
        match failure {
            StepFailure::CorrectorStuck => PathError::CorrectorStuck,
            StepFailure::TooManyTries => PathError::TooManyTries,
        }
    }
}

/// Follow the path of `homotopy` from `(x, 0)` to `λ = 1`.
///
/// `homotopy` has `λ` as its last unknown, see [`crate::NewtonHomotopy`] and
/// [`crate::FixpointHomotopy`]. `x` holds the start point on entry and the
/// solution of `H(x, 1) = 0` on success, it is unchanged otherwise.
/// `x_scaling` scales the unknowns, `λ` last. The path starts in `direction`.
///
/// Each step computes the tangent at the latest point on the path, predicts
/// along it, corrects back to the path and lets the step control accept the
/// step or retry it with a smaller step size. Once a predictor reaches
/// `λ = 1` exactly, the corrector keeps `λ` fixed and uses the residual
/// tolerance `f_tol` of Newton's method, for the rest of the path like the C
/// runtime. Returns the number of accepted steps. From `homotopyAlgorithm`.
///
/// # Errors
///
/// Returns a [`PathError`] if a step fails or the path leaves its limits.
///
/// # Panics
///
/// Panics if `x` and `x_scaling` don't match the homotopy.
pub fn follow_path(
    homotopy: &mut impl ModelFunctions,
    x_scaling: &[f64],
    direction: Direction,
    f_tol: f64,
    settings: &Settings,
    x: &mut [f64],
) -> Result<usize, PathError> {
    let n = homotopy.info().num_residuals;
    let m = homotopy.info().num_unknowns;
    assert_eq!(m, n + 1, "Homotopy needs λ as last unknown");
    assert_eq!(x.len(), n, "Start point has wrong length");
    assert_eq!(x_scaling.len(), m, "Scaling has wrong length");

    let mut y0: Vec<f64> = x.iter().copied().chain([0.0]).collect();
    let mut y1 = vec![0.0; m];
    let mut yt = vec![0.0; m];
    let mut t = vec![0.0; m];
    let mut previous = vec![0.0; m];
    previous[n] = direction.sign();
    let mut h = vec![0.0; n];
    let mut jacobian = Mat::zeros(n, m);
    let mut work = CorrectorWorkspace::new(x_scaling);

    homotopy.eval(&y0, &mut h)?;
    let mut tau = settings.tau_start;
    let mut tries = 0;
    let mut steps = 0;
    let mut fixed = n;
    let mut h_eps = settings.h_eps;

    while y0[n] < 1.0 {
        if y0[n] < -1.0 {
            return Err(PathError::LambdaBelowMinusOne);
        }
        if steps >= settings.max_lambda_steps {
            return Err(PathError::TooManySteps);
        }

        if tries == 0 {
            fixed = tangent(homotopy, &y0, &h, x_scaling, &mut jacobian, &mut t)?;
            let first_step = steps == 0;
            orient(
                &mut t,
                &previous,
                first_step && direction == Direction::Negative,
            );
            tau = limit_to_end(tau, y0[n], t[n]);
        }

        tau = predictor(homotopy, &y0, &t, tau, settings, &mut y1, &mut h)?;
        yt.copy_from_slice(&y1);
        if yt[n] == 1.0 {
            fixed = n;
            h_eps = f_tol;
        }

        let correction = corrector(
            homotopy,
            fixed,
            h_eps,
            settings.max_corrector_steps,
            &mut work,
            &mut y1,
            &mut h,
        );
        let (outcome, bend) = match correction {
            Ok(Correction::Converged) => (Outcome::Converged, bend(&y0, &yt, &y1)),
            Ok(Correction::NotConverged) => (Outcome::NotConverged, bend(&y0, &yt, &y1)),
            Err(CorrectorError::ModelError | CorrectorError::Singular) => {
                (Outcome::Failed, f64::NAN)
            }
        };

        match control(outcome, bend, tau, tries, settings) {
            Control::Accept { tau: next_tau } => {
                trace!(target: "NLS_HOMOTOPY", "Homotopy step {steps}: λ {}, τ {tau}, bend {bend}", y1[n]);
                tau = next_tau;
                tries = 0;
                steps += 1;
                y0.copy_from_slice(&y1);
                previous.copy_from_slice(&t);
            }
            Control::Retry {
                tau: next_tau,
                tries: next_tries,
            } => {
                tau = next_tau;
                tries = next_tries;
            }
            Control::Fail(failure) => return Err(failure.into()),
        }
    }

    x.copy_from_slice(&y0[..n]);
    Ok(steps)
}

#[cfg(test)]
mod path_follower_tests {
    use super::*;
    use crate::NewtonHomotopy;
    use nls_model::ClosureModel;

    /// Follow the Newton homotopy of the one-dimensional `f` from `x0`
    fn follow(
        f: impl FnMut(&[f64], &mut [f64]) -> Result<(), EvalError>,
        x0: f64,
        direction: Direction,
        settings: &Settings,
    ) -> (Result<usize, PathError>, f64) {
        let mut model = ClosureModel::new(1, 1, f);
        let mut fx0 = [0.0];
        model.eval(&[x0], &mut fx0).unwrap();
        let mut homotopy = NewtonHomotopy::new(&mut model, &fx0);
        let mut x = [x0];
        let x_scaling = [x0.abs().max(1.0), 1.0];

        let result = follow_path(
            &mut homotopy,
            &x_scaling,
            direction,
            1e-12,
            settings,
            &mut x,
        );
        (result, x[0])
    }

    #[test]
    /// The Newton homotopy of `x − 2` from 0 is the line `x = 2λ`
    fn linear_test() {
        let (result, x) = follow(
            |x, f| {
                f[0] = x[0] - 2.0;
                Ok(())
            },
            0.0,
            Direction::Positive,
            &Settings::new(1),
        );

        assert!(result.is_ok(), "got {result:?}");
        assert!((x - 2.0).abs() < 1e-5, "{x}");
    }

    #[test]
    /// Newton's method diverges for `arctan(x) = 0` from `x = 2`, the path of
    /// the Newton homotopy leads to the root 0
    fn arctan_test() {
        let (result, x) = follow(
            |x, f| {
                f[0] = x[0].atan();
                Ok(())
            },
            2.0,
            Direction::Positive,
            &Settings::new(1),
        );

        assert!(result.is_ok(), "got {result:?}");
        assert!(x.abs() < 1e-5, "{x}");
    }

    /// `F(x) = x³ − 3x − 3` with its only real root at `x ≈ 2.1038`. The path
    /// `F(x) = −3·(1 − λ)` of its Newton homotopy from 0 turns at `x = ±1`.
    fn cubic(x: &[f64], f: &mut [f64]) -> Result<(), EvalError> {
        f[0] = x[0].powi(3) - 3.0 * x[0] - 3.0;
        Ok(())
    }

    #[test]
    /// In the negative direction the path goes down to `λ = −2/3` at `x = 1`,
    /// turns and reaches the root at `λ = 1`
    fn turning_point_test() {
        let (result, x) = follow(cubic, 0.0, Direction::Negative, &Settings::new(1));

        assert!(result.is_ok(), "got {result:?}");
        assert!((x - 2.103_803_402_735_535).abs() < 1e-5, "{x}");
    }

    #[test]
    /// In the positive direction the path turns at `λ = 2/3` and falls below
    /// `λ = −1`
    fn lambda_below_minus_one_test() {
        let (result, _) = follow(cubic, 0.0, Direction::Positive, &Settings::new(1));

        assert!(
            matches!(result, Err(PathError::LambdaBelowMinusOne)),
            "got {result:?}"
        );
    }

    #[test]
    /// A path needs more than one step to reach `λ = 1`
    fn too_many_steps_test() {
        let settings = Settings::new(1).with_max_lambda_steps(1);
        let (result, x) = follow(
            |x, f| {
                f[0] = x[0].atan();
                Ok(())
            },
            2.0,
            Direction::Positive,
            &settings,
        );

        assert!(
            matches!(result, Err(PathError::TooManySteps)),
            "got {result:?}"
        );
        assert_eq!(x, 2.0, "x is unchanged");
    }

    #[test]
    /// A homotopy that can't be evaluated at the start fails at once
    fn model_error_test() {
        let mut model = ClosureModel::new(1, 1, |x: &[f64], f: &mut [f64]| {
            if x[0] < 1.0 {
                return Err(EvalError::Assert);
            }
            f[0] = x[0];
            Ok(())
        });
        let mut homotopy = NewtonHomotopy::new(&mut model, &[0.0]);
        let mut x = [0.0];

        let result = follow_path(
            &mut homotopy,
            &[1.0, 1.0],
            Direction::Positive,
            1e-12,
            &Settings::new(1),
            &mut x,
        );

        assert!(
            matches!(result, Err(PathError::ModelError(EvalError::Assert))),
            "got {result:?}"
        );
    }
}
