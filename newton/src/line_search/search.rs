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

// TODO: Remove once the port is complete
#![allow(dead_code, unused_imports, unused_variables)]

use super::TrialPoint;
use super::cubic_interpolation::cubic_interpolation;

use nls_model::EvalError;

/// Smallest step factor the line search takes
const LAMBDA_MIN: f64 = 1e-4;

/// Evaluation of the residual `F(x)` along the search direction
///
/// Writes `F(x)` into the second slice, or returns the model's [`EvalError`].
/// Every closure with this signature implements it, for example
/// `|x, f| model.eval(x, f)`.
// A trait alias once `trait_alias` is stable
pub(crate) trait EvalFunction: FnMut(&[f64], &mut [f64]) -> Result<(), EvalError> {}
impl<T: FnMut(&[f64], &mut [f64]) -> Result<(), EvalError>> EvalFunction for T {}

/// Why the line search found no step
#[derive(Debug)]
pub(crate) enum LineSearchError {
    /// The model asserted for every step factor down to [`LAMBDA_MIN`]
    BackedOffBelowLambdaMin,
    /// The model failed at a damped step, where the line search doesn't back
    /// off
    ModelError(EvalError),
}

impl From<EvalError> for LineSearchError {
    fn from(value: EvalError) -> Self {
        LineSearchError::ModelError(value)
    }
}

/// Line search
///
/// Iterative approach to find a local minimum `x1` of the error function with
/// the residual `f(x)`.
///
/// `x`, the current iterate `dy`, the Newton step, already scaled back to the
/// original variables by `solve` `res_scaling` `err` and `err_scaled` from the
/// last iteration
/// Writes accepted trial point into `x1`  and residual `F(x1)` into `f`
pub(crate) fn line_search(
    eval: &mut impl EvalFunction,
    starting_x: &[f64],
    direction_dy: &[f64],
    phi_0: f64,
    g: f64,
) -> Result<(), LineSearchError> {
    // lambda = back_off_failed_eval(eval, starting_x, direction_dy)?;

    // 2. Sufficient decrease test (Armijo, α = 0.1). The slope is g = −2·err. If
    //    the new error err1 > err + α·λ·g, and the scaled error fails the same
    //    test, and both errors are above 1e-12, the step is too long.
    // 3. Quadratic interpolation. Compute λ₂ = max(−λ²·g / (2·(err1 − err − λ·g)),
    //    1e-4) and evaluate at λ₂.
    // 4. Cubic interpolation. If the test still fails, fit a cubic through both
    //    trial points, take its minimum as λ, evaluate once more and accept that
    //    step without testing it again.

    todo!()
}

/// Evaluate the trial point `x1 = x + λ·dy` and write `F(x1)` into `f`
fn eval_at(
    eval: &mut impl EvalFunction,
    x: &[f64],
    dy: &[f64],
    lambda: f64,
    x1: &mut [f64],
    f: &mut [f64],
) -> Result<(), EvalError> {
    for ((x1, x), dy) in x1.iter_mut().zip(x).zip(dy) {
        *x1 = x + lambda * dy;
    }
    eval(x1, f)
}

/// Factor by which `λ` shrinks after a failed evaluation.
const LAMBDA_SCALE_DOWN: f64 = 0.655;

/// Back off from the full step `x + dy` until the model can be evaluated.
///
/// Tries `x1 = x + λ·dy` starting at `λ = 1`. While `eval` reports an
/// [`EvalError::Assert`], `λ` is multiplied by [`LAMBDA_SCALE_DOWN`].
///
/// Returns the first `λ` at which `eval` succeeds. `x1` then holds the trial
/// point and `f` its residual `F(x1)`.
///
/// Unlike the C runtime, which evaluates once more at a `λ` below
/// [`LAMBDA_MIN`] and then gives up anyway, it stops before that evaluation.
///
/// # Errors
///
/// * [`LineSearchError::ModelError`] at once if `eval` reports a violated
///   constraint, without backing off
/// * [`LineSearchError::BackedOffBelowLambdaMin`] if `eval` still fails at
///   the smallest `λ` not below [`LAMBDA_MIN`]
///
/// After an error `x1` and `f` are unspecified. `x` is never changed.
fn back_off_failed_eval(
    eval: &mut impl EvalFunction,
    x: &[f64],
    dy: &[f64],
    x1: &mut [f64],
    f: &mut [f64],
) -> Result<f64, LineSearchError> {
    debug_assert_eq!(x1.len(), x.len());
    debug_assert_eq!(x1.len(), dy.len());

    let mut lambda = 1.0;
    while lambda >= LAMBDA_MIN {
        match eval_at(eval, x, dy, lambda, x1, f) {
            Ok(()) => return Ok(lambda),
            Err(EvalError::Assert) => lambda *= LAMBDA_SCALE_DOWN,
            Err(e @ EvalError::ConstraintViolated) => return Err(e.into()),
        }
    }
    Err(LineSearchError::BackedOffBelowLambdaMin)
}

/// Factor `α` of the sufficient decrease in the Armijo test
pub(crate) const ALPHA: f64 = 1e-1;

/// Whether the trial step `λ₁` fails the Armijo test: its error `err_1`
/// didn't decrease by at least the fraction [`ALPHA`] of the slope,
/// `err_1 > err_0 + α·λ₁·g_0`.
fn insufficient_decrease(err_0: f64, err_1: f64, lambda_1: f64, g_0: f64) -> bool {
    err_1 > err_0 + ALPHA * lambda_1 * g_0
}

fn quadratic_step() {}

/// Step factor from the cubic interpolation through the trial points `p1` and
/// `p2`, where `p1` is the first trial point of the line search.
///
/// Like the C runtime it takes half the first step if the cubic has no local
/// minimum, which also covers a non-finite minimum. The result is at least
/// [`LAMBDA_MIN`].
fn cubic_step(err: f64, g: f64, p1: &TrialPoint, p2: &TrialPoint) -> f64 {
    cubic_interpolation(err, g, p1, p2)
        .filter(|lambda| lambda.is_finite())
        .unwrap_or(0.5 * p1.lambda)
        .max(LAMBDA_MIN)
}

#[cfg(test)]
mod line_search_tests {
    use std::f32::consts::E;

    use super::*;

    #[test]
    /// `φ(λ) = 1 − 2λ` has neither a cubic nor a quadratic term, so the
    /// minimum of the fit is `+∞`. The step is half the first one instead.
    fn cubic_step_not_finite_test() {
        let err = 1.0;
        let g = -2.0 * err;
        let p1 = TrialPoint {
            lambda: 1.0,
            error: -1.0,
        };
        let p2 = TrialPoint {
            lambda: 0.5,
            error: 0.0,
        };
        assert_eq!(cubic_interpolation(err, g, &p1, &p2), Some(f64::INFINITY));

        assert_eq!(cubic_step(err, g, &p1, &p2), 0.5);
    }

    /// Result, number of evaluations, `x1` and `f` of a back off
    type BackOff = (Result<f64, LineSearchError>, usize, [f64; 2], [f64; 2]);

    /// Back off from `x = [1, 2]` along `dy = [2, -4]` with the model
    /// `F(x) = x`. `fail(call)` is the error of the evaluation number `call`,
    /// counted from 1, or `None` if it succeeds.
    fn back_off(fail: impl Fn(usize) -> Option<EvalError>) -> BackOff {
        let x = [1.0, 2.0];
        let dy = [2.0, -4.0];
        let mut x1 = [0.0; 2];
        let mut f = [0.0; 2];
        let mut calls = 0;
        let mut eval = |x: &[f64], f: &mut [f64]| {
            calls += 1;
            f.copy_from_slice(x);
            fail(calls).map_or(Ok(()), Err)
        };

        let result = back_off_failed_eval(&mut eval, &x, &dy, &mut x1, &mut f);
        (result, calls, x1, f)
    }

    /// `λ` after `k` back offs, shrunk in the same order as the line search
    fn shrunk(k: usize) -> f64 {
        (0..k).fold(1.0, |lambda, _| lambda * LAMBDA_SCALE_DOWN)
    }

    #[test]
    /// The full step is accepted after one evaluation, `x1 = x + dy`
    fn back_off_full_step_test() {
        let (result, calls, x1, f) = back_off(|_| None);

        assert!(matches!(result, Ok(1.0)), "got {result:?}");
        assert_eq!(calls, 1);
        assert_eq!(x1, [3.0, -2.0]);
        assert_eq!(f, x1, "f must be F(x1)");
    }

    #[test]
    /// Three asserts shrink `λ` three times, the fourth evaluation succeeds
    fn back_off_asserts_test() {
        let (result, calls, x1, f) = back_off(|call| (call <= 3).then_some(EvalError::Assert));

        let lambda = shrunk(3);
        assert!(matches!(result, Ok(l) if l == lambda), "got {result:?}");
        assert_eq!(calls, 4);
        assert_eq!(x1, [1.0 + lambda * 2.0, 2.0 + lambda * -4.0]);
        assert_eq!(f, x1, "f must be F(x1)");
    }

    #[test]
    /// A success at the smallest `λ = 0.655²¹ ≥ λ_min` is accepted
    fn back_off_smallest_lambda_test() {
        let (result, calls, _, _) = back_off(|call| (call <= 21).then_some(EvalError::Assert));

        let lambda = shrunk(21);
        assert!(lambda >= LAMBDA_MIN && lambda * LAMBDA_SCALE_DOWN < LAMBDA_MIN);
        assert!(matches!(result, Ok(l) if l == lambda), "got {result:?}");
        assert_eq!(calls, 22);
    }

    #[test]
    /// A model that always asserts is tried at the 22 step factors from 1 to
    /// `0.655²¹`, then the line search gives up
    fn back_off_below_lambda_min_test() {
        let (result, calls, _, _) = back_off(|_| Some(EvalError::Assert));

        assert!(
            matches!(result, Err(LineSearchError::BackedOffBelowLambdaMin)),
            "got {result:?}"
        );
        assert_eq!(calls, 22);
    }

    #[test]
    /// A violated constraint aborts at once, without backing off
    fn back_off_constraint_test() {
        let (result, calls, _, _) = back_off(|_| Some(EvalError::ConstraintViolated));

        assert!(
            matches!(
                result,
                Err(LineSearchError::ModelError(EvalError::ConstraintViolated))
            ),
            "got {result:?}"
        );
        assert_eq!(calls, 1);
    }

    #[test]
    /// A violated constraint after some asserts still aborts at once
    fn back_off_constraint_after_asserts_test() {
        let (result, calls, _, _) = back_off(|call| match call {
            1 | 2 => Some(EvalError::Assert),
            _ => Some(EvalError::ConstraintViolated),
        });

        assert!(
            matches!(
                result,
                Err(LineSearchError::ModelError(EvalError::ConstraintViolated))
            ),
            "got {result:?}"
        );
        assert_eq!(calls, 3);
    }
}
