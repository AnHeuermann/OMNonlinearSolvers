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

use crate::NewtonWorkspace;
use crate::cubic_interpolation::{TrialPoint, cubic_interpolation};

use nls_model::EvalError;

///
///
/// 1. Back off on a failed evaluation. Try x + λ·dy starting at λ = 1. While the model throws an assert, multiply λ by 0.655. Give up once λ < 1e-4. A constraint violation aborts at once instead of backing off.
/// 2. Sufficient decrease test (Armijo, α = 0.1). The slope is g = −2·err. If the new error err1 > err + α·λ·g, and the scaled error fails the same test, and both errors are above 1e-12, the step is too long.
/// 3. Quadratic interpolation. Compute λ₂ = max(−λ²·g / (2·(err1 − err − λ·g)), 1e-4) and evaluate at λ₂.
/// 4. Cubic interpolation. If the test still fails, fit a cubic through both trial points, take its minimum as λ, evaluate once more and accept that step without testing it again.
pub(crate) fn line_search(workspace: &NewtonWorkspace) {
    todo!()
}

/// Armijo inequality err_1 > err_0 + α·λ₁·g with `α = 1e-1`.
fn armijo_inequality(err_0: f64, err_1: f64, lambda_1: f64, g_0: f64) -> bool {
    const ALPHA: f64 = 1e-1;
    err_1 > err_0 + ALPHA * lambda_1 * g_0
}

/// Smallest step factor the line search takes
const LAMBDA_MIN: f64 = 1e-4;

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
}
