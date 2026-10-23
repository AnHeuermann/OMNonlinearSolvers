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

//! Tests of the strategy that combines Newton's method and the homotopies.

use homotopy_core::settings::Settings;
use homotopy_core::{Direction, HomotopyError, HomotopyRun, HomotopySolver, Solved};
use newton::settings::Settings as NewtonSettings;
use nls_model::{ClosureModel, EvalError, ModelFunctions};

/// A solver with the default settings for `size` unknowns
fn solver(size: usize) -> HomotopySolver {
    HomotopySolver::new(size, NewtonSettings::new(size), Settings::new(size))
}

#[test]
/// A start point that solves the system is taken as it is
fn start_point_test() {
    let mut model = ClosureModel::new(1, 1, |x: &[f64], f: &mut [f64]| {
        f[0] = x[0] - 3.0;
        Ok(())
    });
    let mut x = [3.0];

    let result = solver(1).solve(&mut model, &mut x);

    assert!(matches!(result, Ok(Solved::StartPoint)), "got {result:?}");
    assert_eq!(x, [3.0]);
}

#[test]
/// A system Newton's method solves needs no homotopy
fn newton_test() {
    let mut model = ClosureModel::new(2, 2, |x: &[f64], f: &mut [f64]| {
        f[0] = x[0] + 2.0 * x[1] - 5.0;
        f[1] = 3.0 * x[0] - x[1] - 1.0;
        Ok(())
    });
    let mut x = [0.0, 0.0];

    let result = solver(2).solve(&mut model, &mut x);

    assert!(matches!(result, Ok(Solved::Newton)), "got {result:?}");
    assert!(
        (x[0] - 1.0).abs() < 1e-10 && (x[1] - 2.0).abs() < 1e-10,
        "{x:?}"
    );
}

#[test]
/// `F(x) = x³ − 2x + 2` from 0: Newton's method circles around the local
/// minimum of `|F|` at `x = √(2/3)` until it runs out of iterations. The
/// path of the Newton homotopy turns
/// there and falls below `λ = −1`, in the negative direction it turns at
/// `x = −√(2/3)` and reaches the root, which Newton then refines.
fn homotopy_test() {
    let mut model = ClosureModel::new(1, 1, |x: &[f64], f: &mut [f64]| {
        f[0] = x[0].powi(3) - 2.0 * x[0] + 2.0;
        Ok(())
    });
    let mut x = [0.0];

    let result = solver(1).solve(&mut model, &mut x);

    assert!(
        matches!(result, Ok(Solved::Homotopy(HomotopyRun::NewtonNegative))),
        "got {result:?}"
    );
    assert!((x[0] + 1.769_292_354_238_631).abs() < 1e-10, "{x:?}");
}

#[test]
/// `x² + 1 = 0` has no real root: everything fails and `x` is unchanged
fn not_solved_test() {
    let mut model = ClosureModel::new(1, 1, |x: &[f64], f: &mut [f64]| {
        f[0] = x[0] * x[0] + 1.0;
        Ok(())
    });
    let mut x = [0.5];

    let result = solver(1).solve(&mut model, &mut x);

    assert!(
        matches!(result, Err(HomotopyError::NotSolved)),
        "got {result:?}"
    );
    assert_eq!(x, [0.5]);
}

/// `F = (x₀ − 1, x₁ − 2)`, which asserts at `x₁ = 0`
fn asserts_at_zero(x: &[f64], f: &mut [f64]) -> Result<(), EvalError> {
    if x[1] == 0.0 {
        return Err(EvalError::Assert);
    }
    f[0] = x[0] - 1.0;
    f[1] = x[1] - 2.0;
    Ok(())
}

#[test]
/// A model that asserts at the start point starts from a moved one: the
/// second unknown moves by `1/2·1%` of its scaling
fn moved_start_point_test() {
    let mut model = ClosureModel::new(2, 2, asserts_at_zero);
    let mut x = [0.0, 0.0];

    let result = solver(2).solve(&mut model, &mut x);

    assert!(matches!(result, Ok(Solved::Newton)), "got {result:?}");
    assert!(
        (x[0] - 1.0).abs() < 1e-10 && (x[1] - 2.0).abs() < 1e-10,
        "{x:?}"
    );
}

#[test]
/// A casual tearing set doesn't move its start point
fn casual_tearing_set_test() {
    let mut model = ClosureModel::new(2, 2, asserts_at_zero).with_casual_tearing_set();
    let mut x = [0.0, 0.0];

    let result = solver(2).solve(&mut model, &mut x);

    assert!(
        matches!(result, Err(HomotopyError::ModelError(EvalError::Assert))),
        "got {result:?}"
    );
    assert_eq!(x, [0.0, 0.0]);
}

#[test]
/// A violated constraint at the start point gives up at once
fn constraint_violated_test() {
    let mut model = ClosureModel::new(1, 1, |_: &[f64], _: &mut [f64]| {
        Err(EvalError::ConstraintViolated)
    })
    .with_casual_tearing_set();
    let mut x = [1.0];

    let result = solver(1).solve(&mut model, &mut x);

    assert!(
        matches!(
            result,
            Err(HomotopyError::ModelError(EvalError::ConstraintViolated))
        ),
        "got {result:?}"
    );
}

#[test]
/// A model that can't be evaluated anywhere fails after the moved start points
fn start_point_model_error_test() {
    let mut model = ClosureModel::new(1, 1, |_: &[f64], _: &mut [f64]| Err(EvalError::Assert));
    let mut x = [1.0];

    let result = solver(1).solve(&mut model, &mut x);

    assert!(
        matches!(result, Err(HomotopyError::ModelError(EvalError::Assert))),
        "got {result:?}"
    );
}

/// The model's own homotopy `H(x, λ) = x − 2λ` with the solution 2 at `λ = 1`
fn own_homotopy() -> impl ModelFunctions {
    ClosureModel::new(1, 2, |y: &[f64], h: &mut [f64]| {
        h[0] = y[0] - 2.0 * y[1];
        Ok(())
    })
}

#[test]
/// The end of the path of the model's own homotopy is the solution
fn own_homotopy_test() {
    let mut x = [0.0];

    let result = solver(1).solve_own_homotopy(&mut own_homotopy(), &mut x, Direction::Positive);

    assert!(matches!(result, Ok(Direction::Positive)), "got {result:?}");
    assert!((x[0] - 2.0).abs() < 1e-5, "{x:?}");
}

#[test]
/// In the negative direction the path falls below `λ = −1`, the positive
/// direction is tried next
fn own_homotopy_opposite_test() {
    let mut x = [0.0];

    let result = solver(1).solve_own_homotopy(&mut own_homotopy(), &mut x, Direction::Negative);

    assert!(matches!(result, Ok(Direction::Positive)), "got {result:?}");
    assert!((x[0] - 2.0).abs() < 1e-5, "{x:?}");
}
