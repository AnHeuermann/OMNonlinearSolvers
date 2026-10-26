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

//! Tests of the C interface, with callbacks as a C host would pass them.

use homotopy_capi::*;
use std::ffi::{c_int, c_void};
use std::ptr::{null, null_mut};

/// `F(x) = x³ − 2x + 2`, which needs a homotopy from 0
unsafe extern "C" fn cubic(_: *mut c_void, x: *const f64, f: *mut f64) -> c_int {
    unsafe {
        let x = *x;
        *f = x * x * x - 2.0 * x + 2.0;
    }
    HOM_EVAL_OK
}

/// `J(x) = 3x² − 2`, counting its calls in `user`
unsafe extern "C" fn cubic_jacobian(
    user: *mut c_void,
    x: *const f64,
    _: *const f64,
    jac: *mut f64,
) -> c_int {
    unsafe {
        *(user as *mut usize) += 1;
        *jac = 3.0 * *x * *x - 2.0;
    }
    HOM_EVAL_OK
}

/// `F(x) = (x₀ + 2x₁ − 5, 3x₀ − x₁ − 1)`, which Newton solves at once
unsafe extern "C" fn linear(_: *mut c_void, x: *const f64, f: *mut f64) -> c_int {
    unsafe {
        let x = std::slice::from_raw_parts(x, 2);
        let f = std::slice::from_raw_parts_mut(f, 2);
        f[0] = x[0] + 2.0 * x[1] - 5.0;
        f[1] = 3.0 * x[0] - x[1] - 1.0;
    }
    HOM_EVAL_OK
}

/// Asserts everywhere
unsafe extern "C" fn asserts(_: *mut c_void, _: *const f64, _: *mut f64) -> c_int {
    HOM_EVAL_ASSERT
}

/// Violates its constraint everywhere
unsafe extern "C" fn violates(_: *mut c_void, _: *const f64, _: *mut f64) -> c_int {
    HOM_EVAL_CONSTRAINT_VIOLATED
}

/// The model's own homotopy `H(x, λ) = x − 2λ`
unsafe extern "C" fn own_homotopy(_: *mut c_void, y: *const f64, h: *mut f64) -> c_int {
    unsafe {
        let y = std::slice::from_raw_parts(y, 2);
        *h = y[0] - 2.0 * y[1];
    }
    HOM_EVAL_OK
}

/// A model of `size` unknowns with the `residual` callback only
fn model(size: usize, residual: HomResidual) -> HomModel {
    HomModel {
        user: null_mut(),
        num_residuals: size,
        num_unknowns: size,
        nominal: null(),
        min: null(),
        max: null(),
        casual_tearing_set: false,
        residual,
        jacobian: None,
    }
}

/// Status, solved by and `x` of solving `model` from `x` with the default
/// settings
fn solve(model: &HomModel, x: &mut [f64]) -> (HomStatus, Option<HomSolvedBy>) {
    unsafe {
        let solver = hom_solver_new(x.len(), null());
        assert!(!solver.is_null());
        let mut solved_by = HomSolvedBy::StartPoint;
        let status = hom_solve(solver, model, x.as_mut_ptr(), &mut solved_by);
        hom_solver_free(solver);
        (status, (status == HomStatus::Solved).then_some(solved_by))
    }
}

#[test]
/// A system Newton's method solves
fn newton_test() {
    let mut x = [0.0, 0.0];

    let (status, solved_by) = solve(&model(2, Some(linear)), &mut x);

    assert_eq!(
        (status, solved_by),
        (HomStatus::Solved, Some(HomSolvedBy::Newton))
    );
    assert!(
        (x[0] - 1.0).abs() < 1e-10 && (x[1] - 2.0).abs() < 1e-10,
        "{x:?}"
    );
}

#[test]
/// A system that needs the Newton homotopy in the negative direction
fn homotopy_test() {
    let mut x = [0.0];

    let (status, solved_by) = solve(&model(1, Some(cubic)), &mut x);

    assert_eq!(
        (status, solved_by),
        (HomStatus::Solved, Some(HomSolvedBy::HomotopyNewtonNegative))
    );
    assert!((x[0] + 1.769_292_354_238_631).abs() < 1e-10, "{x:?}");
}

#[test]
/// An analytical Jacobian is used instead of forward differences
fn jacobian_test() {
    let mut calls = 0_usize;
    let model = HomModel {
        user: (&raw mut calls).cast(),
        jacobian: Some(cubic_jacobian),
        ..model(1, Some(cubic))
    };
    let mut x = [0.0];

    let (status, _) = solve(&model, &mut x);

    assert_eq!(status, HomStatus::Solved);
    assert!((x[0] + 1.769_292_354_238_631).abs() < 1e-10, "{x:?}");
    assert!(calls > 0, "the Jacobian callback was used");
}

#[test]
/// A model that can't be evaluated is a model error, `x` is unchanged
fn model_error_test() {
    let mut x = [1.0];

    let (status, _) = solve(&model(1, Some(asserts)), &mut x);

    assert_eq!(status, HomStatus::ModelError);
    assert_eq!(x, [1.0]);
}

#[test]
/// A violated constraint of a casual tearing set is a model error
fn constraint_violated_test() {
    let model = HomModel {
        casual_tearing_set: true,
        ..model(1, Some(violates))
    };
    let mut x = [1.0];

    let (status, _) = solve(&model, &mut x);

    assert_eq!(status, HomStatus::ModelError);
}

#[test]
/// The model's own homotopy with `λ` as additional unknown
fn own_homotopy_test() {
    let model = HomModel {
        num_unknowns: 2,
        ..model(1, Some(own_homotopy))
    };
    let mut x = [0.0];
    let mut solved_by = HomSolvedBy::StartPoint;

    let status = unsafe {
        let solver = hom_solver_new(1, null());
        let status = hom_solve_own_homotopy(solver, &model, x.as_mut_ptr(), true, &mut solved_by);
        hom_solver_free(solver);
        status
    };

    assert_eq!(status, HomStatus::Solved);
    assert_eq!(
        solved_by,
        HomSolvedBy::OwnHomotopyPositive,
        "the negative start falls below λ = −1"
    );
    assert!((x[0] - 2.0).abs() < 1e-5, "{x:?}");
}

/// `F(x) = (x₀ − 1, x₁ − 2)`, which asserts at `x₁ = 0`. Records `x₁` of the
/// first evaluation without an assert in `user`, a NaN-initialized `f64`.
unsafe extern "C" fn record_first_x1(user: *mut c_void, x: *const f64, f: *mut f64) -> c_int {
    unsafe {
        let x = std::slice::from_raw_parts(x, 2);
        if x[1] == 0.0 {
            return HOM_EVAL_ASSERT;
        }
        let first = &mut *(user as *mut f64);
        if first.is_nan() {
            *first = x[1];
        }
        let f = std::slice::from_raw_parts_mut(f, 2);
        f[0] = x[0] - 1.0;
        f[1] = x[1] - 2.0;
    }
    HOM_EVAL_OK
}

/// `x₁` of the first start point the model can be evaluated at, from `(0, 0)`
/// with the `nominal` values
fn first_start_point(nominal: *const f64) -> f64 {
    let mut first = f64::NAN;
    let model = HomModel {
        user: (&raw mut first).cast(),
        nominal,
        ..model(2, Some(record_first_x1))
    };
    let mut x = [0.0, 0.0];

    let (status, _) = solve(&model, &mut x);
    assert_eq!(status, HomStatus::Solved);
    first
}

#[test]
/// The nominal values scale the unknowns: the start point moves by
/// `scaling · i/n · 1%`, with the scaling `max(nominal, |x|)`
fn nominal_test() {
    let nominal = [1.0, 100.0];

    assert_eq!(first_start_point(nominal.as_ptr()), 100.0 * 0.5 * 0.01);
    assert_eq!(
        first_start_point(null()),
        0.5 * 0.01,
        "nominal 1 without values"
    );
}

/// `F(x) = x − 1`, which asserts above 1
unsafe extern "C" fn asserts_above_one(_: *mut c_void, x: *const f64, f: *mut f64) -> c_int {
    unsafe {
        if *x > 1.0 {
            return HOM_EVAL_ASSERT;
        }
        *f = *x - 1.0;
    }
    HOM_EVAL_OK
}

#[test]
/// Just below the upper bound the forward differences step backwards, without
/// the bound they step to where the model asserts
fn max_test() {
    let max = [1.0];
    let with_max = HomModel {
        max: max.as_ptr(),
        ..model(1, Some(asserts_above_one))
    };
    let mut x = [1.0 - 1e-10];
    let (status, _) = solve(&with_max, &mut x);
    assert_eq!(status, HomStatus::Solved);
    assert!((x[0] - 1.0).abs() < 1e-12, "{x:?}");

    let mut x = [1.0 - 1e-10];
    let (status, _) = solve(&model(1, Some(asserts_above_one)), &mut x);
    assert_eq!(status, HomStatus::ModelError);
}

#[test]
/// The lower bounds don't restrict the solution, like in the C runtime
fn min_test() {
    let min = [5.0];
    let model = HomModel {
        min: min.as_ptr(),
        ..model(1, Some(cubic))
    };
    let mut x = [0.0];

    let (status, _) = solve(&model, &mut x);

    assert_eq!(status, HomStatus::Solved);
    assert!((x[0] + 1.769_292_354_238_631).abs() < 1e-10, "{x:?}");
}

#[test]
/// `NULL` pointers, a missing residual and sizes that don't match are
/// rejected
fn invalid_argument_test() {
    let mut x = [0.0];
    unsafe {
        let solver = hom_solver_new(1, null());
        let cubic = model(1, Some(cubic));
        assert_eq!(
            hom_solve(null_mut(), &cubic, x.as_mut_ptr(), null_mut()),
            HomStatus::InvalidArgument
        );
        assert_eq!(
            hom_solve(solver, null(), x.as_mut_ptr(), null_mut()),
            HomStatus::InvalidArgument
        );
        assert_eq!(
            hom_solve(solver, &cubic, null_mut(), null_mut()),
            HomStatus::InvalidArgument
        );
        assert_eq!(
            hom_solve(solver, &model(1, None), x.as_mut_ptr(), null_mut()),
            HomStatus::InvalidArgument
        );
        assert_eq!(
            hom_solve(solver, &model(2, Some(linear)), x.as_mut_ptr(), null_mut()),
            HomStatus::InvalidArgument
        );
        assert_eq!(
            hom_solve_own_homotopy(solver, &cubic, x.as_mut_ptr(), false, null_mut()),
            HomStatus::InvalidArgument,
            "the own homotopy needs λ"
        );
        hom_solver_free(solver);
    }
}

#[test]
/// Settings of 0 keep the defaults, invalid settings and size 0 give `NULL`
fn settings_test() {
    unsafe {
        let defaults = HomSettings::default();
        let solver = hom_solver_new(1, &defaults);
        assert!(!solver.is_null());
        hom_solver_free(solver);

        let custom = HomSettings {
            tau_start: 0.1,
            max_tries: 5,
            ..HomSettings::default()
        };
        let solver = hom_solver_new(1, &custom);
        assert!(!solver.is_null());
        hom_solver_free(solver);

        let invalid = HomSettings {
            tau_decreasing_factor: 0.5,
            ..HomSettings::default()
        };
        assert!(hom_solver_new(1, &invalid).is_null(), "a factor below 1");
        assert!(hom_solver_new(0, null()).is_null(), "size 0");
    }
}

#[test]
/// The generated header compiles as C99, if a C compiler is installed
fn header_compiles_test() {
    let crate_dir = env!("CARGO_MANIFEST_DIR");
    let program =
        std::env::temp_dir().join(format!("homotopy_capi_header_{}.c", std::process::id()));
    std::fs::write(
        &program,
        "#include \"homotopy_capi.h\"\n\
         static int residual(void *user, const double *x, double *f) { (void)user; *f = *x; return HOM_EVAL_OK; }\n\
         int main(void) {\n\
           HomSettings settings = {0};\n\
           HomModel model = {0};\n\
           model.residual = residual;\n\
           HomSolvedBy solved_by;\n\
           (void)settings; (void)model; (void)solved_by;\n\
           return 0;\n\
         }\n",
    )
    .unwrap();

    let compiler = std::env::var("CC").unwrap_or_else(|_| "cc".into());
    let output = std::process::Command::new(&compiler)
        .args([
            "-std=c99",
            "-Wall",
            "-Wextra",
            "-Werror",
            "-fsyntax-only",
            "-I",
        ])
        .arg(format!("{crate_dir}/include"))
        .arg(&program)
        .output();
    std::fs::remove_file(&program).ok();

    match output {
        Ok(output) => assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        ),
        Err(_) => eprintln!("No C compiler `{compiler}`, header not checked"),
    }
}
