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

use super::*;
use equator::assert;
use faer::utils::approx::{ApproxEq, CwiseMat};
use faer::{Mat, mat};

#[test]
/// Evaluate a Jacobian of a non-square system with 2 residuals and 3
/// unknowns, shaped like a homotopy `H(x, λ)` with `λ` as additional
/// unknown:
/// ```text
/// H(x, y, λ) = ( x² + y - λ )
///              ( x·y·λ - 2  )
/// ```
/// at `(x, y, λ) = (1, 2, 1)` and `(2, -1, 0.5)`
fn test_numeric_jacobian() {
    /// `H: ℝ³ → ℝ²`
    struct H {
        info: SystemInfo,
    }
    impl ModelFunctions for H {
        fn info(&self) -> &SystemInfo {
            &self.info
        }

        fn eval(&mut self, x_vec: &[f64], res: &mut [f64]) -> Result<(), EvalError> {
            let [x, y, lambda] = [x_vec[0], x_vec[1], x_vec[2]];

            res[0] = x * x + y - lambda;
            res[1] = x * y * lambda - 2.0;
            Ok(())
        }
    }

    let mut model = H {
        info: SystemInfo::new(2, 3),
    };
    let mut f_x = [0.0; 2];
    let mut jac = Mat::<f64>::zeros(2, 3);

    // Forward differences are only accurate to about the step size
    let approx_eq = CwiseMat(ApproxEq {
        abs_tol: 1e-6,
        rel_tol: 1e-6,
    });

    // J(x, y, λ) = [ 2x   1    -1 ]
    //              [ yλ   xλ   xy ]
    let x = [1.0, 2.0, 1.0];
    model.eval(&x, &mut f_x).unwrap();
    model.derivative(&x, &f_x, jac.as_mut()).unwrap();
    assert!(jac ~ mat![[2.0, 1.0, -1.0], [2.0, 1.0, 2.0]]);

    let x = [2.0, -1.0, 0.5];
    model.eval(&x, &mut f_x).unwrap();
    model.derivative(&x, &f_x, jac.as_mut()).unwrap();
    assert!(jac ~ mat![[4.0, 1.0, -1.0], [-0.5, 1.0, -2.0]]);
}

#[test]
/// A `ClosureModel` evaluates its closure and reports its sizes and
/// settings.
fn closure_model() {
    // F(x) = x² - 4 with x ≥ 0, as casual tearing set
    let mut model = ClosureModel::new(1, 1, |x: &[f64], res: &mut [f64]| {
        if x[0] < 0.0 {
            return Err(EvalError::ConstraintViolated);
        }
        res[0] = x[0] * x[0] - 4.0;
        Ok(())
    })
    .with_casual_tearing_set()
    .with_nominal([2.0])
    .with_min([0.0]);

    assert!(model.info().num_residuals == 1);
    assert!(model.info().num_unknowns == 1);
    assert!(model.info().casual_tearing_set);

    let mut res = [0.0];
    model.eval(&[3.0], &mut res).unwrap();
    assert!(res == [5.0]);
    let result = model.eval(&[-1.0], &mut res);
    assert!(matches!(result, Err(EvalError::ConstraintViolated)));

    assert_eq!(*model.info().nominal, [2.0]);
    assert_eq!(*model.info().min, [0.0]);
    // Not set: the default, no upper bound
    assert_eq!(*model.info().max, [f64::INFINITY]);
}

#[test]
/// The closure can capture and change state, like a struct implementing
/// the trait with fields
fn closure_model_with_state() {
    let mut counter = 0;
    let mut model = ClosureModel::new(1, 1, |x: &[f64], res: &mut [f64]| {
        counter += 1;
        res[0] = x[0] - 1.0;
        Ok(())
    });

    let mut res = [0.0];
    model.eval(&[1.0], &mut res).unwrap();
    model.eval(&[2.0], &mut res).unwrap();
    drop(model);

    assert!(counter == 2);
}

#[test]
#[should_panic = "Expected one value per unknown"]
/// Settings must have one value per unknown
fn closure_model_wrong_length() {
    let _ = ClosureModel::new(1, 2, |_: &[f64], _: &mut [f64]| Ok(())).with_nominal([1.0]);
}

/// `ClosureModel` with a closure that captures nothing, stored as function
/// pointer
type SimpleClosureModel = ClosureModel<fn(&[f64], &mut [f64]) -> Result<(), EvalError>>;

/// `ClosureModel` for `F(x, y, z) = (x + y + z)` with 1 residual and 3 unknowns
fn three_unknowns() -> SimpleClosureModel {
    ClosureModel::new(1, 3, |x: &[f64], res: &mut [f64]| {
        res[0] = x[0] + x[1] + x[2];
        Ok(())
    })
}

#[test]
/// Without settings the model has the defaults of `SystemInfo::new`: nominal
/// values of 1 and no bounds
fn closure_model_defaults() {
    let model = three_unknowns();

    assert_eq!(*model.info().nominal, [1.0; 3]);
    assert_eq!(*model.info().min, [f64::NEG_INFINITY; 3]);
    assert_eq!(*model.info().max, [f64::INFINITY; 3]);
    assert!(!model.info().casual_tearing_set);
}

#[test]
/// Nominal values, lower and upper bounds set for every unknown
fn closure_model_nominal_min_max() {
    let model = three_unknowns()
        .with_nominal([1.0, 10.0, 1e-3])
        .with_min([0.0, -5.0, f64::NEG_INFINITY])
        .with_max([1.0, 5.0, 2.0]);

    assert_eq!(*model.info().nominal, [1.0, 10.0, 1e-3]);
    assert_eq!(*model.info().min, [0.0, -5.0, f64::NEG_INFINITY]);
    assert_eq!(*model.info().max, [1.0, 5.0, 2.0]);
}

#[test]
/// The values can be given as array, `Vec` or boxed slice
fn closure_model_value_types() {
    let model = three_unknowns()
        .with_nominal([2.0; 3])
        .with_min(vec![-1.0; 3])
        .with_max(vec![1.0; 3].into_boxed_slice());

    assert_eq!(*model.info().nominal, [2.0; 3]);
    assert_eq!(*model.info().min, [-1.0; 3]);
    assert_eq!(*model.info().max, [1.0; 3]);
}

#[test]
/// Setting a value twice keeps the last one, the order of the `with_*` calls
/// doesn't matter otherwise
fn closure_model_set_twice() {
    let model = three_unknowns()
        .with_max([1.0; 3])
        .with_nominal([5.0; 3])
        .with_max([3.0; 3]);

    assert_eq!(*model.info().nominal, [5.0; 3]);
    assert_eq!(*model.info().max, [3.0; 3]);
}

#[test]
#[should_panic = "Expected one value per unknown"]
/// Lower bounds must have one value per unknown
fn closure_model_min_wrong_length() {
    let _ = three_unknowns().with_min([0.0; 2]);
}

#[test]
#[should_panic = "Expected one value per unknown"]
/// Upper bounds must have one value per unknown
fn closure_model_max_wrong_length() {
    let _ = three_unknowns().with_max([0.0; 4]);
}

#[test]
/// `SystemInfo::new` has the defaults: nominal values of 1, no bounds and no
/// casual tearing set
fn system_info_defaults() {
    let info = SystemInfo::new(2, 3);

    assert!(info.num_residuals == 2);
    assert!(info.num_unknowns == 3);
    assert_eq!(*info.nominal, [1.0; 3]);
    assert_eq!(*info.min, [f64::NEG_INFINITY; 3]);
    assert_eq!(*info.max, [f64::INFINITY; 3]);
    assert!(!info.casual_tearing_set);
}

#[test]
/// The `with_*` methods of `SystemInfo` set the values and the tearing kind
fn system_info_with() {
    let info = SystemInfo::new(1, 2)
        .with_nominal([1.0, 1e3])
        .with_min([0.0, -1.0])
        .with_max([1.0, 1.0])
        .with_casual_tearing_set();

    assert_eq!(*info.nominal, [1.0, 1e3]);
    assert_eq!(*info.min, [0.0, -1.0]);
    assert_eq!(*info.max, [1.0, 1.0]);
    assert!(info.casual_tearing_set);
}
