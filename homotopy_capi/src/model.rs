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

//! A system described by C callbacks, as a model of the solvers.

use faer::MatMut;
use nls_model::{EvalError, ModelFunctions, SystemInfo, forward_differences};
use std::ffi::{c_int, c_void};

/// The callback evaluated successfully
pub const HOM_EVAL_OK: c_int = 0;
/// The model asserted, the solver tries another point
pub const HOM_EVAL_ASSERT: c_int = 1;
/// A casual tearing set violated a constraint, the solver gives up
pub const HOM_EVAL_CONSTRAINT_VIOLATED: c_int = 2;

/// Evaluate the residuals `f = F(x)`.
///
/// `x` has `num_unknowns` entries, `f` has `num_residuals`. Returns
/// [`HOM_EVAL_OK`], [`HOM_EVAL_ASSERT`] or [`HOM_EVAL_CONSTRAINT_VIOLATED`],
/// other values count as an assert. The callback must return: it must not
/// leave with `longjmp` or an exception, since that would cross the solver's
/// Rust stack frames.
pub type HomResidual =
    Option<unsafe extern "C" fn(user: *mut c_void, x: *const f64, f: *mut f64) -> c_int>;

/// Evaluate the Jacobian `jacobian = J(x)` with `f = F(x)`.
///
/// `jacobian` is a column-major `num_residuals×num_unknowns` matrix. Returns
/// a value like [`HomResidual`] and must return like it.
pub type HomJacobian = Option<
    unsafe extern "C" fn(
        user: *mut c_void,
        x: *const f64,
        f: *const f64,
        jacobian: *mut f64,
    ) -> c_int,
>;

/// A nonlinear system `F(x) = 0` of a C host.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct HomModel {
    /// Passed to the callbacks as it is
    pub user: *mut c_void,
    /// Number of residuals
    pub num_residuals: usize,
    /// Number of unknowns, with `λ` as last one for a model with its own
    /// homotopy
    pub num_unknowns: usize,
    /// Nominal values of the unknowns, or `NULL` for 1
    pub nominal: *const f64,
    /// Lower bounds of the unknowns, or `NULL` for none. Not used by the
    /// solver, like in the C runtime.
    pub min: *const f64,
    /// Upper bounds of the unknowns, or `NULL` for none. The forward
    /// differences step backwards where they would pass one.
    pub max: *const f64,
    /// The system is the casual set of a dynamic tearing
    pub casual_tearing_set: bool,
    /// Evaluates `F(x)`, required
    pub residual: HomResidual,
    /// Evaluates `J(x)`, or `NULL` for forward differences
    pub jacobian: HomJacobian,
}

/// A [`HomModel`] as a model of the solvers.
pub(crate) struct CModel<'m> {
    model: &'m HomModel,
    residual: unsafe extern "C" fn(user: *mut c_void, x: *const f64, f: *mut f64) -> c_int,
    info: SystemInfo,
    /// Column-major buffer of the Jacobian callback
    jacobian: Box<[f64]>,
}

impl<'m> CModel<'m> {
    /// Wraps `model`, or returns `None` if it has no residual callback or no
    /// unknowns.
    ///
    /// # Safety
    ///
    /// `nominal`, `min` and `max` of `model` must be `NULL` or point to
    /// `num_unknowns` values that stay valid while the solver runs.
    pub(crate) unsafe fn new(model: &'m HomModel) -> Option<Self> {
        let residual = model.residual?;
        let (num_residuals, num_unknowns) = (model.num_residuals, model.num_unknowns);
        if num_unknowns == 0 {
            return None;
        }

        // SAFETY: guaranteed by the caller
        let values = |pointer: *const f64| -> Option<Box<[f64]>> {
            (!pointer.is_null())
                .then(|| unsafe { std::slice::from_raw_parts(pointer, num_unknowns) }.into())
        };
        let mut info = SystemInfo::new(num_residuals, num_unknowns);
        if let Some(nominal) = values(model.nominal) {
            info = info.with_nominal(nominal);
        }
        if let Some(min) = values(model.min) {
            info = info.with_min(min);
        }
        if let Some(max) = values(model.max) {
            info = info.with_max(max);
        }
        if model.casual_tearing_set {
            info = info.with_casual_tearing_set();
        }

        Some(CModel {
            model,
            residual,
            info,
            jacobian: vec![0.0; num_residuals * num_unknowns].into(),
        })
    }
}

/// The result of a callback
fn eval_result(status: c_int) -> Result<(), EvalError> {
    match status {
        HOM_EVAL_OK => Ok(()),
        HOM_EVAL_CONSTRAINT_VIOLATED => Err(EvalError::ConstraintViolated),
        _ => Err(EvalError::Assert),
    }
}

impl ModelFunctions for CModel<'_> {
    fn info(&self) -> &SystemInfo {
        &self.info
    }

    fn eval(&mut self, x: &[f64], res: &mut [f64]) -> Result<(), EvalError> {
        assert_eq!(x.len(), self.info.num_unknowns);
        assert_eq!(res.len(), self.info.num_residuals);
        // SAFETY: the slices have the sizes the callback expects
        let status = unsafe { (self.residual)(self.model.user, x.as_ptr(), res.as_mut_ptr()) };
        eval_result(status)
    }

    fn derivative(
        &mut self,
        x: &[f64],
        fx: &[f64],
        jacobian: MatMut<'_, f64>,
    ) -> Result<(), EvalError> {
        let Some(callback) = self.model.jacobian else {
            return forward_differences(self, x, fx, jacobian);
        };
        assert_eq!(x.len(), self.info.num_unknowns);
        assert_eq!(fx.len(), self.info.num_residuals);
        assert_eq!(
            (jacobian.nrows(), jacobian.ncols()),
            (fx.len(), x.len()),
            "Jacobian has wrong shape"
        );

        // SAFETY: the buffer is column-major with the shape the callback
        // expects, `x` and `fx` have the sizes of the model
        let status = unsafe {
            callback(
                self.model.user,
                x.as_ptr(),
                fx.as_ptr(),
                self.jacobian.as_mut_ptr(),
            )
        };
        eval_result(status)?;

        let rows = fx.len();
        for (j, column) in jacobian.col_iter_mut().enumerate() {
            let values = &self.jacobian[j * rows..(j + 1) * rows];
            for (entry, value) in column.iter_mut().zip(values) {
                *entry = *value;
            }
        }
        Ok(())
    }
}
