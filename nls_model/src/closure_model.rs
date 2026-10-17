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

//! A model built from a closure, for small systems and tests.

use crate::{EvalError, ModelFunctions, SystemInfo};

/// A model built from a closure that evaluates `F(x)`.
///
/// Saves writing a type and a trait implementation for small systems, for
/// example in tests. The closure has the shape of [`ModelFunctions::eval`]:
/// it writes `F(x)` into `res` and may fail. Nominal values, bounds and the
/// casual tearing flag are optional and set with the `with_*` methods, which
/// set them in the model's [`SystemInfo`].
///
/// # Example
///
/// ```
/// use nls_model::{ClosureModel, ModelFunctions};
///
/// // F(x, y) = (x² - y, x + y - 2) with x ≥ 0
/// let mut model = ClosureModel::new(2, 2, |x: &[f64], res: &mut [f64]| {
///     res[0] = x[0] * x[0] - x[1];
///     res[1] = x[0] + x[1] - 2.0;
///     Ok(())
/// })
/// .with_min(vec![0.0, f64::NEG_INFINITY]);
///
/// let mut res = [0.0; 2];
/// model.eval(&[1.0, 1.0], &mut res).unwrap();
/// assert_eq!(res, [0.0, 0.0]);
/// ```
pub struct ClosureModel<F> {
    /// Evaluates `F(x)` into `res`
    f: F,
    /// Facts about the system
    info: SystemInfo,
}

impl<F> ClosureModel<F>
where
    F: FnMut(&[f64], &mut [f64]) -> Result<(), EvalError>,
{
    /// Creates a model with `num_residuals` residuals and `num_unknowns`
    /// unknowns, evaluated by `f`, with the defaults of [`SystemInfo::new`].
    pub fn new(num_residuals: usize, num_unknowns: usize, f: F) -> Self {
        ClosureModel {
            f,
            info: SystemInfo::new(num_residuals, num_unknowns),
        }
    }

    /// Marks the model as the casual set of a dynamic tearing, see
    /// [`SystemInfo::with_casual_tearing_set`].
    pub fn with_casual_tearing_set(mut self) -> Self {
        self.info = self.info.with_casual_tearing_set();
        self
    }

    /// Sets the nominal values of the unknowns, see
    /// [`SystemInfo::with_nominal`].
    ///
    /// # Panics
    ///
    /// Panics if `nominal` doesn't have `num_unknowns` entries.
    pub fn with_nominal(mut self, nominal: impl Into<Box<[f64]>>) -> Self {
        self.info = self.info.with_nominal(nominal);
        self
    }

    /// Sets the lower bounds of the unknowns, see [`SystemInfo::with_min`].
    ///
    /// # Panics
    ///
    /// Panics if `min` doesn't have `num_unknowns` entries.
    pub fn with_min(mut self, min: impl Into<Box<[f64]>>) -> Self {
        self.info = self.info.with_min(min);
        self
    }

    /// Sets the upper bounds of the unknowns, see [`SystemInfo::with_max`].
    ///
    /// # Panics
    ///
    /// Panics if `max` doesn't have `num_unknowns` entries.
    pub fn with_max(mut self, max: impl Into<Box<[f64]>>) -> Self {
        self.info = self.info.with_max(max);
        self
    }
}

impl<F> ModelFunctions for ClosureModel<F>
where
    F: FnMut(&[f64], &mut [f64]) -> Result<(), EvalError>,
{
    fn info(&self) -> &SystemInfo {
        &self.info
    }

    fn eval(&mut self, x: &[f64], res: &mut [f64]) -> Result<(), EvalError> {
        (self.f)(x, res)
    }
}
