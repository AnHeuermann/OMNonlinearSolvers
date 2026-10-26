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

//! The functions a C host calls.
//!
//! Every function checks its pointers and catches panics, so neither a null
//! pointer nor a panic reaches the host as undefined behaviour.

use crate::model::{CModel, HomModel};
use crate::settings::HomSettings;
use homotopy_core::{Direction, HomotopyError, HomotopyRun, HomotopySolver, Solved};
use std::panic::{AssertUnwindSafe, catch_unwind};

/// A solver for systems of one size, created by [`hom_solver_new`] and
/// opaque to C.
pub struct HomSolver {
    solver: HomotopySolver,
    size: usize,
}

/// Whether and how a solve ended.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HomStatus {
    /// `x` holds the solution
    Solved,
    /// Newton's method and all homotopies failed, `x` is unchanged
    NotSolved,
    /// The model can't be evaluated at the start point or violates a
    /// constraint, `x` is unchanged
    ModelError,
    /// A pointer is `NULL` or the sizes of the model, the solver and `x`
    /// don't match
    InvalidArgument,
    /// The solver failed unexpectedly, a bug in the solver
    Panic,
}

/// How a solve found the solution.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HomSolvedBy {
    /// The start point already solves the system
    StartPoint,
    /// Newton's method from the start point
    Newton,
    /// Newton's method after the Newton homotopy towards increasing `λ`
    HomotopyNewtonPositive,
    /// Newton's method after the Newton homotopy towards decreasing `λ`
    HomotopyNewtonNegative,
    /// Newton's method after the fixpoint homotopy
    HomotopyFixpoint,
    /// The model's own homotopy towards increasing `λ`
    OwnHomotopyPositive,
    /// The model's own homotopy towards decreasing `λ`
    OwnHomotopyNegative,
}

impl From<Solved> for HomSolvedBy {
    fn from(solved: Solved) -> Self {
        match solved {
            Solved::StartPoint => HomSolvedBy::StartPoint,
            Solved::Newton => HomSolvedBy::Newton,
            Solved::Homotopy(HomotopyRun::NewtonPositive) => HomSolvedBy::HomotopyNewtonPositive,
            Solved::Homotopy(HomotopyRun::NewtonNegative) => HomSolvedBy::HomotopyNewtonNegative,
            Solved::Homotopy(HomotopyRun::Fixpoint) => HomSolvedBy::HomotopyFixpoint,
        }
    }
}

impl From<HomotopyError> for HomStatus {
    fn from(error: HomotopyError) -> Self {
        match error {
            HomotopyError::ModelError(_) => HomStatus::ModelError,
            HomotopyError::NotSolved => HomStatus::NotSolved,
        }
    }
}

/// Create a solver for systems with `size` residuals and unknowns.
///
/// `settings` may be `NULL` for the defaults, see [`HomSettings`]. Returns
/// `NULL` if `size` is 0 or a setting is invalid. Free the solver with
/// [`hom_solver_free`].
///
/// # Safety
///
/// `settings` must be `NULL` or point to a valid [`HomSettings`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn hom_solver_new(
    size: usize,
    settings: *const HomSettings,
) -> *mut HomSolver {
    if size == 0 {
        return std::ptr::null_mut();
    }
    // SAFETY: guaranteed by the caller
    let settings = unsafe { settings.as_ref() }.copied().unwrap_or_default();
    let created = catch_unwind(|| {
        let (newton, homotopy) = settings.to_settings(size);
        HomSolver {
            solver: HomotopySolver::new(size, newton, homotopy),
            size,
        }
    });
    match created {
        Ok(solver) => Box::into_raw(Box::new(solver)),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Free a solver of [`hom_solver_new`]. `NULL` is ignored.
///
/// # Safety
///
/// `solver` must be `NULL` or come from [`hom_solver_new`] and not be freed
/// already.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn hom_solver_free(solver: *mut HomSolver) {
    if !solver.is_null() {
        // SAFETY: guaranteed by the caller
        drop(unsafe { Box::from_raw(solver) });
    }
}

/// Solve the system `model` from the start point `x`, like `solveHomotopy`
/// of the OpenModelica C runtime: Newton's method, and if that fails Newton's
/// method after each of the homotopies.
///
/// `x` has `size` entries and holds the solution on success. `solved_by` may
/// be `NULL`, otherwise it receives how the solution was found.
///
/// # Safety
///
/// `solver` must come from [`hom_solver_new`], `model` must point to a valid
/// [`HomModel`] and `x` to `size` values. The callbacks of `model` are called
/// from this thread only, while the function runs.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn hom_solve(
    solver: *mut HomSolver,
    model: *const HomModel,
    x: *mut f64,
    solved_by: *mut HomSolvedBy,
) -> HomStatus {
    // SAFETY: guaranteed by the caller
    let (Some(solver), Some(model)) = (unsafe { solver.as_mut() }, unsafe { model.as_ref() })
    else {
        return HomStatus::InvalidArgument;
    };
    let size = solver.size;
    if x.is_null() || model.num_residuals != size || model.num_unknowns != size {
        return HomStatus::InvalidArgument;
    }
    // SAFETY: guaranteed by the caller
    let Some(mut model) = (unsafe { CModel::new(model) }) else {
        return HomStatus::InvalidArgument;
    };
    // SAFETY: guaranteed by the caller
    let x = unsafe { std::slice::from_raw_parts_mut(x, size) };

    let result = catch_unwind(AssertUnwindSafe(|| solver.solver.solve(&mut model, x)));
    match result {
        Ok(Ok(solved)) => {
            // SAFETY: guaranteed by the caller
            if let Some(solved_by) = unsafe { solved_by.as_mut() } {
                *solved_by = solved.into();
            }
            HomStatus::Solved
        }
        Ok(Err(error)) => error.into(),
        Err(_) => HomStatus::Panic,
    }
}

/// Solve a system `model` that has its own homotopy, with `λ` as last of its
/// `size + 1` unknowns, from the start point `x`.
///
/// The path starts towards decreasing `λ` if `negative_start` is set,
/// otherwise towards increasing `λ`, and the opposite direction is tried if
/// it fails. `x` has `size` entries and holds the solution on success.
/// `solved_by` may be `NULL`, otherwise it receives the direction of the
/// path.
///
/// # Safety
///
/// Like [`hom_solve`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn hom_solve_own_homotopy(
    solver: *mut HomSolver,
    model: *const HomModel,
    x: *mut f64,
    negative_start: bool,
    solved_by: *mut HomSolvedBy,
) -> HomStatus {
    // SAFETY: guaranteed by the caller
    let (Some(solver), Some(model)) = (unsafe { solver.as_mut() }, unsafe { model.as_ref() })
    else {
        return HomStatus::InvalidArgument;
    };
    let size = solver.size;
    if x.is_null() || model.num_residuals != size || model.num_unknowns != size + 1 {
        return HomStatus::InvalidArgument;
    }
    // SAFETY: guaranteed by the caller
    let Some(mut model) = (unsafe { CModel::new(model) }) else {
        return HomStatus::InvalidArgument;
    };
    // SAFETY: guaranteed by the caller
    let x = unsafe { std::slice::from_raw_parts_mut(x, size) };
    let first = if negative_start {
        Direction::Negative
    } else {
        Direction::Positive
    };

    let result = catch_unwind(AssertUnwindSafe(|| {
        solver.solver.solve_own_homotopy(&mut model, x, first)
    }));
    match result {
        Ok(Ok(direction)) => {
            // SAFETY: guaranteed by the caller
            if let Some(solved_by) = unsafe { solved_by.as_mut() } {
                *solved_by = match direction {
                    Direction::Positive => HomSolvedBy::OwnHomotopyPositive,
                    Direction::Negative => HomSolvedBy::OwnHomotopyNegative,
                };
            }
            HomStatus::Solved
        }
        Ok(Err(error)) => error.into(),
        Err(_) => HomStatus::Panic,
    }
}
