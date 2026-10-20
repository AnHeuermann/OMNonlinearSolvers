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
#![allow(dead_code)]

//! Detection of a Newton iteration that makes no progress, from the counters
//! of `newtonAlgorithm` in the C runtime.

use crate::error_history::ErrorHistory;

/// The error grew strongly if it is more than this factor above the previous
/// error
const GROWTH_FACTOR: f64 = 10.0;
/// The error hovers while it stays below this factor times the hover level
const HOVER_WINDOW: f64 = 10.0;
/// A step is heavily damped if its step factor is below this value
const DAMPED_STEP_FACTOR: f64 = 1e-3;
/// An error below this factor times `f_tol²` is inside the "less accurate"
/// band, close enough to a solution that slow progress is no divergence
const LESS_ACCURATE_BAND: f64 = 1e6;
/// The best error improved if the error is below this factor times the best
/// error, by more than 1%
const IMPROVEMENT_FACTOR: f64 = 0.99;

/// Allowed steps where the error grew strongly, over the whole iteration
const MAX_ERROR_GROWTHS: usize = 20;
/// Allowed limit cycles in a row, far from a solution
const MAX_CYCLES_IN_A_ROW: usize = 20;
/// Allowed steps without improving the best error, far from a solution
const MAX_STEPS_WITHOUT_IMPROVEMENT: usize = 400;
/// Allowed heavily damped steps in a row, far from a solution
const MAX_DAMPED_STEPS_IN_A_ROW: usize = 400;

/// Why the Newton iteration moves away from any solution.
#[derive(Debug, PartialEq)]
pub(crate) enum Divergence {
    /// The error grew strongly in more than [`MAX_ERROR_GROWTHS`] steps
    ErrorGrew,
    /// The error cycled in more than [`MAX_CYCLES_IN_A_ROW`] steps in a row
    Cycling,
    /// No improvement for more than [`MAX_STEPS_WITHOUT_IMPROVEMENT`] steps, or
    /// more than [`MAX_DAMPED_STEPS_IN_A_ROW`] heavily damped steps in a row
    NoProgress,
}

/// Tracks the progress of a Newton iteration to detect divergence.
///
/// Updated with the error of every iteration, see [`Tracker::update`].
pub(crate) struct Tracker {
    /// Steps where the error grew strongly, never reset
    num_error_growths: usize,
    /// Limit cycles in a row, far from a solution
    num_cycles_in_a_row: usize,
    /// Steps since the best error last improved by 1%
    num_steps_without_improvement: usize,
    /// Steps since the error started hovering around `hover_error`
    num_steps_hovering: usize,
    /// Heavily damped steps in a row
    num_damped_steps_in_a_row: usize,

    /// The level the error hovers around
    hover_error: f64,
    /// Best error so far
    best_error: f64,

    /// The last errors, to detect limit cycles
    error_history: ErrorHistory,
    /// Error of the previous step
    previous_error: f64,

    /// The latest plain and scaled error are outside the less accurate band
    is_far_from_solution: bool,
    /// Squared residual tolerance, defines the less accurate band
    f_tol_squared: f64,
}

impl Tracker {
    /// Creates a tracker for a new Newton iteration with the residual
    /// tolerance `f_tol`, which defines the less accurate band.
    ///
    /// # Panics
    ///
    /// Panics if `f_tol` is not positive.
    pub(crate) fn new(f_tol: f64) -> Self {
        assert!(f_tol > 0.0, "Invalid tolerance {f_tol}");
        Tracker {
            num_error_growths: 0,
            num_cycles_in_a_row: 0,
            num_steps_without_improvement: 0,
            num_steps_hovering: 0,
            num_damped_steps_in_a_row: 0,

            hover_error: f64::INFINITY,
            best_error: f64::INFINITY,

            error_history: ErrorHistory::new(),
            previous_error: f64::INFINITY,

            is_far_from_solution: true,
            f_tol_squared: f_tol * f_tol,
        }
    }

    /// Updates the counters with the squared error `err` and the squared
    /// scaled error `err_scaled` of the current step, and the step factor
    /// `λ` the line search used.
    pub(crate) fn update(&mut self, err: f64, err_scaled: f64, step_factor: f64) {
        let band = self.f_tol_squared * LESS_ACCURATE_BAND;
        self.is_far_from_solution = err >= band && err_scaled >= band;

        self.count_error_growth(err);
        self.count_cycles(err);
        self.count_progress(err);
        self.count_damped_steps(step_factor);

        self.previous_error = err;
    }

    /// Counts a strong growth of the error since the previous step.
    fn count_error_growth(&mut self, err: f64) {
        if err > GROWTH_FACTOR * self.previous_error {
            self.num_error_growths += 1;
        }
    }

    /// Counts limit cycles in a row, far from a solution. A cycle inside the
    /// less accurate band is left to the convergence check.
    fn count_cycles(&mut self, err: f64) {
        let has_cycle = self.error_history.push(err);
        if has_cycle && self.is_far_from_solution {
            self.num_cycles_in_a_row += 1;
        } else {
            self.num_cycles_in_a_row = 0;
        }
    }

    /// Counts the steps without improvement and the steps hovering.
    ///
    /// An improvement of the best error by 1% resets both, otherwise the error
    /// hovers within 10 times the hover level or jumps away from it.
    fn count_progress(&mut self, err: f64) {
        if err < IMPROVEMENT_FACTOR * self.best_error {
            self.best_error = err;
            self.num_steps_without_improvement = 0;
            self.num_steps_hovering = 0;
        } else if err < HOVER_WINDOW * self.hover_error {
            self.num_steps_without_improvement += 1;
            self.num_steps_hovering += 1;
        } else {
            self.num_steps_without_improvement += 1;
            self.num_steps_hovering = 0;
        }

        // A new hover level after an improvement or a jump, a lower one while
        // hovering
        if self.num_steps_hovering == 0 || err < self.hover_error {
            self.hover_error = err;
        }
    }

    /// Counts heavily damped steps in a row.
    fn count_damped_steps(&mut self, step_factor: f64) {
        if step_factor < DAMPED_STEP_FACTOR {
            self.num_damped_steps_in_a_row += 1;
        } else {
            self.num_damped_steps_in_a_row = 0;
        }
    }

    /// Why the iteration diverges, or `None` while it makes progress.
    pub(crate) fn divergence(&self) -> Option<Divergence> {
        if self.num_error_growths > MAX_ERROR_GROWTHS {
            return Some(Divergence::ErrorGrew);
        }
        if self.num_cycles_in_a_row > MAX_CYCLES_IN_A_ROW {
            return Some(Divergence::Cycling);
        }
        let no_progress = self.num_steps_without_improvement > MAX_STEPS_WITHOUT_IMPROVEMENT
            || self.num_damped_steps_in_a_row > MAX_DAMPED_STEPS_IN_A_ROW;
        if no_progress && self.is_far_from_solution {
            return Some(Divergence::NoProgress);
        }

        None
    }
}

#[cfg(test)]
mod tracker_tests {
    use super::*;

    /// `newtonFTol` of the C runtime. The "less accurate" band is
    /// `e < F_TOL²·10⁶ = 1e-18`.
    const F_TOL: f64 = 1e-12;

    /// An error far outside the less accurate band
    const FAR: f64 = 1.0;
    /// An error inside the less accurate band
    const NEAR: f64 = 1e-20;

    /// Run `update` for each `(err, step_factor)` in `steps`, with the error
    /// of the previous step as `previous_error` and the same value for the plain and
    /// the scaled error. Returns `divergence()` after each step.
    fn run(steps: &[(f64, f64)]) -> Vec<Option<Divergence>> {
        let mut tracker = Tracker::new(F_TOL);
        // C starts with error_f_old = 1e100
        steps
            .iter()
            .map(|&(err, step_factor)| {
                tracker.update(err, err, step_factor);
                tracker.divergence()
            })
            .collect()
    }

    /// `n` full steps (step factor 1) with the errors `err(k)`, `k = 0..n`
    fn full_steps(n: usize, err: impl Fn(usize) -> f64) -> Vec<(f64, f64)> {
        (0..n).map(|k| (err(k), 1.0)).collect()
    }

    #[test]
    /// A steadily decreasing error never diverges
    fn decreasing_error_test() {
        let diverging = run(&full_steps(100, |k| 0.5_f64.powi(k as i32)));

        assert!(diverging.iter().all(Option::is_none));
    }

    #[test]
    /// More than 20 steps where the error grows by more than 10 times diverge
    fn error_growth_test() {
        // The first error has no predecessor, then 21 steps growing 11 times
        let diverging = run(&full_steps(22, |k| 1e-3 * 11_f64.powi(k as i32)));

        assert_eq!(diverging[20], None, "20 error growths are allowed");
        assert_eq!(
            diverging[21],
            Some(Divergence::ErrorGrew),
            "the 21st error growth diverges"
        );
    }

    #[test]
    /// Error growths are counted over all iterations, decreasing steps in
    /// between don't reset the counter
    fn error_growth_not_reset_test() {
        // Alternately grow 11 times and shrink 12 times, 21 times each. The
        // values never repeat within 1%, so this is no cycle.
        let mut steps = vec![(1.0, 1.0)];
        for _ in 0..21 {
            let err = steps.last().unwrap().0;
            steps.push((err * 11.0, 1.0));
            steps.push((err * 11.0 / 12.0, 1.0));
        }

        let diverging = run(&steps);

        assert_eq!(
            diverging[diverging.len() - 3],
            None,
            "20 error growths so far"
        );
        assert_eq!(
            diverging.last().unwrap(),
            &Some(Divergence::ErrorGrew),
            "21 error growths in total"
        );
    }

    #[test]
    /// More than 20 limit cycles in a row, far from a solution, diverge
    fn cycles_far_test() {
        // Period 2: 4 steps fill the history, then every step is a cycle
        let diverging = run(&full_steps(25, |k| {
            if k.is_multiple_of(2) { FAR } else { 2.0 * FAR }
        }));

        assert_eq!(diverging[23], None, "20 cycles in a row are allowed");
        assert_eq!(
            diverging[24],
            Some(Divergence::Cycling),
            "the 21st cycle in a row diverges"
        );
    }

    #[test]
    /// Limit cycles inside the less accurate band are left to the stopping
    /// decision, they don't count as diverging
    fn cycles_near_test() {
        let diverging = run(&full_steps(40, |k| {
            if k.is_multiple_of(2) {
                NEAR
            } else {
                2.0 * NEAR
            }
        }));

        assert!(diverging.iter().all(Option::is_none));
    }

    #[test]
    /// A cycle only counts if the scaled error is outside the less accurate
    /// band as well
    fn cycles_scaled_error_near_test() {
        let mut tracker = Tracker::new(F_TOL);
        for k in 0..40_usize {
            let err = if k.is_multiple_of(2) { FAR } else { 2.0 * FAR };
            tracker.update(err, NEAR, 1.0);

            assert_eq!(tracker.divergence(), None, "step {k}");
        }
    }

    #[test]
    /// A step that is no cycle resets the count of cycles in a row
    fn cycles_reset_test() {
        let alternate = |k: usize| if k.is_multiple_of(2) { FAR } else { 2.0 * FAR };
        // 4 steps fill the history, then 20 cycles
        let mut steps = full_steps(24, alternate);
        // 1.5 returns to no older error: no cycle, the count is reset
        steps.push((1.5 * FAR, 1.0));
        // 20 more cycles, starting with 1.0 after 1.5
        steps.extend(full_steps(20, alternate));

        let diverging = run(&steps);

        assert!(diverging.iter().all(Option::is_none));
    }

    #[test]
    /// More than 400 steps without improving the best error by 1%, far from a
    /// solution, diverge
    fn no_improvement_far_test() {
        // The first step sets the best error, then 401 steps don't improve it
        let diverging = run(&full_steps(402, |_| FAR));

        assert_eq!(
            diverging[400], None,
            "400 steps without improvement are allowed"
        );
        assert_eq!(
            diverging[401],
            Some(Divergence::NoProgress),
            "the 401st diverges"
        );
    }

    #[test]
    /// Steps without improvement inside the less accurate band don't diverge
    fn no_improvement_near_test() {
        let diverging = run(&full_steps(500, |_| NEAR));

        assert!(diverging.iter().all(Option::is_none));
    }

    #[test]
    /// Improving the best error by more than 1% resets the count
    fn no_improvement_reset_test() {
        // 400 steps without improvement, an improvement by 2%, then 400 more
        let mut steps = full_steps(401, |_| FAR);
        steps.extend(full_steps(401, |_| 0.98 * FAR));

        let diverging = run(&steps);

        assert!(diverging.iter().all(Option::is_none));
    }

    #[test]
    /// More than 400 heavily damped steps (step factor below `1e-3`) in a row,
    /// far from a solution, diverge
    fn damped_steps_far_test() {
        // The error improves by 2% per step, so every step is an improvement
        let steps: Vec<_> = (0..401).map(|k| (0.98_f64.powi(k), 1e-4)).collect();

        let diverging = run(&steps);

        assert_eq!(diverging[399], None, "400 heavily damped steps are allowed");
        assert_eq!(
            diverging[400],
            Some(Divergence::NoProgress),
            "the 401st heavily damped step diverges"
        );
    }

    #[test]
    /// A step that isn't heavily damped resets the count of damped steps
    fn damped_steps_reset_test() {
        let mut steps: Vec<_> = (0..400).map(|k| (0.98_f64.powi(k), 1e-4)).collect();
        steps.push((0.98_f64.powi(400), 1.0));
        steps.extend((401..801).map(|k| (0.98_f64.powi(k), 1e-4)));

        let diverging = run(&steps);

        assert!(diverging.iter().all(Option::is_none));
    }

    #[test]
    /// Heavily damped steps inside the less accurate band don't diverge
    fn damped_steps_near_test() {
        let steps: Vec<_> = (0..500).map(|k| (NEAR * 0.98_f64.powi(k), 1e-4)).collect();

        let diverging = run(&steps);

        assert!(diverging.iter().all(Option::is_none));
    }

    /// Run `update` for each error with full steps, return the tracker
    fn track(errors: &[f64]) -> Tracker {
        let mut tracker = Tracker::new(F_TOL);
        for &err in errors {
            tracker.update(err, err, 1.0);
        }
        tracker
    }

    #[test]
    /// An improvement of the best error by 1% resets the hover count, and the
    /// error becomes the new hover level
    fn hover_reset_by_improvement_test() {
        let tracker = track(&[1.0, 2.0, 3.0, 0.5]);

        assert_eq!(tracker.num_steps_hovering, 0);
        assert_eq!(tracker.hover_error, 0.5);
    }

    #[test]
    /// Errors within 10 times of the hover level count as hovering, the level
    /// stays while the error is above it
    fn hover_count_test() {
        let tracker = track(&[1.0, 2.0, 3.0, 9.0]);

        assert_eq!(tracker.num_steps_hovering, 3);
        assert_eq!(tracker.hover_error, 1.0);
    }

    #[test]
    /// While hovering, an error below the hover level lowers it. 0.995 is no
    /// improvement of the best error 1.0 by 1%.
    fn hover_level_lowered_test() {
        let tracker = track(&[1.0, 2.0, 0.995]);

        assert_eq!(tracker.num_steps_hovering, 2);
        assert_eq!(tracker.hover_error, 0.995);
    }

    #[test]
    /// A jump above 10 times the hover level resets the hover count and sets
    /// a new hover level
    fn hover_jump_test() {
        let tracker = track(&[1.0, 2.0, 20.0]);
        assert_eq!(tracker.num_steps_hovering, 0);
        assert_eq!(tracker.hover_error, 20.0);

        let tracker = track(&[1.0, 2.0, 20.0, 30.0]);
        assert_eq!(tracker.num_steps_hovering, 1);
        assert_eq!(tracker.hover_error, 20.0);
    }

    #[test]
    #[should_panic = "Invalid tolerance"]
    /// The residual tolerance must be positive
    fn invalid_tolerance_test() {
        let _ = Tracker::new(0.0);
    }
}
