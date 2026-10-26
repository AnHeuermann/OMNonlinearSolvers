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

/* Smoke test of the C interface: one function per way a C host uses it. */

#include <math.h>
#include <stdio.h>

#include "homotopy_capi.h"

static int failures = 0;

#define CHECK(condition)                                                    \
  do {                                                                      \
    if (!(condition)) {                                                     \
      fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__,      \
              #condition);                                                  \
      failures++;                                                           \
    }                                                                       \
  } while (0)

/* F(x) = x³ − 2x + 2, Newton's method alone gets stuck from x = 0 */
static int cubic(void *user, const double *x, double *f) {
  (void)user;
  f[0] = x[0] * x[0] * x[0] - 2.0 * x[0] + 2.0;
  return HOM_EVAL_OK;
}

/* J(x) = 3x² − 2, counting its calls in user */
static int cubic_jacobian(void *user, const double *x, const double *f,
                          double *jacobian) {
  (void)f;
  (*(int *)user)++;
  jacobian[0] = 3.0 * x[0] * x[0] - 2.0;
  return HOM_EVAL_OK;
}

/* A model that can't be evaluated anywhere */
static int asserts(void *user, const double *x, double *f) {
  (void)user;
  (void)x;
  (void)f;
  return HOM_EVAL_ASSERT;
}

/* The model's own homotopy H(x, λ) = x − 2λ, with λ as last unknown */
static int own_homotopy(void *user, const double *y, double *h) {
  (void)user;
  h[0] = y[0] - 2.0 * y[1];
  return HOM_EVAL_OK;
}

/* The root of the cubic */
static const double ROOT = -1.769292354238631;

static void solve_with_homotopy(void) {
  HomSolver *solver = hom_solver_new(1, NULL);
  HomModel model = {0};
  model.num_residuals = model.num_unknowns = 1;
  model.residual = cubic;
  double x[1] = {0.0};
  HomSolvedBy solved_by;

  CHECK(hom_solve(solver, &model, x, &solved_by) == HOM_STATUS_SOLVED);
  CHECK(solved_by == HOM_SOLVED_BY_HOMOTOPY_NEWTON_NEGATIVE);
  CHECK(fabs(x[0] - ROOT) < 1e-10);
  hom_solver_free(solver);
}

static void solve_with_jacobian(void) {
  HomSolver *solver = hom_solver_new(1, NULL);
  int calls = 0;
  HomModel model = {0};
  model.user = &calls;
  model.num_residuals = model.num_unknowns = 1;
  model.residual = cubic;
  model.jacobian = cubic_jacobian;
  double x[1] = {0.0};

  CHECK(hom_solve(solver, &model, x, NULL) == HOM_STATUS_SOLVED);
  CHECK(fabs(x[0] - ROOT) < 1e-10);
  CHECK(calls > 0);
  hom_solver_free(solver);
}

static void model_error(void) {
  HomSolver *solver = hom_solver_new(1, NULL);
  HomModel model = {0};
  model.num_residuals = model.num_unknowns = 1;
  model.residual = asserts;
  double x[1] = {1.0};

  CHECK(hom_solve(solver, &model, x, NULL) == HOM_STATUS_MODEL_ERROR);
  CHECK(x[0] == 1.0);
  hom_solver_free(solver);
}

static void solve_own_homotopy(void) {
  HomSolver *solver = hom_solver_new(1, NULL);
  HomModel model = {0};
  model.num_residuals = 1;
  model.num_unknowns = 2;
  model.residual = own_homotopy;
  double x[1] = {0.0};
  HomSolvedBy solved_by;

  CHECK(hom_solve_own_homotopy(solver, &model, x, false, &solved_by) ==
        HOM_STATUS_SOLVED);
  CHECK(solved_by == HOM_SOLVED_BY_OWN_HOMOTOPY_POSITIVE);
  CHECK(fabs(x[0] - 2.0) < 1e-5);
  hom_solver_free(solver);
}

static void settings(void) {
  HomSettings settings = {0};
  settings.tau_start = 0.1;
  settings.max_tries = 5;
  HomSolver *solver = hom_solver_new(1, &settings);
  CHECK(solver != NULL);
  hom_solver_free(solver);

  settings.tau_decreasing_factor = 0.5; /* must be greater than 1 */
  CHECK(hom_solver_new(1, &settings) == NULL);
  CHECK(hom_solver_new(0, NULL) == NULL);
}

static void invalid_arguments(void) {
  HomSolver *solver = hom_solver_new(1, NULL);
  HomModel model = {0};
  model.num_residuals = model.num_unknowns = 1;
  double x[1] = {0.0};

  CHECK(hom_solve(solver, &model, x, NULL) == HOM_STATUS_INVALID_ARGUMENT);
  model.residual = cubic;
  CHECK(hom_solve(NULL, &model, x, NULL) == HOM_STATUS_INVALID_ARGUMENT);
  CHECK(hom_solve(solver, NULL, x, NULL) == HOM_STATUS_INVALID_ARGUMENT);
  CHECK(hom_solve(solver, &model, NULL, NULL) == HOM_STATUS_INVALID_ARGUMENT);
  hom_solver_free(solver);
  hom_solver_free(NULL);
}

int main(void) {
  solve_with_homotopy();
  solve_with_jacobian();
  model_error();
  solve_own_homotopy();
  settings();
  invalid_arguments();

  if (failures > 0) {
    fprintf(stderr, "%d checks failed\n", failures);
    return 1;
  }
  printf("All checks passed\n");
  return 0;
}
