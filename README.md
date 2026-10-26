# Non-linear Homotopy Solver

[![CI][ci-badge]][ci]

Non-linear solver Homotopy, based on C file
[nonlinearSolverHomotopy.c][c-source] originally created by Bernhard Bachmann
for the OpenModelica project.

## Solve NLS

A nonlinear system `F(x) = 0` implements `ModelFunctions` of the `nls_model`
crate: its sizes and nominal values, and how to evaluate `F`. The Jacobian
defaults to forward differences. `HomotopySolver` of the `homotopy_core`
crate solves it like the C runtime: Newton's method from the start point, and
if that fails Newton's method from the end of a homotopy path.

```rust
use homotopy_core::HomotopySolver;
use homotopy_core::settings::Settings;
use newton::settings::Settings as NewtonSettings;

let mut solver = HomotopySolver::new(size, NewtonSettings::new(size), Settings::new(size));
let how = solver.solve(&mut model, &mut x)?;
```

[homotopy_core/examples/solve.rs](homotopy_core/examples/solve.rs) is a
complete example, run it with

```sh
cargo run --example solve -p homotopy_core
```

A system Newton's method solves on its own can also use `newton::solve`
directly.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for how to build, test and develop this package.
The design decisions behind the code are collected in [DESIGN.md](DESIGN.md).

## C Interface

The `homotopy_capi` crate builds a static library with a generic C interface.
Its header `homotopy_capi/include/homotopy_capi.h` is generated when the crate
is built. A host describes its system as a `HomModel` of callbacks:

```c
#include <stdio.h>

#include "homotopy_capi.h"

/* F(x) = x³ − 2x + 2 */
static int residual(void *user, const double *x, double *f) {
  (void)user; /* the host's data, from model.user */
  f[0] = x[0] * x[0] * x[0] - 2.0 * x[0] + 2.0;
  return HOM_EVAL_OK; /* or HOM_EVAL_ASSERT, HOM_EVAL_CONSTRAINT_VIOLATED */
}

int main(void) {
  HomSettings settings = {0}; /* 0 keeps a default */
  HomSolver *solver = hom_solver_new(1, &settings);

  HomModel model = {0};
  model.num_residuals = model.num_unknowns = 1;
  model.residual = residual;

  double x[1] = {0.0}; /* the start point */
  HomSolvedBy solved_by;
  HomStatus status = hom_solve(solver, &model, x, &solved_by);
  if (status == HOM_STATUS_SOLVED) {
    printf("x = %g\n", x[0]);
  }

  hom_solver_free(solver);
  return status == HOM_STATUS_SOLVED ? 0 : 1;
}
```

The callbacks must return and must not leave with `longjmp`. Link the static
library together with the system libraries it needs, on Linux
`-lpthread -ldl -lm`. CI compiles and runs this example, see
[homotopy_capi/tests/c](homotopy_capi/tests/c), which also has a test of the
rest of the interface. The glue to the OpenModelica C runtime lives in
OpenModelica and uses this interface.

## License

This package is part of the OpenModelica Run-Time System and is distributed, at your choice,
under the terms of one of the following licenses:

* [BSD 3-Clause License (BSD New)][bsd-3-clause]
* [GNU Affero General Public License (AGPL) version 3][agpl-3]
* [OSMC Public License (OSMC-PL) version 1.8][osmc-pl]

This is the same choice as for the run-time system of [OpenModelica][openmodelica],
see the [OpenModelica Run-time System Public License][osmc-runtime].

[agpl-3]: AGPL-License.txt
[bsd-3-clause]: BSD-License.txt
[c-source]: https://github.com/OpenModelica/OpenModelica/blob/18d4ff455d1beceb4709d770ca4faf53d65732da/OMCompiler/SimulationRuntime/c/simulation/solver/nonlinearSolverHomotopy.c
[ci-badge]: https://github.com/AnHeuermann/homotopy_nls_solver/actions/workflows/ci.yml/badge.svg
[ci]: https://github.com/AnHeuermann/homotopy_nls_solver/actions/workflows/ci.yml
[openmodelica]: https://openmodelica.org/
[osmc-pl]: OSMC-License.txt
[osmc-runtime]: OSMC-Runtime-License.txt
