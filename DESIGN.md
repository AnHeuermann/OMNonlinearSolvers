# Design Decisions

The core decisions behind this project and the reasons for them. Add new
decisions at the end of their section, and change an entry when a decision
changes, in the same commit as the code.

## Goals

### Learning project and drop-in replacement

The project ports [nonlinearSolverHomotopy.c][c-source] to Rust. It is a
project to learn Rust, and it shall later replace both the C solver and the
existing Rust port in `OMCompiler/SimulationRuntime/rust/openmodelica_nls` of
OpenModelica.

### Same behaviour as the C solver

The port reproduces the C algorithm, including its constants, thresholds, the
order of its checks and the tie-breaking of the pivot search. Changing which
systems converge would change simulation results in OpenModelica.

* Expected values in tests come from the C algorithm, computed with a port of
  the C function where they are not obvious.
* Deviations from C are deliberate and documented where they are made, for
  example rejecting a matrix with NaN before the elimination.

## Structure

### A workspace of small crates

Each crate has one job and can be used on its own:

| Crate                | Contents                                             |
| -------------------- | ---------------------------------------------------- |
| `nls_model`          | The interface to a nonlinear system, scaling         |
| `permuted_matrix`    | A matrix view with permuted rows and columns         |
| `total_pivot_solver` | Linear solver with complete pivoting                 |
| `newton`             | Damped Newton-Raphson solver                         |
| `homotopy_core`      | The homotopy algorithm                               |
| `homotopy_capi`      | The C interface for the OpenModelica C runtime       |

Dependencies point one way, from the C interface down to the matrix types.
The linear solver crates know nothing about models. Newton knows nothing
about the homotopy.

### Pure Rust core, all `unsafe` in the C interface

All crates except `homotopy_capi` are safe Rust without C types, testable
without OpenModelica. `homotopy_capi` exports the functions of
`nonlinearSolverHomotopy.h`, so `nonlinearSystem.c` doesn't change. It holds
the bindings, the pointers into the runtime and a C wrapper that catches the
`longjmp` of a failed model assert before it crosses Rust stack frames.

Runtime state that is not numerics, like the relations of mixed systems and
`discreteCall`, stays in `homotopy_capi`.

### The line search is a module of `newton`, not a crate yet

Only Newton damps its steps, so `newton/src/line_search/` stays a module until
a second solver needs it. It only depends on `nls_model` and takes slices and
the slope `φ'(0)` instead of Newton's workspace and errors, so it can move into
its own crate without changes.

## Model Interface

### A trait implemented by each system

A nonlinear system implements `ModelFunctions`. Its fields hold whatever it
needs to evaluate itself: parameters of a test system, or pointers into the C
runtime. There is no generic model struct with user data and function
pointers.

* Evaluation takes `&mut self`, because evaluating a C model writes into the
  runtime and models may keep buffers or counters.
* Results are written into output slices, so the solvers reuse their buffers.
* Evaluation returns a `Result`. `EvalError` distinguishes a model assert,
  after which Newton shrinks the step, from a violated constraint of a casual
  tearing set, after which it gives up. Constraints are checked in `eval`
  itself, there is no second evaluation function.
* The Jacobian defaults to forward differences with the step of the C
  runtime, stepping backwards at an upper bound. Systems with an analytical
  Jacobian override it.

### Facts about a system in `SystemInfo`

Sizes, nominal values, bounds and whether the system is a casual tearing set
are a `SystemInfo` struct with public fields, returned by
`ModelFunctions::info`. A trait can't require fields, and one struct checks
the lengths of all values in one place. Its constructor sets the defaults:
nominal values of 1, no bounds, no casual tearing set.

### `ClosureModel` for small systems

`ClosureModel` builds a model from a closure, mainly for tests. It is planned
to move into a test-support crate, see [TODO.md](TODO.md).

## Linear Algebra

### faer for matrices

Matrices are [faer][faer] matrices, the version pinned to the one the
OpenModelica Rust runtime uses. faer is column-major like the C runtime, so C
buffers can be viewed without copying. Functions take views, `MatRef` and
`MatMut`, not owned matrices.

### Linear systems as `[J | f]`

A linear system is the `n×(n+1)` matrix `[J | f]`. The total pivot solver
solves `A·y = 0` with one entry of `y` fixed to 1: fixing the last column
solves `J·dy = -f`, letting the solver choose the column computes the tangent
of a homotopy path. `f` itself is stored in the last column.

### Permutations without moving data

The pivot search swaps rows and columns of a `PermutedMatrix`, a view that
only swaps entries of two index permutations. Algorithms use logical
positions, the data in memory never moves.

## Errors and Results

### Enums instead of status codes

Each crate has its own error enum, a successful result carries a status enum
where C has several return codes, for example `SolveStatus::Solved` and
`SolvedDeterminantZero`. Errors of a lower crate are wrapped by the crate
above with a `From` implementation, so `?` converts them.

An error variant only exists if a caller would act on it differently. The
line search backing off below `λ_min` is a `NewtonError::ModelError(Assert)`:
the model asserted at every point tried, and C reports all Newton failures
alike.

### Panics only for programming errors

Wrong shapes or indices passed by the caller panic. Bad values from a model,
like NaN in a Jacobian, are a normal failure in a simulation and return an
error, also because a panic must not unwind into the C runtime.

### No allocation in iterations

Buffers are allocated once and reused: workspace structs like
`NewtonWorkspace`, output slices and matrix views as parameters. Buffers of a
fixed size are `Box<[f64]>`, not `Vec`.

### Inputs first, outputs last

Like the standard library, for example `Read::read(&mut self, buf)`,
parameters that are only read come first and output slices last, as in
`ModelFunctions::eval(&mut self, x, res)` and
`scaling::variables(nominal, x, factor)`. A value that is modified in place,
like the matrix of `total_pivot_solver::solve` or `x` of `newton::solve`,
comes first, as `self` would. Small results are returned, not written into
parameters.

## Settings and Logging

### Settings with C defaults

Settings have the defaults of the C runtime, for example `newtonXTol`,
`newtonFTol` and `100·n` Newton iterations for `n` unknowns. They are set with
`with_*` methods, which check the values. Tolerances are stored as plain
values, as users know them; code that compares squared norms squares them
itself.

### `log` for messages

The libraries write messages with the [`log`][log] crate and don't choose a
logger. The target is the name of the C runtime's stream without the
`OMC_LOG_` prefix, for example `NLS_V` for the Newton iteration, so
`homotopy_capi` can forward them to the `OMC_LOG_*` streams. The Newton
iteration logs each step with `trace!` and its result with `debug!`.

## Tests

### Unit tests next to the code, integration tests for the public API

Unit tests are a `#[cfg(test)]` module at the end of the file they test, so
they can test private items. Tests of the public API of a crate are
integration tests in its `tests/` directory.

### Floating-point comparisons with named tolerances

Computed values are compared with a tolerance, using the comparators of faer
and [equator][equator] with a named function in the assertion, like
`assert!(a :approx_eq_eps(): b)`. The name says which tolerance is used:
`approx_eq_eps` for rounding errors, `approx_eq_fd_tol` for finite
differences.

## Project

### License

The code is part of the OpenModelica Run-Time System, distributed under the
BSD 3-Clause License, the AGPL version 3 or the OSMC-PL version 1.8, at the
user's choice. Every Rust file starts with the runtime license header, checked
by CI and added by the pre-commit hook.

### Checks

CI enforces the default formatting, Clippy and rustdoc without warnings, the
tests and the license headers. Versions of external dependencies are defined
once for the workspace.

[c-source]: https://github.com/OpenModelica/OpenModelica/blob/18d4ff455d1beceb4709d770ca4faf53d65732da/OMCompiler/SimulationRuntime/c/simulation/solver/nonlinearSolverHomotopy.c
[equator]: https://crates.io/crates/equator
[faer]: https://docs.rs/faer
[log]: https://crates.io/crates/log
