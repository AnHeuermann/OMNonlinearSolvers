# TODO

* Split DataHomotopy into multiple structs
  * Newton data
  * Homotopy parameters
  * Scaling vectors
  * Tolerances and settings
  * Statistics
* Check how to include hybrj_ from Minpack
* Move `max_number_of_iteartions` to configuration
* Replace our implementation of a newton algorithm with
 [newton_rootfinder][newton_rootfinder] or [nrfind][nrfind] (or similar)
* Move `ClosureModel` from `nls_model` into a test-support crate, e.g.
  `nls_test_problems`, used only as dev-dependency, so it is no longer part of
  the public API of `nls_model`. The crate can also collect standard test
  problems for nonlinear solvers (e.g. Rosenbrock's function, Powell's
  singular function, Broyden's tridiagonal function) and small systems from
  OpenModelica models, shared by the tests of `newton`, `homotopy_core` and
  regression tests against the C solver. Do this when the Newton tests need
  the same test systems in several crates.
* `SystemInfo` has public members, so a caller could input wrong values. Add
  getters.
* `total_pivot_solver::solve` allocates its solution and the permutations on
  every call, so each Newton iteration allocates. Let it write into buffers of
  the caller, as DESIGN.md asks.
* Statistics of a solve, e.g. the number of iterations the C runtime adds to
  `numberOfIterations`. Add them together with the homotopy.
* Maybe add the hybrd fallback of Minpack, which `solveHomotopy` tries when
  Newton fails, before the homotopy runs. Left out for now, there is no
  Minpack in Rust. It could also be a fallback solver passed in by
  `homotopy_capi`.
* Port corrector strategy 2 of `homotopyAlgorithm` (`homBacktraceStrategy`):
  back to the path orthogonally to the tangent, with the extended system of
  `orthogonalBacktraceMatrix`. Only strategy 1, the default, fixing one
  coordinate, is ported.
* Fix the scaling of the fixpoint homotopy's Jacobian once there is a test
  oracle against the C runtime of OpenModelica. C scales the columns of `J`
  before adding the identity, so the path following solves with
  `λ·J·S + (1 − λ)·I` instead of `(λ·J + (1 − λ)·I)·S`. `FixpointHomotopy`
  reproduces that, see its doc comment.
* Implement OpenModelica-like logging. The target is the C stream without the
  `OMC_LOG_` prefix and the `_V` suffix, and a verbose stream becomes `trace!`
  unless the message is a clear error. So
  `warningStreamPrint(OMC_LOG_NLS_V, ...)` usually becomes
  `trace!(target: "NLS", ...)`. Update the log calls in `newton` and
  `nls_model::scaling` (target `NLS_V` so far) and the `log` entry in
  DESIGN.md.

[newton_rootfinder]: https://github.com/Nateckert/newton_rootfinder
[nrfind]: https://github.com/NoraCodes/nrfind
