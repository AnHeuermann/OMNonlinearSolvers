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

[newton_rootfinder]: https://github.com/Nateckert/newton_rootfinder
[nrfind]: https://github.com/NoraCodes/nrfind
