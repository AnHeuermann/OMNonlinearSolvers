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
