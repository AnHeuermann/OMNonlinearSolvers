# Instructions for coding agents

This file is for AI coding agents (and the humans steering them) working on
this repository. [CONTRIBUTING.md](CONTRIBUTING.md) describes the setup, the
checks and CI, [DESIGN.md](DESIGN.md) the design decisions. Both apply to
agent-written changes too.

## The project

The project ports `nonlinearSolverHomotopy.c` of the OpenModelica C runtime to
Rust. It shall replace the C solver and the Rust port in OpenModelica, so it
reproduces the C algorithm, see "Same behaviour as the C solver" in
DESIGN.md.

If you are unsure about a design decision, ask the human you work with
instead of guessing.

## Clean code, but not dogmatic

Follow the principles of clean code, for example meaningful names, small
functions with one job, no duplication, and code that says what it does
without comments. They are a means to readable code, not rules in their own
right:

* Prefer the version that is easier to read over the one that follows a rule
  more strictly.
* Don't split a function or add an abstraction only to satisfy a principle.
  Three similar lines can be clearer than a helper for them.
* Don't introduce structs, traits or generics before there is a second user,
  unless they prevent real mistakes, such as swapping two `f64` arguments.

## Conventions

* The decisions in DESIGN.md apply, for example inputs first and outputs last,
  no allocation in iterations, and error enums instead of status codes.
  Update DESIGN.md in the same commit when a decision changes or a new one is
  made.
* Keep the constants, thresholds and order of checks of the C code. Document
  a deliberate deviation where it is made, and add it to DESIGN.md.
* Every item has a doc comment. Name the C function or variable a port comes
  from, so the two can be compared.
* Comments in code explain only what the code can't say itself.
* Unit tests live next to the code, integration tests of the public API in
  `tests/`. Derive expected values by hand or from the C algorithm and write
  the derivation into the test's doc comment.

## Before committing

Run the checks of [CONTRIBUTING.md](CONTRIBUTING.md): `cargo fmt`, Clippy with
`-D warnings`, `cargo test`, `cargo doc` and `.ci/check_license.py`. CI fails
on any warning.

## Git

* Commit or push only when asked.
* Stage files explicitly and commit only files that are tracked or belong to
  the change. Untracked files in the working tree may be the human's private
  notes.
* The first line of a commit message is a summary of at most 50 characters,
  the body explains what changed and why, wrapped at 72 characters.
* Say in the commit message which parts an agent wrote, and end it with a
  `Co-Authored-By` trailer naming the model.
