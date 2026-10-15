# Contributing

## Setup

Install a stable Rust toolchain with [rustup][rustup].

### Git hook

The repository ships a pre-commit hook in [.githooks/](.githooks/) that runs on
every commit. For staged Rust files it adds a missing license header, updates an
outdated copyright year and formats the file with `rustfmt`. Enable it once
after cloning:

```sh
git config core.hooksPath .githooks
```

This only changes the configuration of your local clone. To disable the hook
again, run `git config --unset core.hooksPath`.

If a staged file also has unstaged changes, the hook only checks it, so that the
unstaged changes don't end up in the commit. When such a file has no valid
license header or is not formatted, the commit is aborted. Run
`.ci/check_license.py --fix` and `cargo fmt`, stage the result and commit again.

## Build

```sh
cargo build
```

## Test

```sh
cargo test
```

## Formatting and Linting

The code uses the default `rustfmt` style and is checked with Clippy:

```sh
cargo fmt
cargo clippy --all-targets
```

## License Headers

Every Rust file must start with the OpenModelica Run-Time System license header,
with the current year in the copyright line. Add missing headers and update the
year with:

```sh
.ci/check_license.py --fix
```

Without `--fix`, the script only checks the files. Pass file names to check only
those files instead of all Rust files in the repository.

## Coverage

Coverage is measured with [cargo-llvm-cov][cargo-llvm-cov]. Install it once:

```sh
rustup component add llvm-tools-preview
cargo install cargo-llvm-cov --locked
```

Run the tests with coverage and open the HTML report:

```sh
cargo llvm-cov --open
```

## Documentation

Build the documentation, including private items, and open it in the browser:

```sh
cargo doc --no-deps --document-private-items --open
```

## Continuous Integration

Every pull request and push to `main` runs [ci.yml](.github/workflows/ci.yml),
which checks the formatting and license headers, builds the crate, runs Clippy
and the tests, and builds the documentation. The coverage job prints a summary
on the run's page and uploads the HTML report as the `coverage` artifact.

[cargo-llvm-cov]: https://github.com/taiki-e/cargo-llvm-cov
[rustup]: https://rustup.rs/
