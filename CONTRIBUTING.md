# Contributing

## Setup

Install a stable Rust toolchain with [rustup][rustup].

### Git hook

The repository ships a pre-commit hook in [.githooks/](.githooks/) that formats
staged Rust files with `rustfmt` on every commit. Enable it once after cloning:

```sh
git config core.hooksPath .githooks
```

This only changes the configuration of your local clone. To disable the hook
again, run `git config --unset core.hooksPath`.

If a staged file also has unstaged changes, the hook only checks it, so that the
unstaged changes don't end up in the commit. When such a file is not formatted,
the commit is aborted. Run `cargo fmt`, stage the result and commit again.

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
which checks the formatting, builds the crate, runs Clippy, the tests and the
coverage, and builds the documentation. Pushes to `main` also publish the
documentation and the coverage report to [GitHub Pages][pages] with
[docs.yml](.github/workflows/docs.yml).

[cargo-llvm-cov]: https://github.com/taiki-e/cargo-llvm-cov
[pages]: https://anheuermann.github.io/homotopy_nls_solver/
[rustup]: https://rustup.rs/
