# The C tests of the C interface of homotopy_capi, for CMake and for IDEs that
# open the repository root. The Rust crates are built with cargo, see
# CONTRIBUTING.md.

cmake_minimum_required(VERSION 3.16)
project(homotopy_nls_solver C)

# Like enable_testing(), and writes the DartConfiguration.tcl that `ctest -T`
# of IDEs reads
include(CTest)
add_subdirectory(homotopy_capi/tests/c)
