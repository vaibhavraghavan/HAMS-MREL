# Development Documentation

This page contains information for the developers and mantainers of HAMS-MREL:

- [Tests](https://github.com/vaibhavraghavan/HAMS-MREL/edit/develop/docs/developer_documentation.md#tests)
- [Continous Integration](https://github.com/vaibhavraghavan/HAMS-MREL/edit/develop/docs/developer_documentation.md#continous-integration-in-github-actions)
- [Profiling Results](https://github.com/vaibhavraghavan/HAMS-MREL/edit/develop/docs/developer_documentation.md#profiling-results)

## Tests

HAMS-MREL includes both unit tests and integration tests written using the [Test-Drive](https://github.com/fortran-lang/test-drive) framework. The tests are designed to ensure that both individual components (subroutines, functions) and full simulation workflows are correct.

The tests are defined in [fpm.toml](../fpm.toml) as two separate executables: unit-tester and integration-tester. The unit-tester runs small-scale tests on individual subroutines and functions. The integration-tester runs full simulations for pre-defined test configurations and compares outputs with reference results. All of the data needed for the tests in included in the tests/data folder.

### Running the tests

HAMS-MREL can be built using multiple methods, but tests must always be run via **fpm** because they rely on the fpm build/test environment configured in `fpm.toml`. fpm will automatically install and link the Test-Drive framework to the test executables.

```bash
# FPM: set the compiler to IFX
export FPM_FC=IFX

# FPM: set the compiler flags
# `-heap-arrays 0` is required to keep large multi-body cases from overflowing the per-thread
# stack — see the installation docs for the full explanation.
export FPM_FFLAGS="-qmkl -qopenmp -heap-arrays 0"
export FPM_LDFLAGS="-qmkl -qopenmp"

# Compilation step. Build is handled by FPM (HAMS-MREL and test executables will be placed in build directory)
fpm build

# Place the executable in 'install' directory (an exact locaion is needed by integration tests)
fpm install --prefix install

# Run unit tests
fpm test --target unit-tester

# Run integration tests (they expect the path to the executable as an argument)
fpm test --target integration-tester -- install/bin/hamsmrel
```

After executing these commands you should see the following output (or similar):

<img width="557" height="280" alt="hams-tests" src="https://github.com/user-attachments/assets/1bb39a74-593d-498b-8a6d-40f20fcc36f9" />

## Continous Integration in GitHub Actions

The project has [Continuous Integration (CI)](https://www.software.ac.uk/guide/how-continuous-integration-can-help-you-regularly-test-and-release-your-software) workflows in place to verify the program compiles and runs correctly. The workflows are configured with [GitHub Actions](https://docs.github.com/en/actions). They automatically compile HAMS-MREL and run the tests.

There are two workflows: one for Windows and another for Linux. Both run independently of each other and results can be seen in the [Actions tab](https://github.com/vaibhavraghavan/HAMS-MREL-DCC/actions) of the GitHub repository. The workflows are automatically triggered whenever a new commit is pushed to the repository (main or develop branch) or a pull request is opened. Running the workflows when a pull request is opened helps us verify that the changes that will be integrated compile correctly and don't cause tests to fail.

Both workflows define a single job each: *Build and Test*. The job:
- Configures the GitHub runner (server) where the workflow is executed. This means installing the Intel compilers and Intel OneAPI Math Kernel Library, Fortran Package Manager and Test-Drive framework.
- Compiles HAMS-MREL
- Runs the tests

To simplify the installation of the dependencies mentioned above, we use Actions from the [GitHub Actions marketplace](https://github.com/marketplace?type=actions). Actions can be thought of as mini-apps, extensions or plugins that do one thing specifically. In this case, the Actions install different software on the GitHub Actions runners. The following Actions are used:
- [setup-fortran](https://github.com/marketplace/actions/setup-fortran): install the Fortran Intel compiler and Intel OneAPI Math Kernel Library (Linux only)
  Although the setup-fortran action can install the Intel Fortran compiler on both Windows and Linux, installation of the oneMKL library is only available for Linux systems. Intel dependencies are manually installed on the Windows runner.
- [setup-fpm](https://github.com/marketplace/actions/setup-fpm): install the Fortran Package Manager

**The Windows workflow configuration is defined in [build_and_test_windows.yml](../.github/build_and_test_windows.yml).**

**The configuration for Linux is in [build_and_test_linux.yml](../.github/build_and_test_linux.yml).**

> [!NOTE]
> **Windows runner pinned to `windows-2022`.** Around mid-June 2026, GitHub migrated
> `windows-latest` to `windows-2025` / `windows-2025-vs2026`. The new image places
> Git for Windows' coreutils `link.exe` on PATH ahead of MSVC's linker, which breaks
> fpm + ifx's link step (`link: unknown option -- s`). The workflow is therefore
> pinned to `windows-2022` until upstream fpm or setup-fortran handles the new image
> cleanly. When you eventually migrate, the build step will need to strip
> `C:\Program Files\Git\usr\bin` (and possibly `C:\Program Files\Git\bin`) from
> PATH before `fpm build`, OR you'll need to use a different invocation that
> ensures MSVC tools are resolved first.

### About GitHub Actions

A GitHub Actions workflow is defined in a `<workflow-name>.yml` file. There can be multiple workflows, each one defined in its own file. All workflow files must be located in a `.github/workflows` folder, at the top level of the repository. 

Workflow files are written in YAML. [See the GitHub specification for workflow files](https://docs.github.com/en/actions/writing-workflows/workflow-syntax-for-github-actions) (**important link**). 

Workflows run on virtual machines provided by GitHub. These machines are called `runners`. GitHub allows us to choose among different pre-configured runners, each with a different operating system and hardware. All runners come with pre-installed compilers, tools and software. Unfortunately, none of the HAMS-MREL dependencies come by default with the runners so we must install them ourselves. Runners are created when a workflow is 

Workflows define `triggers` that tell GitHub when the workflow needs to be executed. A workflow defines one or more `jobs`, and each job is made up of one or more `steps`. Steps within a job are executed sequentially. By default, all jobs run in parallel unless dependencies are defined between them.

## Iterative Solvers

HAMS-MREL supports two linear-solver paths for the BEM system, dispatched at runtime by the `ISOLV` parameter in `ControlFile.in` (1 = direct LU, default; 2 = GMRES). The full user-facing description is in [iterative_solvers.md](iterative_solvers.md); this section covers the developer-side layout.

**Source files:**

| File | Contents |
|---|---|
| [`src/IterativeSolvers.f90`](../src/IterativeSolvers.f90) | Restarted GMRES(m), block-diagonal LU preconditioner, Block GMRES with deflated initial guesses. Public entry points: `ZGMRES_SOLVE`, `ZBLOCK_GMRES_SOLVE`, `ZITER_SOLVE_MULTI`, `SETUP_BLOCK_PRECONDITIONER`. Uses inline BLAS (`ZDOTC`/`DZNRM2`/`ZGEMV` unrolled) to sidestep known ifort/MKL calling-convention issues with the compiled BLAS bindings. |
| [`src/HMatrix.f90`](../src/HMatrix.f90) | Hierarchical-matrix compression: cluster tree, block cluster tree with admissibility test, Adaptive Cross Approximation (ACA), fast H-matvec. Public entry points: `HMATRIX_BUILD`, `HMATVEC`, `HMATRIX_DESTROY`. |

**Dispatch pattern:**

Each of the four `AssbMatx*.f90` files gates the LAPACK path on `ISOLV`. `ASSB_LEFT*` skips the `ZGETRF` factorization when `ISOLV > 1`. Each `RADIATION_SOLVER*` / `DIFFRACTION_SOLVER*` branches: `ISOLV = 1` uses the existing in-place `ZGETRS` path (Phase 2.4 — passes `AMAT`/`CMAT` directly, no full-size workspace copy); `ISOLV = 2` allocates a scratch matrix, builds the H-matrix + block preconditioner, and calls `ZITER_SOLVE_MULTI`.

**Interaction with `GRN_KIND=4`:**

The `ACA_TOL` and `ITER_TOL` constants are both `1e-4`, matched to the ~`1e-7` per-entry precision floor of an `AMAT` assembled from single-precision Green's-function values. Iterative results agree with the direct-LU reference at 4–5 significant figures on physically meaningful entries. For bit-identity against WAMIT references, users should stay on the default (`GRN_KIND = 8`, `ISOLV = 1`).

**Adding a new iterative test:**

Small multi-body test cases with `ISOLV = 2` can be added to `test/integration_tester.f90`. The recommended check is that `ISOLV = 2` output agrees with an `ISOLV = 1` baseline to `1e-3` relative on the WAMIT `.1` (`AmssDamp.1`) file peak-magnitude entries.

## Profiling Results

This section documents the primary sources of memory consumption in HAMS-MREL and how they scale with problem size. It is intended to guide future development and optimization efforts aimed at reducing RAM usage.

Profiling of HAMS-MREL shows that memory usage is dominated by a small number of large, multi-dimensional arrays whose sizes scale with the number of bodies and the number of panels per body. These arrays are dynamically allocated before the main solver iteration loop begins and remain allocated for the duration of the simulation.

The table below lists the highest-memory-consumption variables and their sizes for systems involving 2, 3, and 6 bodies. For all cases shown, each body is discretized using 234 panels. For all cases shown, each body is discretized using 234 panels.

| System | Variable Name | Type | Size | Number of Elements |
| :--- | :--- | :--- | :--- | :--- |
| 2-bodies (234 panels each) | `CGRN_MULTI_COMB` & `RKBN_MULTI_COMB` | `COMPLEX(8)` | (1:468,1:468,1,4) | 876096 |
| 2-bodies (234 panels each) | `AMAT_MULTI` | `COMPLEX(8)` | (1:468,1:468,1:1) | 219024 |
| 3-bodies (234 panels each) | `CGRN_MULTI_COMB` & `RKBN_MULTI_COMB` | `COMPLEX(8)` | (1:702,1:702,1:1,1:4) | 1971216 |
| 3-bodies (234 panels each) | `AMAT_MULTI` | `COMPLEX(8)` | (1:702,1:702,1:1) | 492804 |
| 6-bodies (234 panels each) | `CGRN_MULTI_COMB` & `RKBN_MULTI_COMB` | `COMPLEX(8)` | (1:1404,1:1404,1:1,1:4) | 7884864 |
| 6-bodies (234 panels each) | `AMAT_MULTI` | `COMPLEX(8)` | (1:1404,1:1404,1:1) | 1971216 |

### Observations

- **Dominant memory consumers**
  - The first two dimensions of these arrays scale linearly with the total number of panels across all bodies, resulting in a quadratic growth in memory usage as the number of bodies increases.
  - The variables are allocated once per program execution, before the solver iteration begins.
- **Overall memory behavior**
  - For a given simulation, memory usage is largely constant throughout execution.
  - Peak memory consumption is approximately 2 GB for the cases studied.
  - Memory does not accumulate across solver iterations.



