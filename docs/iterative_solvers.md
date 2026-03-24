# Iterative Solvers

HAMS-MREL supports three linear solver options for the BEM system of equations. The solver is selected via the `ISOLV` parameter in `ControlFile.in`.

| ISOLV | Solver | Description |
| :--- | :--- | :--- |
| 1 | Direct LU | LAPACK ZGETRF/ZGETRS factorization and back-substitution (default) |
| 2 | GMRES | Restarted Generalized Minimal Residual method with Jacobi preconditioner |
| 3 | BiCGSTAB | Bi-Conjugate Gradient Stabilized method with Jacobi preconditioner |

When `ISOLV` is omitted from `ControlFile.in`, the solver defaults to direct LU (ISOLV=1), preserving full backward compatibility.

## Background

The BEM formulation in HAMS-MREL produces a dense complex linear system `A * x = b` at each wave frequency, where `A` is the influence coefficient matrix of size `NELEM x NELEM` (number of body panels squared). The default direct LU solver has O(N^3) computational cost for the factorization step, which becomes the dominant bottleneck for large panel counts.

Iterative solvers offer an O(N^2 * k) alternative, where k is the number of iterations required for convergence. For well-conditioned BEM systems, k is typically much smaller than N, yielding significant speedups for large meshes (N > 2000 panels).

## Solver Details

### Direct LU (ISOLV=1)

The default solver uses LAPACK routines:
- **ZGETRF**: LU factorization with partial pivoting of the complex coefficient matrix
- **ZGETRS**: Triangular back-substitution to solve for the radiation (6 or 6*NBODY right-hand sides) and diffraction (1 right-hand side) potentials

The LU factorization is computed once per frequency in the `ASSB_LEFT` subroutine and reused for both the radiation and diffraction solves. This is the most robust option and is recommended for small to moderate panel counts.

### GMRES (ISOLV=2)

The Generalized Minimal Residual method (Saad & Schultz, 1986) builds an orthonormal basis for the Krylov subspace and minimizes the residual norm over that subspace.

**Implementation details:**
- **Restart parameter**: m = 50 (GMRES is restarted every 50 iterations to limit memory usage)
- **Preconditioner**: Jacobi (diagonal) preconditioner using `M = diag(A)`. The diagonal of the BEM influence matrix contains the `2*pi` self-influence term, providing effective scaling.
- **Convergence tolerance**: Relative residual `||r|| / ||b|| < 1e-7`
- **Maximum iterations**: 500 (across all restarts)
- **Memory**: Requires storage for the Krylov basis vectors: O(N * m) complex values per right-hand side

GMRES is guaranteed to converge for non-singular systems (in exact arithmetic) and is generally the more robust iterative option.

### BiCGSTAB (ISOLV=3)

The Bi-Conjugate Gradient Stabilized method (van der Vorst, 1992) is a Krylov method for non-symmetric systems that avoids the irregular convergence behavior of standard Bi-CG.

**Implementation details:**
- **Preconditioner**: Jacobi (diagonal), same as GMRES
- **Convergence tolerance**: `||r|| / ||b|| < 1e-7`
- **Maximum iterations**: 500
- **Memory**: Requires only O(N) work vectors per right-hand side (lower than GMRES)

BiCGSTAB uses less memory than GMRES but can occasionally exhibit breakdown or stagnation for ill-conditioned systems.

## Enabling Iterative Solvers

### ControlFile.in

Add the following optional line immediately after `Number of threads` in the Solver Options section:

```
    Solver_type (ISOLV):     2
```

| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| (optional) | `26x,i16` | `ISOLV` | Linear solver type: 1 = Direct LU, 2 = GMRES, 3 = BiCGSTAB |

If this line is absent, ISOLV defaults to 1 (direct LU). The line is detected by the presence of the keyword `Solver_type`.

**Example ControlFile.in excerpt:**
```
    Reference length          1.000000000000000
    Solution type (ISOL):    1
    If_remove_irr_freq        0
    Number of threads        32
    Solver_type (ISOLV):      2
```

The selected solver type is printed to the console at startup:
```
 Linear solver type:             GMRES
```

## Scope

The iterative solver option applies to all four solver variants in HAMS-MREL:

| File | Variant | Subroutines Modified |
| :--- | :--- | :--- |
| [AssbMatx.f90](../src/AssbMatx.f90) | Single-body, no IRR | ASSB_LEFT, RADIATION_SOLVER, DIFFRACTION_SOLVER |
| [AssbMatxMulti.f90](../src/AssbMatxMulti.f90) | Multi-body, no IRR | ASSB_LEFT_MULTI, RADIATION_SOLVER_MULTI, DIFFRACTION_SOLVER_MULTI |
| [AssbMatx_irr.f90](../src/AssbMatx_irr.f90) | Single-body, with IRR | ASSB_LEFT_IRR, RADIATION_SOLVER_IRR, DIFFRACTION_SOLVER_IRR |
| [AssbMatx_irrMulti.f90](../src/AssbMatx_irrMulti.f90) | Multi-body, with IRR | ASSB_LEFT_IRR_MULTI, RADIATION_SOLVER_IRR_MULTI, DIFFRACTION_SOLVER_IRR_MULTI |

For all variants:
- When `ISOLV=1`: the coefficient matrix is LU-factorized in `ASSB_LEFT*`, and `ZGETRS` is called in the solver subroutines (existing behavior).
- When `ISOLV>1`: the coefficient matrix is preserved unfactored (LU factorization is skipped), and the iterative solver is called in the solver subroutines with the unfactored matrix for matrix-vector products.

For the irregular frequency removal variants (`*_irr*`), the iterative solver operates on the normal equations matrix `C = A^T * A`, which is the system actually solved by these routines.

## Source Files

| File | Description |
| :--- | :--- |
| [IterativeSolvers.f90](../src/IterativeSolvers.f90) | Self-contained GMRES and BiCGSTAB implementations for COMPLEX*16 systems. Uses only standard BLAS routines (ZGEMV, ZDOTC, DZNRM2). No additional dependencies. |
| [WavDynMods.f90](../src/WavDynMods.f90) | `ISOLV` declared in `HAMS_mod` alongside existing solver parameters |
| [InputFiles.f90](../src/InputFiles.f90) | `ISOLV` parsed from `ControlFile.in` in `ReadControlFile` |

## Performance Guidance

| Panel Count | Recommended Solver | Rationale |
| :--- | :--- | :--- |
| < 1000 | Direct LU (ISOLV=1) | LU is fast for small systems; overhead of iterative setup not worthwhile |
| 1000 - 3000 | Either | Performance is comparable; GMRES may be slightly faster |
| > 3000 | GMRES (ISOLV=2) | O(N^2 * k) iterative cost significantly beats O(N^3) direct factorization |

**Notes:**
- The iterative solvers solve each right-hand side independently. For the radiation problem, this means 6 (or 6*NBODY for multi-body) independent solves per frequency.
- If the iterative solver fails to converge, a warning is printed with the final residual norm. The computation continues with the approximate solution rather than aborting.
- At irregular frequencies (when `IRSP=0`), the BEM matrix becomes singular and iterative solvers will fail to converge. Use `IRSP=1` (irregular frequency removal) when using iterative solvers.

## Convergence Troubleshooting

If the iterative solver reports convergence warnings:

1. **Switch to direct LU**: Set `ISOLV=1` for reliable results at any frequency.
2. **Enable irregular frequency removal**: Set `IRSP=1` to improve matrix conditioning.
3. **Try the other iterative solver**: GMRES and BiCGSTAB have different convergence characteristics. If one fails, the other may succeed.
4. **Check mesh quality**: Poorly shaped panels (very elongated or with near-zero area) degrade matrix conditioning.

## References

- Saad, Y. & Schultz, M.H. (1986). GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems. SIAM Journal on Scientific and Statistical Computing, 7(3), 856-869.
- van der Vorst, H.A. (1992). Bi-CGSTAB: A fast and smoothly converging variant of Bi-CG for the solution of nonsymmetric linear systems. SIAM Journal on Scientific and Statistical Computing, 13(2), 631-644.
- Ancellin, M. & Dias, F. (2019). Capytaine: a Python-based linear potential flow solver. Journal of Open Source Software, 4(36), 1341.
