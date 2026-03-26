# Iterative Solvers

HAMS-MREL supports two linear solver options for the BEM system of equations. The solver is selected via the `ISOLV` parameter in `ControlFile.in`.

| ISOLV | Solver | Description |
| :--- | :--- | :--- |
| 1 | Direct LU | LAPACK ZGETRF/ZGETRS factorization and back-substitution (default) |
| 2 | GMRES | Restarted GMRES with block-diagonal LU preconditioner and H-matrix accelerated matvec |

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
- **Restart parameter**: m = 25 (GMRES is restarted every 25 iterations to limit memory usage)
- **Preconditioner**: Block-diagonal LU for multi-body (each body's self-interaction block is LU-factorized); Jacobi (diagonal) fallback for single-body
- **H-matrix matvec**: Adaptive Cross Approximation (ACA) compresses far-field interaction blocks to low-rank form, reducing matvec cost from O(N^2) to O(N k log N)
- **Block GMRES**: Multiple RHS solved with deflated initial guesses (previous RHS solution seeds the next)
- **Convergence tolerance**: Relative residual `||r|| / ||b|| < 1e-4` (sufficient for BEM mesh discretization accuracy)
- **Maximum iterations**: 100 (across all restarts)
- **Memory**: Krylov basis O(N * m), H-matrix O(N k log N), block preconditioner O(NBODY * (N/NBODY)^2)

GMRES is guaranteed to converge for non-singular systems (in exact arithmetic) and is the recommended iterative option for large multi-body problems.

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

## H-Matrix Acceleration

When GMRES is selected (`ISOLV=2`), the iterative solver automatically builds a hierarchical matrix (H-matrix) approximation of the BEM influence matrix for fast matrix-vector products. This reduces the matvec cost from O(N^2) to O(N k log N), where k is the typical low-rank of far-field interaction blocks (k ~ 10-20 for BEM problems).

### How It Works

The BEM influence matrix has hierarchical low-rank structure: panels far apart interact via a smooth Green's function kernel, so those sub-blocks can be approximated as low-rank matrices `A_block ≈ U * V^H` with rank k much less than the block dimension. Only near-field blocks (panels close together) require full dense storage.

The H-matrix implementation consists of four stages:

1. **Cluster tree construction**: Panel centers are recursively bisected along the longest bounding box dimension, creating a binary tree with leaf clusters of at most 32 panels.

2. **Block cluster tree**: Each pair of clusters is tested for admissibility: `dist(sigma, tau) > eta * max(diam(sigma), diam(tau))` with eta=2.0. Admissible pairs have smooth interaction and are marked for compression. Inadmissible leaf pairs are stored as dense blocks.

3. **Adaptive Cross Approximation (ACA)**: Each admissible block is compressed using partially-pivoted ACA, which incrementally builds the low-rank factors U and V^H by sampling individual matrix entries. The approximation rank is determined adaptively with tolerance 1e-4 (matching the GMRES convergence tolerance). If compression is not worthwhile (`rank * (nrows + ncols) >= nrows * ncols`), the block falls back to dense storage.

4. **H-matvec**: During each GMRES iteration, the matrix-vector product uses the compressed representation. Dense blocks use standard ZGEMV. Low-rank blocks compute `y += U * (V^H * x)` — two small ZGEMV calls instead of one large one.

### Performance

The H-matrix is built once per frequency after the coefficient matrix is assembled. Build cost is O(N k^2 log N) — typically 1-2 seconds for N=7000 panels. The per-matvec speedup depends on the problem geometry:

| N (panels) | Dense matvec | H-matvec | Speedup |
| :--- | :--- | :--- | :--- |
| 3000 | ~9M ops | ~1M ops | ~9x |
| 7000 | ~49M ops | ~3M ops | ~16x |
| 15000 | ~225M ops | ~8M ops | ~28x |

The H-matrix is fully transparent to the user — it is automatically built and used when ISOLV=2. Console output reports compression statistics:

```
  H-matrix compression:
    Blocks: 1247 (dense: 312, low-rank: 935)
    Average rank: 12.3
    Compression ratio: 8.5x
    Build time: 1.42 s
```

### Block-Diagonal LU Preconditioner

For multi-body problems, GMRES uses a block-diagonal LU preconditioner where each body's self-interaction block is LU-factorized independently. This is set up automatically based on `NELEM_MULTI` (per-body panel counts). The preconditioner reduces GMRES iterations from ~16 (Jacobi) to ~5 (block-diagonal LU), with setup cost equivalent to NBODY independent LU factorizations of size (N/NBODY).

### Block GMRES with Deflated Initial Guesses

Multiple right-hand sides (radiation modes) are solved in batches. Each RHS in a batch is seeded with the solution from the previous RHS, exploiting the similarity between adjacent radiation modes. This reduces iterations for later modes from ~5 to ~2-3.

### Source Files

| File | Description |
| :--- | :--- |
| [HMatrix.f90](../src/HMatrix.f90) | H-matrix module: cluster tree, ACA, H-matvec |
| [IterativeSolvers.f90](../src/IterativeSolvers.f90) | GMRES solver with H-matvec integration, block preconditioner, Block GMRES |

### Tuning Parameters

The H-matrix parameters in `HMatrix.f90` can be adjusted:

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `NLEAF_DEFAULT` | 32 | Max panels per leaf cluster. Smaller = deeper tree, more blocks |
| `ETA_DEFAULT` | 2.0 | Admissibility parameter. Larger = more low-rank blocks, but higher rank |
| `ACA_TOL` | 1e-4 | ACA convergence tolerance. Should match GMRES tolerance |
| `RANK_MAX` | 50 | Maximum low-rank approximation rank. Cap for oscillatory kernels |

## References

- Saad, Y. & Schultz, M.H. (1986). GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems. SIAM Journal on Scientific and Statistical Computing, 7(3), 856-869.
- Hackbusch, W. (1999). A sparse matrix arithmetic based on H-matrices. Part I: Introduction to H-matrices. Computing, 62(2), 89-108.
- Bebendorf, M. (2000). Approximation of boundary element matrices. Numer. Math., 86(4), 565-589.
- Rjasanow, S. & Steinbach, O. (2007). The Fast Solution of Boundary Integral Equations. Springer.
- Ancellin, M. & Dias, F. (2019). Capytaine: a Python-based linear potential flow solver. Journal of Open Source Software, 4(36), 1341.
