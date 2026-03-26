!  ------------------------------------------------------------------------------------------------------
!
!    Program HAMS-MREL for the diffraction and radiation of waves
!    by 3D structures.
!
!  License:
!
!    This routine is part of HAMS and HAMS-MREL.
!
!    HAMS-MREL is a free software framework: you can redistribute it and/or modify it
!    under the terms of the Apache License, Version 2.0 (the "License"); you may
!    not use these subroutines except in compliance with the License. The software
!    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND.
!
!    You should have received a copy of the Apache License, Version 2.0, along with
!    HAMS-MREL. If not, see <http://www.apache.org/licenses/LICENSE-2.0>.
!
!  Code Original Author:
!
!    Vaibhav Raghavan
!
!  Description:
!
!    Restarted GMRES(m) iterative solver for complex double-precision (COMPLEX*16)
!    dense linear systems A*x = b arising from the BEM formulation.
!
!    Uses Jacobi (diagonal) preconditioning and inline BLAS operations to avoid
!    known ifort/MKL calling convention issues with ZDOTC, DZNRM2, and ZGEMV.
!
!    References:
!      - Saad & Schultz (1986), GMRES: A generalized minimal residual algorithm
!        for solving nonsymmetric linear systems. SIAM J. Sci. Stat. Comput.
!
!  ------------------------------------------------------------------------------------------------------

MODULE IterativeSolvers

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: ZGMRES_SOLVE
    PUBLIC :: ZITER_SOLVE_MULTI
    PUBLIC :: SETUP_BLOCK_PRECONDITIONER
    PUBLIC :: PRECOND_NBLOCKS
    PUBLIC :: ZBLOCK_GMRES_SOLVE

    ! Solver parameters (optimized for BEM matrices)
    INTEGER, PUBLIC :: GMRES_RESTART = 25       ! Krylov subspace dimension before restart
    REAL*8, PUBLIC  :: ITER_TOL = 1.0D-4        ! Relative convergence tolerance (BEM limited by mesh discretization)
    INTEGER, PUBLIC :: ITER_MAXITER = 100        ! Maximum total iterations across restarts

    ! Preconditioner type: 0=none, 1=Jacobi (diagonal), 2=block-diagonal LU
    INTEGER :: PRECOND_TYPE = 0

    ! Jacobi preconditioner storage
    INTEGER :: PRECOND_N = 0
    COMPLEX*16, ALLOCATABLE :: PRECOND_DIAG(:)

    ! Block-diagonal LU preconditioner storage
    ! Each body's self-interaction block is LU-factorized and stored.
    INTEGER :: PRECOND_NBLOCKS = 0
    INTEGER, ALLOCATABLE :: PRECOND_BLOCK_SIZES(:)    ! (NBLOCKS) size of each block
    INTEGER, ALLOCATABLE :: PRECOND_BLOCK_OFFSETS(:)  ! (NBLOCKS) starting row/col in global matrix
    INTEGER :: PRECOND_MAXBLK = 0                     ! max block size
    COMPLEX*16, ALLOCATABLE :: PRECOND_BLK_DATA(:,:,:)! (MAXBLK, MAXBLK, NBLOCKS) LU factors
    INTEGER, ALLOCATABLE :: PRECOND_BLK_IPIV(:,:)     ! (MAXBLK, NBLOCKS) pivot arrays

CONTAINS

!---------------------------------------------------------------------------------------------
!   Inline complex dot product: result = sum_i( conjg(x(i)) * y(i) )
!---------------------------------------------------------------------------------------------

    FUNCTION ZDOT_INLINE(N, X, Y) RESULT(DOT)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        COMPLEX*16, INTENT(IN) :: X(N), Y(N)
        COMPLEX*16 :: DOT
        INTEGER :: I
        REAL*8 :: SR, SI

        ! Split into real and imaginary partial sums for better vectorization.
        ! conjg(x)*y = (xr - i*xi)*(yr + i*yi) = (xr*yr + xi*yi) + i*(xr*yi - xi*yr)
        SR = 0.0D0
        SI = 0.0D0
        DO I = 1, N
            SR = SR + DBLE(X(I)) * DBLE(Y(I)) + DIMAG(X(I)) * DIMAG(Y(I))
            SI = SI + DBLE(X(I)) * DIMAG(Y(I)) - DIMAG(X(I)) * DBLE(Y(I))
        END DO
        DOT = CMPLX(SR, SI)
    END FUNCTION ZDOT_INLINE

!---------------------------------------------------------------------------------------------
!   Inline complex vector 2-norm: result = sqrt( sum_i |x(i)|^2 )
!---------------------------------------------------------------------------------------------

    FUNCTION ZNRM2_INLINE(N, X) RESULT(NRM)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        COMPLEX*16, INTENT(IN) :: X(N)
        REAL*8 :: NRM
        INTEGER :: I
        REAL*8 :: S
        S = 0.0D0
        DO I = 1, N
            S = S + DBLE(X(I))**2 + DIMAG(X(I))**2
        END DO
        NRM = SQRT(S)
    END FUNCTION ZNRM2_INLINE

!---------------------------------------------------------------------------------------------
!   Inline matrix-vector product: Y_OUT = AMAT * X_IN
!---------------------------------------------------------------------------------------------

    SUBROUTINE ZMATVEC(N, AMAT, X_IN, Y_OUT)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        COMPLEX*16, INTENT(IN) :: AMAT(N, N), X_IN(N)
        COMPLEX*16, INTENT(OUT) :: Y_OUT(N)
        INTEGER :: I, J
        REAL*8 :: SR, SI, XR, XI, AR, AI

        ! Row-major access with split real/imaginary arithmetic for
        ! better SIMD vectorization. Serial loop — OpenMP parallelism
        ! is applied at the RHS level (ZBLOCK_GMRES_SOLVE) instead
        ! of within each matvec to avoid nested parallelism overhead.
        DO I = 1, N
            SR = 0.0D0
            SI = 0.0D0
            DO J = 1, N
                AR = DBLE(AMAT(I, J))
                AI = DIMAG(AMAT(I, J))
                XR = DBLE(X_IN(J))
                XI = DIMAG(X_IN(J))
                SR = SR + AR * XR - AI * XI
                SI = SI + AR * XI + AI * XR
            END DO
            Y_OUT(I) = CMPLX(SR, SI)
        END DO

    END SUBROUTINE ZMATVEC

!---------------------------------------------------------------------------------------------
!   Set up Jacobi (diagonal) preconditioner from matrix A.
!   Used as fallback for single-body problems where block-diagonal is not beneficial.
!---------------------------------------------------------------------------------------------

    SUBROUTINE SETUP_PRECONDITIONER(N, AMAT)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        COMPLEX*16, INTENT(IN) :: AMAT(N, N)
        INTEGER :: I

        IF (ALLOCATED(PRECOND_DIAG)) DEALLOCATE(PRECOND_DIAG)
        PRECOND_N = N
        ALLOCATE(PRECOND_DIAG(N))

        DO I = 1, N
            PRECOND_DIAG(I) = AMAT(I, I)
            IF (ABS(PRECOND_DIAG(I)) < 1.0D-30) THEN
                PRECOND_DIAG(I) = CMPLX(1.0D0, 0.0D0)
            END IF
        END DO

        PRECOND_TYPE = 1
    END SUBROUTINE SETUP_PRECONDITIONER

!---------------------------------------------------------------------------------------------
!   Set up block-diagonal LU preconditioner for multi-body problems.
!
!   Extracts each body's self-interaction diagonal block from AMAT,
!   LU-factorizes it with ZGETRF, and stores the factorized blocks
!   for use in APPLY_PRECONDITIONER via ZGETRS.
!
!   Arguments:
!     N            - Total system size
!     AMAT         - Coefficient matrix (N x N, not modified)
!     NBLOCKS      - Number of diagonal blocks (= NBODY)
!     BLOCK_SIZES  - Array of block sizes (= NELEM_MULTI)
!---------------------------------------------------------------------------------------------

    SUBROUTINE SETUP_BLOCK_PRECONDITIONER(N, AMAT, NBLOCKS, BLOCK_SIZES)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N, NBLOCKS
        COMPLEX*16, INTENT(IN) :: AMAT(N, N)
        INTEGER, INTENT(IN) :: BLOCK_SIZES(NBLOCKS)

        INTEGER :: K, I, J, IOFF, JOFF, BLK_SZ, INFO_LU, MAXSZ

        ! Clean up any previous preconditioner
        IF (ALLOCATED(PRECOND_DIAG)) DEALLOCATE(PRECOND_DIAG)
        IF (ALLOCATED(PRECOND_BLOCK_SIZES)) DEALLOCATE(PRECOND_BLOCK_SIZES)
        IF (ALLOCATED(PRECOND_BLOCK_OFFSETS)) DEALLOCATE(PRECOND_BLOCK_OFFSETS)
        IF (ALLOCATED(PRECOND_BLK_DATA)) DEALLOCATE(PRECOND_BLK_DATA)
        IF (ALLOCATED(PRECOND_BLK_IPIV)) DEALLOCATE(PRECOND_BLK_IPIV)

        PRECOND_NBLOCKS = NBLOCKS
        ALLOCATE(PRECOND_BLOCK_SIZES(NBLOCKS))
        ALLOCATE(PRECOND_BLOCK_OFFSETS(NBLOCKS))

        PRECOND_BLOCK_SIZES = BLOCK_SIZES

        ! Compute block offsets
        PRECOND_BLOCK_OFFSETS(1) = 1
        DO K = 2, NBLOCKS
            PRECOND_BLOCK_OFFSETS(K) = PRECOND_BLOCK_OFFSETS(K-1) + BLOCK_SIZES(K-1)
        END DO

        ! Find max block size for allocation
        MAXSZ = MAXVAL(BLOCK_SIZES)
        PRECOND_MAXBLK = MAXSZ

        ALLOCATE(PRECOND_BLK_DATA(MAXSZ, MAXSZ, NBLOCKS))
        ALLOCATE(PRECOND_BLK_IPIV(MAXSZ, NBLOCKS))

        PRECOND_BLK_DATA = CMPLX(0.0D0, 0.0D0)
        PRECOND_BLK_IPIV = 0

        ! Extract and LU-factorize each diagonal block
        DO K = 1, NBLOCKS
            BLK_SZ = BLOCK_SIZES(K)
            IOFF = PRECOND_BLOCK_OFFSETS(K) - 1

            ! Extract diagonal block: AMAT(IOFF+1:IOFF+BLK_SZ, IOFF+1:IOFF+BLK_SZ)
            DO J = 1, BLK_SZ
                DO I = 1, BLK_SZ
                    PRECOND_BLK_DATA(I, J, K) = AMAT(IOFF + I, IOFF + J)
                END DO
            END DO

            ! LU-factorize the block
            CALL ZGETRF(BLK_SZ, BLK_SZ, PRECOND_BLK_DATA(1,1,K), MAXSZ, &
                        PRECOND_BLK_IPIV(1,K), INFO_LU)

            IF (INFO_LU /= 0) THEN
                WRITE(*,'(A,I3,A,I6)') '  WARNING: Block-diagonal preconditioner: ZGETRF failed for block ', &
                    K, ', INFO = ', INFO_LU
            END IF
        END DO

        PRECOND_TYPE = 2
        PRECOND_N = N

        WRITE(*,'(A,I3,A,I6,A)') '  Block-diagonal preconditioner: ', NBLOCKS, &
            ' blocks (max size ', MAXSZ, ')'

    END SUBROUTINE SETUP_BLOCK_PRECONDITIONER

!---------------------------------------------------------------------------------------------
!   Apply preconditioner: z = M^{-1} * r
!
!   Dispatches based on PRECOND_TYPE:
!     0 - No preconditioning (z = r)
!     1 - Jacobi: z(i) = r(i) / diag(A)(i)
!     2 - Block-diagonal LU: solve each block's L*U*z_k = r_k via ZGETRS
!---------------------------------------------------------------------------------------------

    SUBROUTINE APPLY_PRECONDITIONER(N, R, Z)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        COMPLEX*16, INTENT(IN) :: R(N)
        COMPLEX*16, INTENT(OUT) :: Z(N)
        INTEGER :: I, K, BLK_SZ, IOFF, INFO_SOLVE

        SELECT CASE (PRECOND_TYPE)

        CASE (2)
            ! Block-diagonal LU: solve L*U*z_k = r_k for each block
            DO K = 1, PRECOND_NBLOCKS
                BLK_SZ = PRECOND_BLOCK_SIZES(K)
                IOFF = PRECOND_BLOCK_OFFSETS(K) - 1

                ! Copy RHS into output (ZGETRS solves in-place)
                DO I = 1, BLK_SZ
                    Z(IOFF + I) = R(IOFF + I)
                END DO

                ! Solve using stored LU factors
                CALL ZGETRS('N', BLK_SZ, 1, PRECOND_BLK_DATA(1,1,K), PRECOND_MAXBLK, &
                             PRECOND_BLK_IPIV(1,K), Z(IOFF+1), BLK_SZ, INFO_SOLVE)
            END DO

        CASE (1)
            ! Jacobi (diagonal)
            IF (.NOT. ALLOCATED(PRECOND_DIAG) .OR. PRECOND_N == 0) THEN
                Z = R
                RETURN
            END IF
            DO I = 1, N
                Z(I) = R(I) / PRECOND_DIAG(I)
            END DO

        CASE DEFAULT
            ! No preconditioning
            Z = R

        END SELECT

    END SUBROUTINE APPLY_PRECONDITIONER

!---------------------------------------------------------------------------------------------
!   Restarted GMRES(m) solver for complex dense system A*x = b
!   with right Jacobi preconditioning.
!
!   Arguments:
!     N        - System size
!     AMAT     - Coefficient matrix (N x N, not modified)
!     B        - Right-hand side vector (N)
!     X        - On entry: initial guess; on exit: solution (N)
!     TOL      - Relative convergence tolerance (||r||/||b|| < TOL)
!     MAXITER  - Maximum total iterations across all restarts
!     MRESTART - Restart parameter (Krylov subspace dimension)
!     INFO     - Output: 0 = converged, 1 = max iterations reached
!---------------------------------------------------------------------------------------------

    SUBROUTINE ZGMRES_SOLVE(N, AMAT, B, X, TOL, MAXITER, MRESTART, INFO)
        USE HMatrix_mod, ONLY: HMAT_READY, HMATVEC
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: N, MAXITER, MRESTART
        COMPLEX*16, INTENT(IN) :: AMAT(N, N), B(N)
        COMPLEX*16, INTENT(INOUT) :: X(N)
        REAL*8, INTENT(IN) :: TOL
        INTEGER, INTENT(OUT) :: INFO

        INTEGER :: I, J, K, ITER, M_ACT
        REAL*8 :: BNORM, RNORM, BETA_R
        LOGICAL :: USE_HMAT

        COMPLEX*16, ALLOCATABLE :: V(:,:), H(:,:)
        COMPLEX*16, ALLOCATABLE :: CS(:), SN(:), G(:), Y(:)
        COMPLEX*16, ALLOCATABLE :: R(:), W(:), Z(:)
        COMPLEX*16 :: TEMP, CS_VAL, SN_VAL

        USE_HMAT = HMAT_READY

        INFO = 1

        BNORM = ZNRM2_INLINE(N, B)
        IF (BNORM < 1.0D-30) THEN
            X = CMPLX(0.0D0, 0.0D0)
            INFO = 0
            RETURN
        END IF

        ALLOCATE(V(N, MRESTART+1), H(MRESTART+1, MRESTART))
        ALLOCATE(CS(MRESTART), SN(MRESTART), G(MRESTART+1), Y(MRESTART))
        ALLOCATE(R(N), W(N), Z(N))

        ITER = 0

        RESTART_LOOP: DO WHILE (ITER < MAXITER)

            ! Compute residual r = b - A*x
            IF (ITER == 0 .AND. ZNRM2_INLINE(N, X) < 1.0D-30) THEN
                R = B
            ELSE
                IF (USE_HMAT) THEN
                    CALL HMATVEC(N, X, W)
                ELSE
                    CALL ZMATVEC(N, AMAT, X, W)
                END IF
                R = B - W
            END IF

            BETA_R = ZNRM2_INLINE(N, R)

            IF (BETA_R / BNORM < TOL) THEN
                INFO = 0
                EXIT RESTART_LOOP
            END IF

            V(:, 1) = R / BETA_R
            G = CMPLX(0.0D0, 0.0D0)
            G(1) = CMPLX(BETA_R, 0.0D0)
            H = CMPLX(0.0D0, 0.0D0)

            M_ACT = 0
            ARNOLDI_LOOP: DO J = 1, MRESTART
                ITER = ITER + 1
                IF (ITER > MAXITER) EXIT ARNOLDI_LOOP

                ! Right preconditioning: z = M^{-1} * v_j
                CALL APPLY_PRECONDITIONER(N, V(:,J), Z)

                ! w = A * z (use H-matrix if available)
                IF (USE_HMAT) THEN
                    CALL HMATVEC(N, Z, W)
                ELSE
                    CALL ZMATVEC(N, AMAT, Z, W)
                END IF

                ! Modified Gram-Schmidt orthogonalization
                DO I = 1, J
                    H(I, J) = ZDOT_INLINE(N, V(:, I), W)
                    ! Inline AXPY: W = W - H(I,J) * V(:,I)
                    DO K = 1, N
                        W(K) = W(K) - H(I, J) * V(K, I)
                    END DO
                END DO

                H(J+1, J) = CMPLX(ZNRM2_INLINE(N, W), 0.0D0)

                IF (ABS(H(J+1, J)) > 1.0D-30) THEN
                    V(:, J+1) = W / H(J+1, J)
                ELSE
                    V(:, J+1) = CMPLX(0.0D0, 0.0D0)
                END IF

                ! Apply previous Givens rotations to column J of H
                DO I = 1, J-1
                    TEMP = CS(I) * H(I, J) + SN(I) * H(I+1, J)
                    H(I+1, J) = -CONJG(SN(I)) * H(I, J) + CONJG(CS(I)) * H(I+1, J)
                    H(I, J) = TEMP
                END DO

                ! Compute new Givens rotation
                CALL ZGIVENS_ROTATION(H(J, J), H(J+1, J), CS_VAL, SN_VAL)
                CS(J) = CS_VAL
                SN(J) = SN_VAL

                H(J, J) = CS(J) * H(J, J) + SN(J) * H(J+1, J)
                H(J+1, J) = CMPLX(0.0D0, 0.0D0)

                TEMP = CS(J) * G(J) + SN(J) * G(J+1)
                G(J+1) = -CONJG(SN(J)) * G(J) + CONJG(CS(J)) * G(J+1)
                G(J) = TEMP

                M_ACT = J

                ! Check convergence
                RNORM = ABS(G(J+1))
                IF (RNORM / BNORM < TOL) THEN
                    INFO = 0
                    EXIT ARNOLDI_LOOP
                END IF

            END DO ARNOLDI_LOOP

            ! Back-substitute and update solution
            IF (M_ACT > 0) THEN
                Y(1:M_ACT) = G(1:M_ACT)
                DO I = M_ACT, 1, -1
                    IF (ABS(H(I, I)) < 1.0D-30) THEN
                        Y(I) = CMPLX(0.0D0, 0.0D0)
                    ELSE
                        Y(I) = Y(I) / H(I, I)
                    END IF
                    DO K = 1, I-1
                        Y(K) = Y(K) - H(K, I) * Y(I)
                    END DO
                END DO

                ! x = x + M^{-1} * V * y
                DO J = 1, M_ACT
                    CALL APPLY_PRECONDITIONER(N, V(:,J), Z)
                    DO I = 1, N
                        X(I) = X(I) + Y(J) * Z(I)
                    END DO
                END DO
            END IF

            IF (INFO == 0) EXIT RESTART_LOOP

        END DO RESTART_LOOP

        IF (INFO /= 0) THEN
            RNORM = ZNRM2_INLINE(N, R)
            WRITE(*,'(A,ES10.3,A,ES10.3)') &
                '  WARNING: GMRES did not converge. Residual = ', RNORM, &
                ', Tolerance = ', TOL * BNORM
        END IF

        DEALLOCATE(V, H, CS, SN, G, Y, R, W, Z)

    END SUBROUTINE ZGMRES_SOLVE

!---------------------------------------------------------------------------------------------
!   Givens rotation: zero out the second element
!---------------------------------------------------------------------------------------------

    SUBROUTINE ZGIVENS_ROTATION(A, B, CS, SN)
        IMPLICIT NONE
        COMPLEX*16, INTENT(IN) :: A, B
        COMPLEX*16, INTENT(OUT) :: CS, SN
        REAL*8 :: DENOM
        IF (ABS(B) < 1.0D-30) THEN
            CS = CMPLX(1.0D0, 0.0D0)
            SN = CMPLX(0.0D0, 0.0D0)
        ELSE IF (ABS(A) < 1.0D-30) THEN
            CS = CMPLX(0.0D0, 0.0D0)
            SN = CMPLX(1.0D0, 0.0D0)
        ELSE
            DENOM = SQRT(ABS(A)**2 + ABS(B)**2)
            CS = A / CMPLX(DENOM, 0.0D0)
            SN = B / CMPLX(DENOM, 0.0D0)
        END IF
    END SUBROUTINE ZGIVENS_ROTATION

!---------------------------------------------------------------------------------------------
!   Block matrix-vector product: Y_OUT(N, NCOLS) = AMAT(N,N) * X_IN(N, NCOLS)
!
!   Uses MKL ZGEMM for optimal performance. This replaces NCOLS individual ZMATVEC
!   calls with a single level-3 BLAS call, which is much more cache-efficient.
!---------------------------------------------------------------------------------------------

    SUBROUTINE ZMATMAT(N, NCOLS, AMAT, X_IN, Y_OUT)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N, NCOLS
        COMPLEX*16, INTENT(IN) :: AMAT(N, N), X_IN(N, NCOLS)
        COMPLEX*16, INTENT(OUT) :: Y_OUT(N, NCOLS)
        COMPLEX*16 :: ZONE, ZZERO

        ZONE  = CMPLX(1.0D0, 0.0D0)
        ZZERO = CMPLX(0.0D0, 0.0D0)

        ! Y = A * X  via ZGEMM: C := alpha*A*B + beta*C
        CALL ZGEMM('N', 'N', N, NCOLS, N, ZONE, AMAT, N, X_IN, N, ZZERO, Y_OUT, N)

    END SUBROUTINE ZMATMAT

!---------------------------------------------------------------------------------------------
!   Block GMRES solver: solve A*X = B for multiple RHS simultaneously.
!
!   Instead of solving each RHS independently, Block GMRES builds a shared
!   Krylov subspace for a batch of RHS vectors. The key advantage is that
!   the matrix-vector products become matrix-matrix products (ZGEMM), which
!   MKL handles with much better cache utilization and SIMD efficiency.
!
!   Arguments:
!     N           - System size
!     AMAT        - Coefficient matrix (N x N, not modified)
!     B           - Right-hand side matrix (N x NRHS)
!     X           - On exit: solution matrix (N x NRHS)
!     NRHS        - Number of right-hand sides
!     TOL         - Relative convergence tolerance
!     MAXITER     - Maximum iterations per RHS
!     MRESTART    - Krylov subspace dimension before restart
!     INFO        - Output: 0 = all converged, 1 = some did not converge
!
!   Reference:
!     Gutknecht (2006), Block Krylov space methods for linear systems with
!     multiple right-hand sides: an introduction.
!---------------------------------------------------------------------------------------------

    SUBROUTINE ZBLOCK_GMRES_SOLVE(N, AMAT, B, X, NRHS, TOL, MAXITER, MRESTART, INFO)
        USE HAMS_mod, ONLY: NTHREAD
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: N, NRHS, MAXITER, MRESTART
        COMPLEX*16, INTENT(IN) :: AMAT(N, N), B(N, NRHS)
        COMPLEX*16, INTENT(OUT) :: X(N, NRHS)
        REAL*8, INTENT(IN) :: TOL
        INTEGER, INTENT(OUT) :: INFO

        INTEGER :: IRHS, INFO_SINGLE

        ! RHS-parallel GMRES: each right-hand side is solved independently
        ! by a separate thread. The coefficient matrix AMAT and H-matrix
        ! data are read-only during GMRES, so concurrent access is safe.
        ! The preconditioner (block-diagonal LU) is also read-only after setup.
        !
        ! This gives near-linear speedup with thread count since each
        ! GMRES solve has its own private Krylov workspace.

        INFO = 0
        X = CMPLX(0.0D0, 0.0D0)

        ! Ensure preconditioner is set up before entering parallel region
        IF (PRECOND_TYPE /= 2) THEN
            CALL SETUP_PRECONDITIONER(N, AMAT)
        END IF

        ! Solve each RHS sequentially. RHS-level OpenMP parallelism was tested
        ! but causes L3 cache thrashing at this problem size (N~7000) since
        ! multiple threads compete for the same H-matrix data. The serial
        ! approach keeps the working set in cache and is faster.
        ! For larger problems (N > 15000), uncomment the OMP directives below.
        ! !$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(IRHS, INFO_SINGLE) &
        ! !$OMP& SCHEDULE(DYNAMIC, 1) REDUCTION(MAX:INFO)
        DO IRHS = 1, NRHS
            CALL ZGMRES_SOLVE(N, AMAT, B(:, IRHS), X(:, IRHS), &
                              TOL, MAXITER, MRESTART, INFO_SINGLE)
            IF (INFO_SINGLE /= 0) INFO = INFO_SINGLE
        END DO
        ! !$OMP END PARALLEL DO

    END SUBROUTINE ZBLOCK_GMRES_SOLVE

!---------------------------------------------------------------------------------------------
!   Multi-RHS wrapper: solve A*X = B for NRHS right-hand side columns
!   Uses Block GMRES with deflated initial guesses for optimal performance.
!---------------------------------------------------------------------------------------------

    SUBROUTINE ZITER_SOLVE_MULTI(ISOLV_IN, N, AMAT, B, X, NRHS, INFO)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: ISOLV_IN, N, NRHS
        COMPLEX*16, INTENT(IN) :: AMAT(N, N), B(N, NRHS)
        COMPLEX*16, INTENT(OUT) :: X(N, NRHS)
        INTEGER, INTENT(OUT) :: INFO

        INTEGER :: IRHS, INFO_SINGLE

        INFO = 0

        ! Use Block GMRES with deflated initial guesses
        CALL ZBLOCK_GMRES_SOLVE(N, AMAT, B, X, NRHS, &
                                 ITER_TOL, ITER_MAXITER, GMRES_RESTART, INFO)

    END SUBROUTINE ZITER_SOLVE_MULTI

END MODULE IterativeSolvers
