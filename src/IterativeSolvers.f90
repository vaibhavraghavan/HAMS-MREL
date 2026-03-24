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
!    Self-contained iterative solvers for complex double-precision (COMPLEX*16)
!    dense linear systems A*x = b arising from the BEM formulation.
!
!    Implements:
!      - Restarted GMRES(m) with Jacobi (diagonal) preconditioner
!      - BiCGSTAB with Jacobi (diagonal) preconditioner
!      - Multi-RHS wrapper for solving multiple right-hand sides
!
!    Uses only standard BLAS routines (ZGEMV, ZDOTC, DZNRM2) already linked
!    through LAPACK. No additional dependencies.
!
!    References:
!      - Saad & Schultz (1986), GMRES: A generalized minimal residual algorithm
!      - van der Vorst (1992), Bi-CGSTAB: A fast and smoothly converging variant
!        of Bi-CG for the solution of nonsymmetric linear systems
!
!  ------------------------------------------------------------------------------------------------------

MODULE IterativeSolvers

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: ZGMRES_SOLVE
    PUBLIC :: ZBICGSTAB_SOLVE
    PUBLIC :: ZITER_SOLVE_MULTI

    ! Default solver parameters
    INTEGER, PUBLIC :: GMRES_RESTART = 50        ! Restart parameter for GMRES(m)
    REAL*8, PUBLIC  :: ITER_TOL = 1.0D-7         ! Relative convergence tolerance
    INTEGER, PUBLIC :: ITER_MAXITER = 500         ! Maximum iterations

CONTAINS

!---------------------------------------------------------------------------------------------
!   Restarted GMRES(m) solver for complex dense system A*x = b
!   with Jacobi (diagonal) preconditioner M = diag(A).
!
!   Arguments:
!     N       - System size
!     AMAT    - Coefficient matrix (N x N, COMPLEX*16, not modified)
!     B       - Right-hand side vector (N)
!     X       - On entry: initial guess; on exit: solution (N)
!     TOL     - Relative convergence tolerance (||r|| / ||b|| < TOL)
!     MAXITER - Maximum total iterations across all restarts
!     MRESTART - Restart parameter (Krylov subspace dimension before restart)
!     INFO    - Output: 0 = converged, 1 = max iterations reached
!---------------------------------------------------------------------------------------------

    SUBROUTINE ZGMRES_SOLVE(N, AMAT, B, X, TOL, MAXITER, MRESTART, INFO)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: N, MAXITER, MRESTART
        COMPLEX*16, INTENT(IN) :: AMAT(N, N), B(N)
        COMPLEX*16, INTENT(INOUT) :: X(N)
        REAL*8, INTENT(IN) :: TOL
        INTEGER, INTENT(OUT) :: INFO

        ! Local variables
        INTEGER :: I, J, K, ITER, M_ACT
        REAL*8 :: BNORM, RNORM, BETA_R

        ! Workspace
        COMPLEX*16, ALLOCATABLE :: V(:,:)       ! Arnoldi basis vectors (N, MRESTART+1)
        COMPLEX*16, ALLOCATABLE :: H(:,:)       ! Hessenberg matrix (MRESTART+1, MRESTART)
        COMPLEX*16, ALLOCATABLE :: CS(:), SN(:) ! Givens rotation parameters
        COMPLEX*16, ALLOCATABLE :: G(:)         ! Residual vector in Hessenberg space
        COMPLEX*16, ALLOCATABLE :: Y(:)         ! Solution in Krylov subspace
        COMPLEX*16, ALLOCATABLE :: R(:), W(:), Z(:)  ! Work vectors
        COMPLEX*16, ALLOCATABLE :: DIAG(:)      ! Jacobi preconditioner (diagonal of A)
        COMPLEX*16 :: TEMP, CS_VAL, SN_VAL
        REAL*8 :: DZNRM2
        COMPLEX*16 :: ZDOTC

        INFO = 1   ! Default: not converged

        ! Compute ||b||
        BNORM = DZNRM2(N, B, 1)
        IF (BNORM < 1.0D-30) THEN
            X = CMPLX(0.0D0, 0.0D0)
            INFO = 0
            RETURN
        END IF

        ! Allocate workspace
        ALLOCATE(V(N, MRESTART+1), H(MRESTART+1, MRESTART))
        ALLOCATE(CS(MRESTART), SN(MRESTART), G(MRESTART+1), Y(MRESTART))
        ALLOCATE(R(N), W(N), Z(N), DIAG(N))

        ! Extract diagonal for Jacobi preconditioner
        DO I = 1, N
            DIAG(I) = AMAT(I, I)
            IF (ABS(DIAG(I)) < 1.0D-30) DIAG(I) = CMPLX(1.0D0, 0.0D0)
        END DO

        ITER = 0

        ! Outer restart loop
        RESTART_LOOP: DO WHILE (ITER < MAXITER)

            ! Compute residual r = b - A*x
            R = B
            CALL ZGEMV('N', N, N, CMPLX(-1.0D0, 0.0D0), AMAT, N, X, 1, &
                        CMPLX(1.0D0, 0.0D0), R, 1)

            BETA_R = DZNRM2(N, R, 1)

            ! Check convergence
            IF (BETA_R / BNORM < TOL) THEN
                INFO = 0
                EXIT RESTART_LOOP
            END IF

            ! Initialize Arnoldi
            V(:, 1) = R / BETA_R
            G = CMPLX(0.0D0, 0.0D0)
            G(1) = CMPLX(BETA_R, 0.0D0)
            H = CMPLX(0.0D0, 0.0D0)

            ! Arnoldi iteration
            M_ACT = 0
            ARNOLDI_LOOP: DO J = 1, MRESTART
                ITER = ITER + 1
                IF (ITER > MAXITER) EXIT ARNOLDI_LOOP

                ! Apply Jacobi preconditioner: z = M^{-1} * v_j
                DO I = 1, N
                    Z(I) = V(I, J) / DIAG(I)
                END DO

                ! Compute w = A * z
                CALL ZGEMV('N', N, N, CMPLX(1.0D0, 0.0D0), AMAT, N, Z, 1, &
                            CMPLX(0.0D0, 0.0D0), W, 1)

                ! Modified Gram-Schmidt orthogonalization
                DO I = 1, J
                    H(I, J) = ZDOTC(N, V(:, I), 1, W, 1)
                    W = W - H(I, J) * V(:, I)
                END DO

                H(J+1, J) = CMPLX(DZNRM2(N, W, 1), 0.0D0)

                IF (ABS(H(J+1, J)) > 1.0D-30) THEN
                    V(:, J+1) = W / H(J+1, J)
                ELSE
                    ! Lucky breakdown
                    V(:, J+1) = CMPLX(0.0D0, 0.0D0)
                END IF

                ! Apply previous Givens rotations to column J of H
                DO I = 1, J-1
                    TEMP = CS(I) * H(I, J) + SN(I) * H(I+1, J)
                    H(I+1, J) = -CONJG(SN(I)) * H(I, J) + CONJG(CS(I)) * H(I+1, J)
                    H(I, J) = TEMP
                END DO

                ! Compute new Givens rotation to eliminate H(J+1, J)
                CALL ZGIVENS_ROTATION(H(J, J), H(J+1, J), CS_VAL, SN_VAL)
                CS(J) = CS_VAL
                SN(J) = SN_VAL

                ! Apply rotation to H and g
                H(J, J) = CS(J) * H(J, J) + SN(J) * H(J+1, J)
                H(J+1, J) = CMPLX(0.0D0, 0.0D0)

                TEMP = CS(J) * G(J) + SN(J) * G(J+1)
                G(J+1) = -CONJG(SN(J)) * G(J) + CONJG(CS(J)) * G(J+1)
                G(J) = TEMP

                M_ACT = J

                ! Check convergence from residual estimate
                RNORM = ABS(G(J+1))
                IF (RNORM / BNORM < TOL) THEN
                    INFO = 0
                    EXIT ARNOLDI_LOOP
                END IF

            END DO ARNOLDI_LOOP

            ! Solve upper triangular system H * y = g
            IF (M_ACT > 0) THEN
                Y(1:M_ACT) = G(1:M_ACT)
                DO I = M_ACT, 1, -1
                    Y(I) = Y(I) / H(I, I)
                    DO K = 1, I-1
                        Y(K) = Y(K) - H(K, I) * Y(I)
                    END DO
                END DO

                ! Update solution: x = x + V * M^{-1} * y
                ! First compute the preconditioned update
                DO J = 1, M_ACT
                    DO I = 1, N
                        X(I) = X(I) + Y(J) * V(I, J) / DIAG(I)
                    END DO
                END DO
            END IF

            IF (INFO == 0) EXIT RESTART_LOOP

        END DO RESTART_LOOP

        ! Print warning if not converged
        IF (INFO /= 0) THEN
            RNORM = DZNRM2(N, R, 1)
            WRITE(*,'(A,ES10.3,A,ES10.3)') &
                '  WARNING: GMRES did not converge. Residual = ', RNORM, &
                ', Tolerance = ', TOL * BNORM
        END IF

        DEALLOCATE(V, H, CS, SN, G, Y, R, W, Z, DIAG)

    END SUBROUTINE ZGMRES_SOLVE

!---------------------------------------------------------------------------------------------
!   Compute Givens rotation to zero out the second element
!   Given (a, b), compute (cs, sn) such that:
!     | cs    sn  | | a |   | r |
!     |-sn*   cs* | | b | = | 0 |
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
!   BiCGSTAB solver for complex dense system A*x = b
!   with Jacobi (diagonal) preconditioner M = diag(A).
!
!   Arguments:
!     N       - System size
!     AMAT    - Coefficient matrix (N x N, COMPLEX*16, not modified)
!     B       - Right-hand side vector (N)
!     X       - On entry: initial guess; on exit: solution (N)
!     TOL     - Relative convergence tolerance
!     MAXITER - Maximum iterations
!     INFO    - Output: 0 = converged, 1 = max iterations, 2 = breakdown
!---------------------------------------------------------------------------------------------

    SUBROUTINE ZBICGSTAB_SOLVE(N, AMAT, B, X, TOL, MAXITER, INFO)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: N, MAXITER
        COMPLEX*16, INTENT(IN) :: AMAT(N, N), B(N)
        COMPLEX*16, INTENT(INOUT) :: X(N)
        REAL*8, INTENT(IN) :: TOL
        INTEGER, INTENT(OUT) :: INFO

        ! Local variables
        INTEGER :: I, ITER
        REAL*8 :: BNORM, RNORM
        COMPLEX*16 :: RHO, RHO_OLD, ALPHA, OMEGA, BETA_C, TEMP

        ! Workspace
        COMPLEX*16, ALLOCATABLE :: R(:), R_HAT(:), P(:), V(:), S(:), T(:)
        COMPLEX*16, ALLOCATABLE :: P_HAT(:), S_HAT(:)
        COMPLEX*16, ALLOCATABLE :: DIAG(:)

        REAL*8 :: DZNRM2
        COMPLEX*16 :: ZDOTC

        INFO = 1

        BNORM = DZNRM2(N, B, 1)
        IF (BNORM < 1.0D-30) THEN
            X = CMPLX(0.0D0, 0.0D0)
            INFO = 0
            RETURN
        END IF

        ALLOCATE(R(N), R_HAT(N), P(N), V(N), S(N), T(N))
        ALLOCATE(P_HAT(N), S_HAT(N), DIAG(N))

        ! Extract diagonal for Jacobi preconditioner
        DO I = 1, N
            DIAG(I) = AMAT(I, I)
            IF (ABS(DIAG(I)) < 1.0D-30) DIAG(I) = CMPLX(1.0D0, 0.0D0)
        END DO

        ! Compute initial residual r = b - A*x
        R = B
        CALL ZGEMV('N', N, N, CMPLX(-1.0D0, 0.0D0), AMAT, N, X, 1, &
                    CMPLX(1.0D0, 0.0D0), R, 1)

        ! Choose shadow residual r_hat = r
        R_HAT = R

        ! Initialize
        RHO_OLD = CMPLX(1.0D0, 0.0D0)
        ALPHA   = CMPLX(1.0D0, 0.0D0)
        OMEGA   = CMPLX(1.0D0, 0.0D0)
        V = CMPLX(0.0D0, 0.0D0)
        P = CMPLX(0.0D0, 0.0D0)

        DO ITER = 1, MAXITER

            ! rho = <r_hat, r>
            RHO = ZDOTC(N, R_HAT, 1, R, 1)

            ! Check for breakdown
            IF (ABS(RHO) < 1.0D-30) THEN
                INFO = 2
                EXIT
            END IF

            ! beta = (rho / rho_old) * (alpha / omega)
            BETA_C = (RHO / RHO_OLD) * (ALPHA / OMEGA)

            ! p = r + beta * (p - omega * v)
            P = R + BETA_C * (P - OMEGA * V)

            ! Apply preconditioner: p_hat = M^{-1} * p
            DO I = 1, N
                P_HAT(I) = P(I) / DIAG(I)
            END DO

            ! v = A * p_hat
            CALL ZGEMV('N', N, N, CMPLX(1.0D0, 0.0D0), AMAT, N, P_HAT, 1, &
                        CMPLX(0.0D0, 0.0D0), V, 1)

            ! alpha = rho / <r_hat, v>
            TEMP = ZDOTC(N, R_HAT, 1, V, 1)
            IF (ABS(TEMP) < 1.0D-30) THEN
                INFO = 2
                EXIT
            END IF
            ALPHA = RHO / TEMP

            ! s = r - alpha * v
            S = R - ALPHA * V

            ! Check if s is small enough
            RNORM = DZNRM2(N, S, 1)
            IF (RNORM / BNORM < TOL) THEN
                X = X + ALPHA * P_HAT
                INFO = 0
                EXIT
            END IF

            ! Apply preconditioner: s_hat = M^{-1} * s
            DO I = 1, N
                S_HAT(I) = S(I) / DIAG(I)
            END DO

            ! t = A * s_hat
            CALL ZGEMV('N', N, N, CMPLX(1.0D0, 0.0D0), AMAT, N, S_HAT, 1, &
                        CMPLX(0.0D0, 0.0D0), T, 1)

            ! omega = <t, s> / <t, t>
            TEMP = ZDOTC(N, T, 1, T, 1)
            IF (ABS(TEMP) < 1.0D-30) THEN
                INFO = 2
                EXIT
            END IF
            OMEGA = ZDOTC(N, T, 1, S, 1) / TEMP

            ! Update solution: x = x + alpha * p_hat + omega * s_hat
            X = X + ALPHA * P_HAT + OMEGA * S_HAT

            ! Update residual: r = s - omega * t
            R = S - OMEGA * T

            ! Check convergence
            RNORM = DZNRM2(N, R, 1)
            IF (RNORM / BNORM < TOL) THEN
                INFO = 0
                EXIT
            END IF

            RHO_OLD = RHO

        END DO

        ! Print warning if not converged
        IF (INFO /= 0) THEN
            RNORM = DZNRM2(N, R, 1)
            IF (INFO == 1) THEN
                WRITE(*,'(A,ES10.3,A,ES10.3)') &
                    '  WARNING: BiCGSTAB did not converge. Residual = ', RNORM, &
                    ', Tolerance = ', TOL * BNORM
            ELSE IF (INFO == 2) THEN
                WRITE(*,'(A,ES10.3)') &
                    '  WARNING: BiCGSTAB breakdown. Residual = ', RNORM
            END IF
        END IF

        DEALLOCATE(R, R_HAT, P, V, S, T, P_HAT, S_HAT, DIAG)

    END SUBROUTINE ZBICGSTAB_SOLVE

!---------------------------------------------------------------------------------------------
!   Multi-RHS wrapper: solve A*X = B for multiple right-hand side columns
!
!   Arguments:
!     ISOLV_IN - Solver type (2 = GMRES, 3 = BiCGSTAB)
!     N        - System size
!     AMAT     - Coefficient matrix (N x N, not modified)
!     B        - Right-hand side matrix (N x NRHS, not modified)
!     X        - Output solution matrix (N x NRHS)
!     NRHS     - Number of right-hand side columns
!     INFO     - Output: 0 = all converged, nonzero = at least one failed
!---------------------------------------------------------------------------------------------

    SUBROUTINE ZITER_SOLVE_MULTI(ISOLV_IN, N, AMAT, B, X, NRHS, INFO)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: ISOLV_IN, N, NRHS
        COMPLEX*16, INTENT(IN) :: AMAT(N, N), B(N, NRHS)
        COMPLEX*16, INTENT(OUT) :: X(N, NRHS)
        INTEGER, INTENT(OUT) :: INFO

        INTEGER :: IRHS, INFO_SINGLE

        INFO = 0

        DO IRHS = 1, NRHS
            ! Zero initial guess
            X(:, IRHS) = CMPLX(0.0D0, 0.0D0)

            IF (ISOLV_IN == 2) THEN
                CALL ZGMRES_SOLVE(N, AMAT, B(:, IRHS), X(:, IRHS), &
                                  ITER_TOL, ITER_MAXITER, GMRES_RESTART, INFO_SINGLE)
            ELSE IF (ISOLV_IN == 3) THEN
                CALL ZBICGSTAB_SOLVE(N, AMAT, B(:, IRHS), X(:, IRHS), &
                                     ITER_TOL, ITER_MAXITER, INFO_SINGLE)
            ELSE
                WRITE(*,*) '  ERROR: Invalid ISOLV value in ZITER_SOLVE_MULTI:', ISOLV_IN
                INFO = -1
                RETURN
            END IF

            IF (INFO_SINGLE /= 0) INFO = INFO_SINGLE
        END DO

    END SUBROUTINE ZITER_SOLVE_MULTI

END MODULE IterativeSolvers
