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
!    Hierarchical matrix (H-matrix) compression for BEM influence matrices.
!    Provides O(N k log N) matrix-vector products by exploiting the low-rank
!    structure of far-field panel interactions in the Green's function kernel.
!
!    The implementation uses:
!      - Geometric bisection for cluster tree construction
!      - Standard admissibility condition for block partitioning
!      - Adaptive Cross Approximation (ACA) for low-rank block compression
!      - OpenMP-parallelized H-matvec
!
!    References:
!      - Hackbusch, W. (1999). A sparse matrix arithmetic based on H-matrices.
!        Part I: Introduction to H-matrices. Computing, 62(2), 89-108.
!      - Bebendorf, M. (2000). Approximation of boundary element matrices.
!        Numer. Math., 86(4), 565-589.
!      - Rjasanow, S. & Steinbach, O. (2007). The Fast Solution of Boundary
!        Integral Equations. Springer.
!
!  ------------------------------------------------------------------------------------------------------

MODULE HMatrix_mod

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: HMATRIX_BUILD
    PUBLIC :: HMATVEC
    PUBLIC :: HMATRIX_DESTROY
    PUBLIC :: HMAT_READY

    ! ---------------------------------------------------------------
    ! Tuning parameters
    ! ---------------------------------------------------------------
    INTEGER, PARAMETER :: NLEAF_DEFAULT = 32     ! Max panels per leaf cluster
    REAL*8, PARAMETER  :: ETA_DEFAULT = 2.0D0    ! Admissibility parameter
    REAL*8, PARAMETER  :: ACA_TOL = 1.0D-4       ! ACA convergence tolerance
    INTEGER, PARAMETER :: RANK_MAX = 50          ! Maximum low-rank approximation rank
    INTEGER, PARAMETER :: MAX_CLUSTERS = 100000  ! Max cluster tree nodes
    INTEGER, PARAMETER :: MAX_BLOCKS = 200000    ! Max H-matrix blocks

    ! ---------------------------------------------------------------
    ! Cluster tree node
    ! ---------------------------------------------------------------
    TYPE :: ClusterNode
        INTEGER :: idx_start, idx_end    ! Range in permutation array
        INTEGER :: nsize                 ! Number of panels in cluster
        INTEGER :: left, right           ! Children (0 = leaf)
        REAL*8 :: center(3)              ! Geometric center
        REAL*8 :: diameter               ! Cluster diameter
    END TYPE ClusterNode

    ! ---------------------------------------------------------------
    ! H-matrix block (dense or low-rank)
    ! ---------------------------------------------------------------
    TYPE :: HBlock
        INTEGER :: row_clust, col_clust  ! Cluster indices
        INTEGER :: row_start, row_end    ! Row range in permuted ordering
        INTEGER :: col_start, col_end    ! Column range in permuted ordering
        INTEGER :: nrows, ncols          ! Block dimensions
        INTEGER :: btype                 ! 0 = dense, 1 = low-rank
        INTEGER :: rank                  ! Low-rank approximation rank (btype=1)
        COMPLEX*16, ALLOCATABLE :: D(:,:)   ! Dense block (nrows x ncols)
        COMPLEX*16, ALLOCATABLE :: U(:,:)   ! Low-rank factor U (nrows x rank)
        COMPLEX*16, ALLOCATABLE :: VH(:,:)  ! Low-rank factor V^H (rank x ncols)
    END TYPE HBlock

    ! ---------------------------------------------------------------
    ! Complete H-matrix structure
    ! ---------------------------------------------------------------
    TYPE :: HMatrixType
        INTEGER :: N                     ! Matrix dimension
        INTEGER :: nclusters             ! Number of cluster tree nodes
        INTEGER :: nblocks               ! Number of H-matrix blocks
        INTEGER :: nleaf                 ! Leaf cluster size
        REAL*8 :: eta                    ! Admissibility parameter
        TYPE(ClusterNode), ALLOCATABLE :: tree(:)   ! Cluster tree
        INTEGER, ALLOCATABLE :: perm(:)             ! Permutation: cluster order -> original
        INTEGER, ALLOCATABLE :: iperm(:)            ! Inverse: original -> cluster order
        TYPE(HBlock), ALLOCATABLE :: blocks(:)      ! Block list
        INTEGER :: n_dense               ! Count of dense blocks
        INTEGER :: n_lowrank             ! Count of low-rank blocks
        INTEGER :: total_rank            ! Sum of all low-rank block ranks
        LOGICAL :: ready                 ! True if H-matrix is built and usable
    END TYPE HMatrixType

    ! Module-level H-matrix instance
    TYPE(HMatrixType), SAVE :: HMAT
    LOGICAL, SAVE :: HMAT_READY = .FALSE.

    ! Module-level pointer to avoid copying the large matrix
    INTEGER, SAVE :: AMAT_N = 0
    COMPLEX*16, POINTER, SAVE :: AMAT_PTR(:,:) => NULL()

CONTAINS

!---------------------------------------------------------------------------------------------
!   Build H-matrix from the dense BEM influence matrix AMAT and panel coordinates.
!
!   Steps:
!     1. Build cluster tree from panel coordinates
!     2. Build block cluster tree (admissibility partitioning)
!     3. Compress admissible blocks via ACA
!     4. Store inadmissible blocks as dense
!
!   Arguments:
!     N     - Matrix dimension
!     AMAT  - Dense influence matrix (N x N), read-only
!     XYZ_P - Panel center coordinates (N x 3)
!---------------------------------------------------------------------------------------------

    SUBROUTINE HMATRIX_BUILD(N, AMAT, XYZ_P)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        COMPLEX*16, INTENT(IN), TARGET :: AMAT(N, N)
        REAL*8, INTENT(IN) :: XYZ_P(N, 3)

        REAL*8 :: T_START, T_END
        INTEGER :: MEM_DENSE, MEM_LR

        CALL CPU_TIME(T_START)

        ! Clean up any previous H-matrix
        CALL HMATRIX_DESTROY()

        ! Point to caller's matrix — NO copy (saves N*N*16 bytes)
        AMAT_N = N
        AMAT_PTR => AMAT

        HMAT%N = N
        HMAT%nleaf = NLEAF_DEFAULT
        HMAT%eta = ETA_DEFAULT
        HMAT%n_dense = 0
        HMAT%n_lowrank = 0
        HMAT%total_rank = 0

        ! Allocate cluster tree
        ALLOCATE(HMAT%tree(MAX_CLUSTERS))
        ALLOCATE(HMAT%perm(N))
        ALLOCATE(HMAT%iperm(N))
        ALLOCATE(HMAT%blocks(MAX_BLOCKS))

        ! Step 1: Build cluster tree (iterative, no recursion)
        HMAT%nclusters = 0
        CALL BUILD_CLUSTER_TREE(N, XYZ_P)

        ! Step 2+3: Build block cluster tree and compress (iterative, uses AMAT_PTR)
        HMAT%nblocks = 0
        CALL BUILD_BLOCK_TREE_ITERATIVE(1, 1)

        ! Trim allocations to actual size
        HMAT%ready = .TRUE.
        HMAT_READY = .TRUE.

        CALL CPU_TIME(T_END)

        ! Report compression statistics
        MEM_DENSE = 0
        MEM_LR = 0
        CALL COMPUTE_MEMORY(MEM_DENSE, MEM_LR)

        WRITE(*,'(A)') '  H-matrix compression:'
        WRITE(*,'(A,I6,A,I6,A,I6)') '    Blocks: ', HMAT%nblocks, &
            ' (dense: ', HMAT%n_dense, ', low-rank: ', HMAT%n_lowrank, ')'
        IF (HMAT%n_lowrank > 0) THEN
            WRITE(*,'(A,F5.1)') '    Average rank: ', &
                DBLE(HMAT%total_rank) / DBLE(MAX(HMAT%n_lowrank, 1))
        END IF
        WRITE(*,'(A,F5.1,A)') '    Compression ratio: ', &
            DBLE(N) * DBLE(N) * 16.0D0 / DBLE(MAX(MEM_DENSE + MEM_LR, 1)), 'x'
        WRITE(*,'(A,F6.2,A)') '    Build time: ', T_END - T_START, ' s'

    END SUBROUTINE HMATRIX_BUILD

!---------------------------------------------------------------------------------------------
!   Build cluster tree via recursive geometric bisection.
!   Splits panel set along longest bounding box dimension.
!---------------------------------------------------------------------------------------------

    SUBROUTINE BUILD_CLUSTER_TREE(N, XYZ_P)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN) :: XYZ_P(N, 3)

        INTEGER :: I
        REAL*8, ALLOCATABLE :: XYZ_WORK(:,:)

        ! Initialize permutation as identity
        DO I = 1, N
            HMAT%perm(I) = I
            HMAT%iperm(I) = I
        END DO

        ! Working copy of coordinates (will be permuted)
        ALLOCATE(XYZ_WORK(N, 3))
        XYZ_WORK = XYZ_P

        ! Build tree iteratively (no recursion, uses explicit stack)
        CALL SPLIT_CLUSTER_ITERATIVE(N, XYZ_WORK)

        DEALLOCATE(XYZ_WORK)

    END SUBROUTINE BUILD_CLUSTER_TREE

!---------------------------------------------------------------------------------------------
!   Recursive cluster splitting. Creates a new cluster node for indices
!   [IDX_START, IDX_END] and splits into two children if the cluster is
!   larger than NLEAF.
!
!   Returns the cluster node index in PARENT_IDX (0 on first call).
!---------------------------------------------------------------------------------------------

    SUBROUTINE SPLIT_CLUSTER_ITERATIVE(N, XYZ_WORK)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(INOUT) :: XYZ_WORK(N, 3)

        ! Explicit stack to avoid recursion (max depth ~20 for millions of panels)
        INTEGER, PARAMETER :: MAX_DEPTH = 64
        INTEGER :: STACK_IS(MAX_DEPTH), STACK_IE(MAX_DEPTH), STACK_PARENT(MAX_DEPTH)
        INTEGER :: STACK_SIDE(MAX_DEPTH)  ! 1=left, 2=right
        INTEGER :: SP  ! stack pointer

        INTEGER :: NODE_IDX, NSIZE, I, J, SPLIT_DIM
        REAL*8 :: BBOX_MIN(3), BBOX_MAX(3), BBOX_LEN(3)
        REAL*8 :: CX, CY, CZ, DIAM, SPLIT_VAL, MAXLEN
        REAL*8 :: TMPXYZ(3)
        INTEGER :: TMPPERM, LEFT_END, IDX_START, IDX_END, PARENT_NODE

        ! Push root onto stack
        SP = 1
        STACK_IS(1) = 1
        STACK_IE(1) = N
        STACK_PARENT(1) = 0
        STACK_SIDE(1) = 0

        DO WHILE (SP > 0)
            ! Pop from stack
            IDX_START = STACK_IS(SP)
            IDX_END = STACK_IE(SP)
            PARENT_NODE = STACK_PARENT(SP)
            I = STACK_SIDE(SP)  ! reuse I temporarily for side
            SP = SP - 1

            NSIZE = IDX_END - IDX_START + 1
            IF (NSIZE <= 0) CYCLE

            ! Create new cluster node
            HMAT%nclusters = HMAT%nclusters + 1
            NODE_IDX = HMAT%nclusters

            ! Link to parent
            IF (PARENT_NODE > 0) THEN
                IF (I == 1) THEN
                    HMAT%tree(PARENT_NODE)%left = NODE_IDX
                ELSE
                    HMAT%tree(PARENT_NODE)%right = NODE_IDX
                END IF
            END IF

            HMAT%tree(NODE_IDX)%idx_start = IDX_START
            HMAT%tree(NODE_IDX)%idx_end = IDX_END
            HMAT%tree(NODE_IDX)%nsize = NSIZE
            HMAT%tree(NODE_IDX)%left = 0
            HMAT%tree(NODE_IDX)%right = 0

            ! Compute bounding box and center
            BBOX_MIN(1) = XYZ_WORK(IDX_START, 1)
            BBOX_MIN(2) = XYZ_WORK(IDX_START, 2)
            BBOX_MIN(3) = XYZ_WORK(IDX_START, 3)
            BBOX_MAX = BBOX_MIN
            CX = 0.0D0; CY = 0.0D0; CZ = 0.0D0

            DO J = IDX_START, IDX_END
                IF (XYZ_WORK(J,1) < BBOX_MIN(1)) BBOX_MIN(1) = XYZ_WORK(J,1)
                IF (XYZ_WORK(J,2) < BBOX_MIN(2)) BBOX_MIN(2) = XYZ_WORK(J,2)
                IF (XYZ_WORK(J,3) < BBOX_MIN(3)) BBOX_MIN(3) = XYZ_WORK(J,3)
                IF (XYZ_WORK(J,1) > BBOX_MAX(1)) BBOX_MAX(1) = XYZ_WORK(J,1)
                IF (XYZ_WORK(J,2) > BBOX_MAX(2)) BBOX_MAX(2) = XYZ_WORK(J,2)
                IF (XYZ_WORK(J,3) > BBOX_MAX(3)) BBOX_MAX(3) = XYZ_WORK(J,3)
                CX = CX + XYZ_WORK(J,1)
                CY = CY + XYZ_WORK(J,2)
                CZ = CZ + XYZ_WORK(J,3)
            END DO

            HMAT%tree(NODE_IDX)%center(1) = CX / DBLE(NSIZE)
            HMAT%tree(NODE_IDX)%center(2) = CY / DBLE(NSIZE)
            HMAT%tree(NODE_IDX)%center(3) = CZ / DBLE(NSIZE)

            BBOX_LEN = BBOX_MAX - BBOX_MIN
            DIAM = SQRT(BBOX_LEN(1)**2 + BBOX_LEN(2)**2 + BBOX_LEN(3)**2)
            HMAT%tree(NODE_IDX)%diameter = DIAM

            ! Leaf node: no splitting
            IF (NSIZE <= HMAT%nleaf) CYCLE

            ! Find longest dimension for splitting
            SPLIT_DIM = 1
            MAXLEN = BBOX_LEN(1)
            IF (BBOX_LEN(2) > MAXLEN) THEN; SPLIT_DIM = 2; MAXLEN = BBOX_LEN(2); END IF
            IF (BBOX_LEN(3) > MAXLEN) THEN; SPLIT_DIM = 3; MAXLEN = BBOX_LEN(3); END IF

            SPLIT_VAL = 0.5D0 * (BBOX_MIN(SPLIT_DIM) + BBOX_MAX(SPLIT_DIM))

            ! Partition
            I = IDX_START
            J = IDX_END

            DO WHILE (I <= J)
                IF (XYZ_WORK(I, SPLIT_DIM) <= SPLIT_VAL) THEN
                    I = I + 1
                ELSE
                    TMPXYZ(1) = XYZ_WORK(I, 1)
                    TMPXYZ(2) = XYZ_WORK(I, 2)
                    TMPXYZ(3) = XYZ_WORK(I, 3)
                    XYZ_WORK(I, :) = XYZ_WORK(J, :)
                    XYZ_WORK(J, 1) = TMPXYZ(1)
                    XYZ_WORK(J, 2) = TMPXYZ(2)
                    XYZ_WORK(J, 3) = TMPXYZ(3)

                    TMPPERM = HMAT%perm(I)
                    HMAT%perm(I) = HMAT%perm(J)
                    HMAT%perm(J) = TMPPERM

                    J = J - 1
                END IF
            END DO

            LEFT_END = I - 1
            IF (LEFT_END < IDX_START) LEFT_END = IDX_START
            IF (LEFT_END >= IDX_END) LEFT_END = IDX_END - 1

            ! Push right child first (so left is processed first = depth-first)
            IF (SP + 2 > MAX_DEPTH) THEN
                WRITE(*,*) 'ERROR: Cluster tree stack overflow'
                STOP
            END IF

            SP = SP + 1
            STACK_IS(SP) = LEFT_END + 1
            STACK_IE(SP) = IDX_END
            STACK_PARENT(SP) = NODE_IDX
            STACK_SIDE(SP) = 2  ! right child

            SP = SP + 1
            STACK_IS(SP) = IDX_START
            STACK_IE(SP) = LEFT_END
            STACK_PARENT(SP) = NODE_IDX
            STACK_SIDE(SP) = 1  ! left child
        END DO

        ! Build inverse permutation once at the end
        DO I = 1, N
            HMAT%iperm(HMAT%perm(I)) = I
        END DO

    END SUBROUTINE SPLIT_CLUSTER_ITERATIVE

!---------------------------------------------------------------------------------------------
!   Build block cluster tree recursively.
!   For each cluster pair (σ, τ), test admissibility:
!     - Admissible: compress via ACA → low-rank block
!     - Inadmissible leaf×leaf: store as dense block
!     - Inadmissible non-leaf: recurse into children
!---------------------------------------------------------------------------------------------

    SUBROUTINE BUILD_BLOCK_TREE_ITERATIVE(ROOT_ROW, ROOT_COL)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ROOT_ROW, ROOT_COL

        ! Explicit stack for iterative traversal
        ! 4-way branching with max depth ~20 => max stack size ~4*20 = 80
        ! Use generous size to be safe
        INTEGER, PARAMETER :: STACK_MAX = 200000
        INTEGER, ALLOCATABLE :: STACK_R(:), STACK_C(:)
        INTEGER :: SP

        INTEGER :: CLUST_ROW, CLUST_COL
        REAL*8 :: DIST, DIAM_MAX
        INTEGER :: RS, RE, CS, CE, NR, NC
        LOGICAL :: IS_LEAF_ROW, IS_LEAF_COL, IS_ADMISSIBLE

        ALLOCATE(STACK_R(STACK_MAX), STACK_C(STACK_MAX))

        ! Push root pair
        SP = 1
        STACK_R(1) = ROOT_ROW
        STACK_C(1) = ROOT_COL

        DO WHILE (SP > 0)
            ! Pop
            CLUST_ROW = STACK_R(SP)
            CLUST_COL = STACK_C(SP)
            SP = SP - 1

            RS = HMAT%tree(CLUST_ROW)%idx_start
            RE = HMAT%tree(CLUST_ROW)%idx_end
            CS = HMAT%tree(CLUST_COL)%idx_start
            CE = HMAT%tree(CLUST_COL)%idx_end
            NR = RE - RS + 1
            NC = CE - CS + 1

            IF (NR <= 0 .OR. NC <= 0) CYCLE

            IS_LEAF_ROW = (HMAT%tree(CLUST_ROW)%left == 0)
            IS_LEAF_COL = (HMAT%tree(CLUST_COL)%left == 0)

            DIST = SQRT( &
                (HMAT%tree(CLUST_ROW)%center(1) - HMAT%tree(CLUST_COL)%center(1))**2 + &
                (HMAT%tree(CLUST_ROW)%center(2) - HMAT%tree(CLUST_COL)%center(2))**2 + &
                (HMAT%tree(CLUST_ROW)%center(3) - HMAT%tree(CLUST_COL)%center(3))**2)

            DIAM_MAX = MAX(HMAT%tree(CLUST_ROW)%diameter, HMAT%tree(CLUST_COL)%diameter)

            IS_ADMISSIBLE = (DIST > HMAT%eta * DIAM_MAX) .AND. (NR > 1) .AND. (NC > 1)

            IF (IS_ADMISSIBLE) THEN
                CALL ADD_LOWRANK_BLOCK(CLUST_ROW, CLUST_COL)
                CYCLE
            END IF

            IF (IS_LEAF_ROW .AND. IS_LEAF_COL) THEN
                CALL ADD_DENSE_BLOCK(CLUST_ROW, CLUST_COL)
                CYCLE
            END IF

            ! Push children pairs onto stack
            IF (IS_LEAF_ROW) THEN
                ! 2 children
                IF (SP + 2 > STACK_MAX) THEN
                    WRITE(*,*) 'ERROR: Block tree stack overflow'
                    STOP
                END IF
                SP = SP + 1
                STACK_R(SP) = CLUST_ROW
                STACK_C(SP) = HMAT%tree(CLUST_COL)%right
                SP = SP + 1
                STACK_R(SP) = CLUST_ROW
                STACK_C(SP) = HMAT%tree(CLUST_COL)%left
            ELSE IF (IS_LEAF_COL) THEN
                IF (SP + 2 > STACK_MAX) THEN
                    WRITE(*,*) 'ERROR: Block tree stack overflow'
                    STOP
                END IF
                SP = SP + 1
                STACK_R(SP) = HMAT%tree(CLUST_ROW)%right
                STACK_C(SP) = CLUST_COL
                SP = SP + 1
                STACK_R(SP) = HMAT%tree(CLUST_ROW)%left
                STACK_C(SP) = CLUST_COL
            ELSE
                ! 4 children
                IF (SP + 4 > STACK_MAX) THEN
                    WRITE(*,*) 'ERROR: Block tree stack overflow'
                    STOP
                END IF
                SP = SP + 1
                STACK_R(SP) = HMAT%tree(CLUST_ROW)%right
                STACK_C(SP) = HMAT%tree(CLUST_COL)%right
                SP = SP + 1
                STACK_R(SP) = HMAT%tree(CLUST_ROW)%right
                STACK_C(SP) = HMAT%tree(CLUST_COL)%left
                SP = SP + 1
                STACK_R(SP) = HMAT%tree(CLUST_ROW)%left
                STACK_C(SP) = HMAT%tree(CLUST_COL)%right
                SP = SP + 1
                STACK_R(SP) = HMAT%tree(CLUST_ROW)%left
                STACK_C(SP) = HMAT%tree(CLUST_COL)%left
            END IF
        END DO

        DEALLOCATE(STACK_R, STACK_C)

    END SUBROUTINE BUILD_BLOCK_TREE_ITERATIVE

!---------------------------------------------------------------------------------------------
!   Add a dense block to the H-matrix.
!   Extracts the sub-matrix AMAT(row_indices, col_indices) using the permutation.
!---------------------------------------------------------------------------------------------

    SUBROUTINE ADD_DENSE_BLOCK(CLUST_ROW, CLUST_COL)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: CLUST_ROW, CLUST_COL

        INTEGER :: BID, RS, RE, CS, CE, NR, NC, I, J, ORIG_I, ORIG_J

        RS = HMAT%tree(CLUST_ROW)%idx_start
        RE = HMAT%tree(CLUST_ROW)%idx_end
        CS = HMAT%tree(CLUST_COL)%idx_start
        CE = HMAT%tree(CLUST_COL)%idx_end
        NR = RE - RS + 1
        NC = CE - CS + 1

        HMAT%nblocks = HMAT%nblocks + 1
        BID = HMAT%nblocks

        HMAT%blocks(BID)%row_clust = CLUST_ROW
        HMAT%blocks(BID)%col_clust = CLUST_COL
        HMAT%blocks(BID)%row_start = RS
        HMAT%blocks(BID)%row_end = RE
        HMAT%blocks(BID)%col_start = CS
        HMAT%blocks(BID)%col_end = CE
        HMAT%blocks(BID)%nrows = NR
        HMAT%blocks(BID)%ncols = NC
        HMAT%blocks(BID)%btype = 0  ! dense
        HMAT%blocks(BID)%rank = 0

        ALLOCATE(HMAT%blocks(BID)%D(NR, NC))

        ! Extract sub-matrix using permutation
        DO J = 1, NC
            ORIG_J = HMAT%perm(CS + J - 1)
            DO I = 1, NR
                ORIG_I = HMAT%perm(RS + I - 1)
                HMAT%blocks(BID)%D(I, J) = AMAT_PTR(ORIG_I, ORIG_J)
            END DO
        END DO

        HMAT%n_dense = HMAT%n_dense + 1

    END SUBROUTINE ADD_DENSE_BLOCK

!---------------------------------------------------------------------------------------------
!   Add a low-rank block via Adaptive Cross Approximation (ACA).
!
!   Partially-pivoted ACA builds A ≈ U * V^H incrementally:
!     1. Select pivot row (row with largest residual norm)
!     2. Compute residual row, find pivot column
!     3. Compute residual column
!     4. Update U and V^H
!     5. Check convergence: ||u_k|| * ||v_k|| < tol * ||U*V^H||_F
!
!   Falls back to dense storage if ACA rank exceeds RANK_MAX or
!   the block is too small for meaningful compression.
!---------------------------------------------------------------------------------------------

    SUBROUTINE ADD_LOWRANK_BLOCK(CLUST_ROW, CLUST_COL)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: CLUST_ROW, CLUST_COL

        INTEGER :: BID, RS, RE, CS, CE, NR, NC
        INTEGER :: I, J, K, PIVOT_ROW, PIVOT_COL, ORIG_I, ORIG_J
        REAL*8 :: MAX_VAL, CURR_VAL, FROB_SQ, UV_NORM
        COMPLEX*16 :: PIVOT_ENTRY
        COMPLEX*16, ALLOCATABLE :: U_WORK(:,:), VH_WORK(:,:)
        COMPLEX*16, ALLOCATABLE :: ROW_VEC(:), COL_VEC(:)
        REAL*8, ALLOCATABLE :: ROW_USED(:)
        INTEGER :: ACTUAL_RANK, KMAX

        RS = HMAT%tree(CLUST_ROW)%idx_start
        RE = HMAT%tree(CLUST_ROW)%idx_end
        CS = HMAT%tree(CLUST_COL)%idx_start
        CE = HMAT%tree(CLUST_COL)%idx_end
        NR = RE - RS + 1
        NC = CE - CS + 1

        KMAX = MIN(RANK_MAX, MIN(NR, NC))

        ! If block is too small, store as dense
        IF (NR <= 4 .OR. NC <= 4) THEN
            CALL ADD_DENSE_BLOCK(CLUST_ROW, CLUST_COL)
            RETURN
        END IF

        ALLOCATE(U_WORK(NR, KMAX), VH_WORK(KMAX, NC))
        ALLOCATE(ROW_VEC(NC), COL_VEC(NR))
        ALLOCATE(ROW_USED(NR))

        U_WORK = CMPLX(0.0D0, 0.0D0)
        VH_WORK = CMPLX(0.0D0, 0.0D0)
        ROW_USED = 0.0D0

        FROB_SQ = 0.0D0
        ACTUAL_RANK = 0

        ! ACA iteration
        DO K = 1, KMAX

            ! Select pivot row: first unused row (simple strategy)
            ! A more sophisticated approach would track residual norms
            PIVOT_ROW = 0
            IF (K == 1) THEN
                PIVOT_ROW = 1
            ELSE
                ! Find row with largest residual contribution from last step
                MAX_VAL = -1.0D0
                DO I = 1, NR
                    IF (ROW_USED(I) > 0.5D0) CYCLE
                    CURR_VAL = ABS(U_WORK(I, K-1))
                    IF (CURR_VAL > MAX_VAL) THEN
                        MAX_VAL = CURR_VAL
                        PIVOT_ROW = I
                    END IF
                END DO
                IF (PIVOT_ROW == 0) THEN
                    ! All rows used, try first unused
                    DO I = 1, NR
                        IF (ROW_USED(I) < 0.5D0) THEN
                            PIVOT_ROW = I
                            EXIT
                        END IF
                    END DO
                END IF
                IF (PIVOT_ROW == 0) EXIT  ! All rows exhausted
            END IF

            ROW_USED(PIVOT_ROW) = 1.0D0

            ! Compute residual row: R(pivot_row, :) = A(pivot_row, :) - sum_{l<k} U(pivot_row,l)*VH(l,:)
            ORIG_I = HMAT%perm(RS + PIVOT_ROW - 1)
            DO J = 1, NC
                ORIG_J = HMAT%perm(CS + J - 1)
                ROW_VEC(J) = AMAT_PTR(ORIG_I, ORIG_J)
            END DO

            ! Subtract previous approximation
            DO I = 1, K-1
                DO J = 1, NC
                    ROW_VEC(J) = ROW_VEC(J) - U_WORK(PIVOT_ROW, I) * VH_WORK(I, J)
                END DO
            END DO

            ! Find pivot column (largest entry in residual row)
            PIVOT_COL = 1
            MAX_VAL = ABS(ROW_VEC(1))
            DO J = 2, NC
                CURR_VAL = ABS(ROW_VEC(J))
                IF (CURR_VAL > MAX_VAL) THEN
                    MAX_VAL = CURR_VAL
                    PIVOT_COL = J
                END IF
            END DO

            PIVOT_ENTRY = ROW_VEC(PIVOT_COL)
            IF (ABS(PIVOT_ENTRY) < 1.0D-30) EXIT  ! Exact zero pivot, block is already approximated

            ! Compute residual column: R(:, pivot_col) = A(:, pivot_col) - sum_{l<k} U(:,l)*VH(l,pivot_col)
            ORIG_J = HMAT%perm(CS + PIVOT_COL - 1)
            DO I = 1, NR
                ORIG_I = HMAT%perm(RS + I - 1)
                COL_VEC(I) = AMAT_PTR(ORIG_I, ORIG_J)
            END DO

            DO I = 1, K-1
                DO J = 1, NR
                    COL_VEC(J) = COL_VEC(J) - U_WORK(J, I) * VH_WORK(I, PIVOT_COL)
                END DO
            END DO

            ! Set U(:,k) = col_vec / pivot_entry
            ! Set VH(k,:) = row_vec
            DO I = 1, NR
                U_WORK(I, K) = COL_VEC(I) / PIVOT_ENTRY
            END DO
            VH_WORK(K, :) = ROW_VEC(:)

            ACTUAL_RANK = K

            ! Check convergence: ||u_k|| * ||v_k|| < tol * ||approx||_F
            UV_NORM = 0.0D0
            DO I = 1, NR
                UV_NORM = UV_NORM + ABS(U_WORK(I, K))**2
            END DO
            UV_NORM = SQRT(UV_NORM)

            CURR_VAL = 0.0D0
            DO J = 1, NC
                CURR_VAL = CURR_VAL + ABS(VH_WORK(K, J))**2
            END DO
            CURR_VAL = SQRT(CURR_VAL)

            UV_NORM = UV_NORM * CURR_VAL

            ! Update Frobenius norm estimate (incremental)
            FROB_SQ = FROB_SQ + UV_NORM**2
            ! Add cross terms 2*Re(sum_{l<k} <u_l,u_k> * <v_k,v_l>)
            DO I = 1, K-1
                CURR_VAL = 0.0D0
                DO J = 1, NR
                    CURR_VAL = CURR_VAL + DBLE(CONJG(U_WORK(J,I)) * U_WORK(J,K)) * &
                                          DBLE(CONJG(VH_WORK(K,J)) * VH_WORK(I,J))
                END DO
                FROB_SQ = FROB_SQ + 2.0D0 * CURR_VAL
            END DO

            ! Convergence check
            IF (UV_NORM < ACA_TOL * SQRT(MAX(FROB_SQ, 1.0D-30))) EXIT

        END DO

        ! Check if compression is worthwhile: rank * (NR + NC) < NR * NC
        IF (ACTUAL_RANK * (NR + NC) >= NR * NC .OR. ACTUAL_RANK == 0) THEN
            ! Not worth compressing, store as dense
            DEALLOCATE(U_WORK, VH_WORK, ROW_VEC, COL_VEC, ROW_USED)
            CALL ADD_DENSE_BLOCK(CLUST_ROW, CLUST_COL)
            RETURN
        END IF

        ! Store as low-rank block
        HMAT%nblocks = HMAT%nblocks + 1
        BID = HMAT%nblocks

        HMAT%blocks(BID)%row_clust = CLUST_ROW
        HMAT%blocks(BID)%col_clust = CLUST_COL
        HMAT%blocks(BID)%row_start = RS
        HMAT%blocks(BID)%row_end = RE
        HMAT%blocks(BID)%col_start = CS
        HMAT%blocks(BID)%col_end = CE
        HMAT%blocks(BID)%nrows = NR
        HMAT%blocks(BID)%ncols = NC
        HMAT%blocks(BID)%btype = 1  ! low-rank
        HMAT%blocks(BID)%rank = ACTUAL_RANK

        ALLOCATE(HMAT%blocks(BID)%U(NR, ACTUAL_RANK))
        ALLOCATE(HMAT%blocks(BID)%VH(ACTUAL_RANK, NC))

        HMAT%blocks(BID)%U(:, 1:ACTUAL_RANK) = U_WORK(:, 1:ACTUAL_RANK)
        HMAT%blocks(BID)%VH(1:ACTUAL_RANK, :) = VH_WORK(1:ACTUAL_RANK, :)

        HMAT%n_lowrank = HMAT%n_lowrank + 1
        HMAT%total_rank = HMAT%total_rank + ACTUAL_RANK

        DEALLOCATE(U_WORK, VH_WORK, ROW_VEC, COL_VEC, ROW_USED)

    END SUBROUTINE ADD_LOWRANK_BLOCK

!---------------------------------------------------------------------------------------------
!   H-matrix - vector product: y = H * x
!
!   Permutes x to cluster ordering, accumulates block contributions,
!   then permutes y back to original ordering.
!
!   For dense blocks: y_σ += D * x_τ
!   For low-rank blocks: y_σ += U * (V^H * x_τ)  (two small matvecs)
!---------------------------------------------------------------------------------------------

    SUBROUTINE HMATVEC(N, X_IN, Y_OUT)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        COMPLEX*16, INTENT(IN) :: X_IN(N)
        COMPLEX*16, INTENT(OUT) :: Y_OUT(N)

        COMPLEX*16, ALLOCATABLE :: X_PERM(:), Y_PERM(:)
        COMPLEX*16 :: ZONE, ZZERO
        INTEGER :: I, BID, RS, NR, CS, NC, K
        COMPLEX*16 :: VTMP(RANK_MAX)

        ZONE  = CMPLX(1.0D0, 0.0D0)
        ZZERO = CMPLX(0.0D0, 0.0D0)

        ALLOCATE(X_PERM(N), Y_PERM(N))

        ! Permute input to cluster ordering
        DO I = 1, N
            X_PERM(I) = X_IN(HMAT%perm(I))
        END DO

        Y_PERM = CMPLX(0.0D0, 0.0D0)

        ! Serial block loop — the individual ZGEMV calls within each block
        ! are already fast (small dense blocks, tiny low-rank factors).
        ! OpenMP overhead for 12,000+ tiny blocks exceeds the parallel gain
        ! at this problem size. The outer GMRES and assembly loops provide
        ! sufficient parallelism.
        DO BID = 1, HMAT%nblocks
            RS = HMAT%blocks(BID)%row_start
            CS = HMAT%blocks(BID)%col_start
            NR = HMAT%blocks(BID)%nrows
            NC = HMAT%blocks(BID)%ncols

            IF (HMAT%blocks(BID)%btype == 0) THEN
                ! Dense block: y(RS:RS+NR-1) += D * x(CS:CS+NC-1)
                CALL ZGEMV('N', NR, NC, ZONE, HMAT%blocks(BID)%D, NR, &
                            X_PERM(CS), 1, ZONE, Y_PERM(RS), 1)

            ELSE
                ! Low-rank block: y += U * (V^H * x_τ)
                K = HMAT%blocks(BID)%rank

                ! Step 1: vtmp(1:K) = V^H * x_τ  (stack-allocated, no heap alloc)
                CALL ZGEMV('N', K, NC, ZONE, HMAT%blocks(BID)%VH, K, &
                            X_PERM(CS), 1, ZZERO, VTMP, 1)

                ! Step 2: y(RS:RS+NR-1) += U * vtmp
                CALL ZGEMV('N', NR, K, ZONE, HMAT%blocks(BID)%U, NR, &
                            VTMP, 1, ZONE, Y_PERM(RS), 1)
            END IF
        END DO

        ! Permute output back to original ordering
        DO I = 1, N
            Y_OUT(HMAT%perm(I)) = Y_PERM(I)
        END DO

        DEALLOCATE(X_PERM, Y_PERM)

    END SUBROUTINE HMATVEC

!---------------------------------------------------------------------------------------------
!   Compute memory usage of the H-matrix (in bytes)
!---------------------------------------------------------------------------------------------

    SUBROUTINE COMPUTE_MEMORY(MEM_DENSE, MEM_LR)
        IMPLICIT NONE
        INTEGER, INTENT(OUT) :: MEM_DENSE, MEM_LR
        INTEGER :: BID, NR, NC, K

        MEM_DENSE = 0
        MEM_LR = 0

        DO BID = 1, HMAT%nblocks
            NR = HMAT%blocks(BID)%nrows
            NC = HMAT%blocks(BID)%ncols

            IF (HMAT%blocks(BID)%btype == 0) THEN
                MEM_DENSE = MEM_DENSE + NR * NC * 16  ! COMPLEX*16 = 16 bytes
            ELSE
                K = HMAT%blocks(BID)%rank
                MEM_LR = MEM_LR + (NR * K + K * NC) * 16
            END IF
        END DO

    END SUBROUTINE COMPUTE_MEMORY

!---------------------------------------------------------------------------------------------
!   Destroy H-matrix and free all memory
!---------------------------------------------------------------------------------------------

    SUBROUTINE HMATRIX_DESTROY()
        IMPLICIT NONE
        INTEGER :: BID

        IF (.NOT. HMAT_READY) RETURN

        ! Free block data
        IF (ALLOCATED(HMAT%blocks)) THEN
            DO BID = 1, HMAT%nblocks
                IF (ALLOCATED(HMAT%blocks(BID)%D)) DEALLOCATE(HMAT%blocks(BID)%D)
                IF (ALLOCATED(HMAT%blocks(BID)%U)) DEALLOCATE(HMAT%blocks(BID)%U)
                IF (ALLOCATED(HMAT%blocks(BID)%VH)) DEALLOCATE(HMAT%blocks(BID)%VH)
            END DO
            DEALLOCATE(HMAT%blocks)
        END IF

        IF (ALLOCATED(HMAT%tree)) DEALLOCATE(HMAT%tree)
        IF (ALLOCATED(HMAT%perm)) DEALLOCATE(HMAT%perm)
        IF (ALLOCATED(HMAT%iperm)) DEALLOCATE(HMAT%iperm)

        HMAT%ready = .FALSE.
        HMAT_READY = .FALSE.
        HMAT%nblocks = 0
        HMAT%nclusters = 0

        ! Nullify module-level matrix pointer (caller owns the memory)
        NULLIFY(AMAT_PTR)
        AMAT_N = 0

    END SUBROUTINE HMATRIX_DESTROY

END MODULE HMatrix_mod
