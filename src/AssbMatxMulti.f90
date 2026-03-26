!  ------------------------------------------------------------------------------------------------------
!                                                               
!    Program HAMS-MREL for the diffraction and radiation of waves 
!    by 3D structures.
! 
!  License:
! 
!    This routine is part of HAMS-MREL.
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
!  ------------------------------------------------------------------------------------------------------
MODULE AssbMatxMulti

   USE HAMS_mod
   USE Body_mod
   USE Const_mod 
   USE WaveDyn_mod
   USE PanelMesh_mod

   USE BodyIntgrMulti
   USE PatcVelct
   USE IterativeSolvers
   USE HMatrix_mod, ONLY: HMATRIX_BUILD
   IMPLICIT NONE

   PRIVATE
   
  ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: ASSB_LEFT_MULTI
   PUBLIC :: ASSB_RBC_MULTI
   
   PUBLIC :: ASSB_DBC_MULTI
   PUBLIC :: RADIATION_SOLVER_MULTI
   PUBLIC :: DIFFRACTION_SOLVER_MULTI
   
CONTAINS
! ------------------------------------------------------------------- 
!    Calculate the element contribution, assembly the left-hand
!    side coefficient matrix [A]
! ------------------------------------------------------------------- 
      SUBROUTINE ASSB_LEFT_MULTI(AMAT,IPIV,NELEM,NSYS)
      IMPLICIT  NONE
!     
      INTEGER BODY_N,NELEM_GLOBAL
      INTEGER,INTENT(IN):: NELEM,NSYS
      INTEGER,INTENT(OUT):: IPIV(NELEM,NSYS)
      COMPLEX*16,INTENT(OUT):: AMAT(NELEM,NELEM,NSYS)
      
      INTEGER IEL,JEL,IS,IP,IRR,MD,FLAG,INFO
      REAL*8 DIST
      COMPLEX*16 TINDP(4)
      
      NELEM_GLOBAL=0
      
      ! Creating combination of the variables required for future computations
      DO BODY_N=1,NBODY
        IF (BODY_N.EQ.1) THEN
         DS_MULTI_COMB(1:NELEM_MULTI(BODY_N))=DS_MULTI(BODY_N,:)
         DXYZ_MULTI_COMB(1:NELEM_MULTI(BODY_N),:)=DXYZ_MULTI_P(BODY_N,:,:)
        ELSEIF (BODY_N.EQ.NBODY) THEN
         DS_MULTI_COMB(NELEM_GLOBAL+1:TNELEM)=DS_MULTI(BODY_N,:)
         DXYZ_MULTI_COMB(NELEM_GLOBAL+1:TNELEM,:)=DXYZ_MULTI_P(BODY_N,:,:)
        ELSE
         DS_MULTI_COMB(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N))=DS_MULTI(BODY_N,:)
         DXYZ_MULTI_COMB(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N),:)=DXYZ_MULTI_P(BODY_N,:,:)
        ENDIF
        NELEM_GLOBAL=NELEM_GLOBAL+NELEM_MULTI(BODY_N)
      ENDDO

      IRR=1
      AMAT=CMPLX(0.0D0,0.0D0)

!$OMP PARALLEL NUM_THREADS(NTHREAD)          
!$OMP DO PRIVATE(IEL,JEL,IP,IS,FLAG,DIST,TINDP) !$OMP REDUCTION(+:AMAT)                            ! The +:AMAT is basically addition operator which will add the value of the coefficiens for all elements
      
      DO  1000 IEL=1,  NELEM
        
       DO 100  IP=1,  NSYS

        AMAT(IEL,IEL,IP)=2.D0*PI
        
100    CONTINUE

       DO 200 JEL=1,  NELEM
            
        DIST=SQRT((XYZ_GLOBAL_MULTI_COMB_P(IEL,1)-XYZ_GLOBAL_MULTI_COMB_P(JEL,1))**2+(XYZ_GLOBAL_MULTI_COMB_P(IEL,2)-XYZ_GLOBAL_MULTI_COMB_P(JEL,2))**2+(XYZ_GLOBAL_MULTI_COMB_P(IEL,3)-XYZ_GLOBAL_MULTI_COMB_P(JEL,3))**2)              ! Obtaining the distance between the source and field points
        IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF
          
        TINDP=(0.0D0, 0.0D0)
        
        DO  200   IS=1,  NSYS
        
         CALL BODINT_MULTI_LEFT(IS,IEL,JEL,TINDP,FLAG)

         DO IP=1, NSYS
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           AMAT(IEL,JEL,IP)=AMAT(IEL,JEL,IP)+RXY_MULTI(IS,IP)*TINDP(IS)
          ELSE
           AMAT(IEL,JEL,IP)=AMAT(IEL,JEL,IP)+RXY_MULTI(IS,IP)*TINDP(IS)
          ENDIF
         ENDDO
        
200    CONTINUE
     
1000   CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF (ISOLV .EQ. 1) THEN
!$omp parallel do private(NTHREAD)
         DO IP=1, NSYS
            CALL ZGETRF( NELEM, NELEM, AMAT(:,:,IP), NELEM, IPIV(:,IP), INFO )
         ENDDO
!$omp end parallel do
      ENDIF

       RETURN
      END SUBROUTINE ASSB_LEFT_MULTI
      
      
! ------------------------------------------------------------------- 
!    Calculate the element contribution, assembly the right-hand
!    side vector[B] and solve the equations
! ------------------------------------------------------------------- 
      SUBROUTINE ASSB_RBC_MULTI(BRMAT,NELEM,NSYS,NELEM_START,NELEM_END)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM,NSYS,NELEM_START,NELEM_END
      COMPLEX*16,INTENT(INOUT):: BRMAT(NELEM,6,NSYS)
      
      INTEGER IEL,JEL,IS,IP,IRR,MD,FLAG
      REAL*8 DIST
      COMPLEX*16 TINRD(4,6,4),BTMP(6,4)

      IRR=1                                 ! This signifies no removal of irregular frequencies
      
!$OMP PARALLEL NUM_THREADS(NTHREAD)           
!$OMP DO PRIVATE(IEL,JEL,MD,IP,IS,FLAG,DIST,TINRD,BTMP) !$OMP REDUCTION(+:BRMAT)
      
      DO  1000 IEL=1,  TNELEM
            
       BTMP=CMPLX(0.0D0,0.0D0)

        DO 200 JEL=NELEM_START,  NELEM_END
            
         DIST=SQRT((XYZ_GLOBAL_MULTI_COMB_P(IEL,1)-XYZ_GLOBAL_MULTI_COMB_P(JEL,1))**2+(XYZ_GLOBAL_MULTI_COMB_P(IEL,2)-XYZ_GLOBAL_MULTI_COMB_P(JEL,2))**2+(XYZ_GLOBAL_MULTI_COMB_P(IEL,3)-XYZ_GLOBAL_MULTI_COMB_P(JEL,3))**2)
         IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
          FLAG=1
         ELSE
          FLAG=0
         ENDIF

         TINRD=CMPLX(0.0D0, 0.0D0)
        
         DO  200   IS=1,  NSYS
        
          CALL RBC_MULTI_RIGHT(IS,IEL,JEL,TINRD,FLAG)
        
          DO MD=1,  6
          DO IP=1, NSYS
            BTMP(MD,IP)=BTMP(MD,IP)+TINRD(IS,MD,IP)
          ENDDO
          ENDDO
        
200     CONTINUE
        
         DO  300  MD=1,  6
         DO  300  IP=1, NSYS

         DO  300  IS=1, NSYS
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           BRMAT(IEL,MD,IP)=BRMAT(IEL,MD,IP)+RXY_MULTI(IP,IS)*BTMP(MD,IS)
          ELSE
           BRMAT(IEL,MD,IP)=BRMAT(IEL,MD,IP)+RXY_MULTI(IP,IS)*BTMP(MD,IS)
          ENDIF
300      CONTINUE
     
1000   CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL

       RETURN
      END SUBROUTINE ASSB_RBC_MULTI
      
! ------------------------------------------------------------------- 
!    Calculate the element contribution, assembly the right-hand
!    side vector[B] and solve the equations - Diffraction equation
! ------------------------------------------------------------------- 
      SUBROUTINE ASSB_DBC_MULTI(BDMAT,NELEM,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM,NSYS
      COMPLEX*16,INTENT(OUT):: BDMAT(NELEM,NSYS)
      
      INTEGER IEL,JEL,IS,IP,IRR,MD,FLAG
      REAL*8 XP,YP,ZP,DIST
      COMPLEX*16 TINRD(4,4),BTMP(4)

      IRR=1
      MD=7
      BDMAT=CMPLX(0.0D0,0.0D0)

!$OMP PARALLEL NUM_THREADS(NTHREAD)           
!$OMP DO PRIVATE(XP,YP,ZP,IEL,JEL,IP,IS,FLAG,DIST,TINRD,BTMP) !$OMP REDUCTION(+:BDMAT)
      
      DO  1000 IEL=1,  NELEM
            
       BTMP=CMPLX(0.0D0,0.0D0)
       
       IF (ISOL.EQ.2) THEN
            
         DO 100  IP=1,  NSYS
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           XP=RX_MULTI(IP,1)*XYZ_GLOBAL_MULTI_COMB_P(IEL,1)
           YP=RX_MULTI(IP,2)*XYZ_GLOBAL_MULTI_COMB_P(IEL,2)
           ZP=         XYZ_GLOBAL_MULTI_COMB_P(IEL,3)
          ELSE
           XP=RY_MULTI(IP,1)*XYZ_GLOBAL_MULTI_COMB_P(IEL,1)
           YP=RY_MULTI(IP,2)*XYZ_GLOBAL_MULTI_COMB_P(IEL,2)
           ZP=         XYZ_GLOBAL_MULTI_COMB_P(IEL,3)
          ENDIF
          BTMP(IP)=4.D0*PI*VINP(XP,YP,ZP,XW(1),XW(2),BETA)
100      CONTINUE
         
       ELSEIF (ISOL.EQ.1) THEN       

        DO 200 JEL=1,  NELEM
            
         DIST=SQRT((XYZ_GLOBAL_MULTI_COMB_P(IEL,1)-XYZ_GLOBAL_MULTI_COMB_P(JEL,1))**2+(XYZ_GLOBAL_MULTI_COMB_P(IEL,2)-XYZ_GLOBAL_MULTI_COMB_P(JEL,2))**2+(XYZ_GLOBAL_MULTI_COMB_P(IEL,3)-XYZ_GLOBAL_MULTI_COMB_P(JEL,3))**2)
         IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
          FLAG=1
         ELSE
          FLAG=0
         ENDIF
         TINRD=CMPLX(0.0D0, 0.0D0)
        
         DO  200   IS=1,  NSYS
        
           CALL DBC_MULTI_RIGHT(IS,IEL,JEL,TINRD,FLAG)
        
           DO IP=1, NSYS
             BTMP(IP)=BTMP(IP)+TINRD(IS,IP)
           ENDDO
        
200     CONTINUE
        
       ELSE
           
        PRINT*,"  Error: The input for ISOL should be either 1 or 2."
        STOP
        
       ENDIF
       
        
         DO  300  IP=1, NSYS

         DO  300  IS=1, NSYS
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           BDMAT(IEL,IP)=BDMAT(IEL,IP)+RXY_MULTI(IP,IS)*BTMP(IS) 
          ELSE
           BDMAT(IEL,IP)=BDMAT(IEL,IP)+RXY_MULTI(IP,IS)*BTMP(IS)  
          ENDIF       
300      CONTINUE
     
1000   CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL

      RETURN
      END SUBROUTINE ASSB_DBC_MULTI
      
! -----------------------------------------------------------------------------
!     Solve the radiation potentials for not removing irregular frequencies
! -----------------------------------------------------------------------------

      SUBROUTINE RADIATION_SOLVER_MULTI(AMAT,BRMAT,IPIV,MXPOT,NELEM,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM,NSYS
      INTEGER,INTENT(IN):: IPIV(NELEM,NSYS)
      COMPLEX*16,INTENT(IN):: AMAT(NELEM,NELEM,NSYS),BRMAT(NELEM,6*NBODY,NSYS)
      COMPLEX*16,INTENT(OUT):: MXPOT(NELEM,6*NBODY+1,NSYS)
      
      INTEGER IEL,IS,IP,MD,INFO
      COMPLEX*16,ALLOCATABLE:: ATMAT(:,:,:),BRTMAT(:,:,:)
      COMPLEX*16,ALLOCATABLE:: AWORK(:,:),BWORK(:,:),XWORK(:,:)

      ALLOCATE(ATMAT(NELEM,NELEM,NSYS),BRTMAT(NELEM,6*NBODY,NSYS))

      IF (ISOLV .EQ. 1) THEN
          ATMAT=AMAT
          BRTMAT=BRMAT
!$omp parallel do private(NTHREAD)
          DO IP=1, NSYS
               CALL ZGETRS( 'No transpose', NELEM, 6*NBODY, ATMAT(:,:,IP), NELEM, IPIV(:,IP), BRTMAT(:,:,IP), NELEM, INFO )
          ENDDO
!$omp end parallel do
      ELSE
          ALLOCATE(AWORK(NELEM,NELEM), BWORK(NELEM,6*NBODY), XWORK(NELEM,6*NBODY))
          BRTMAT = CMPLX(0.0D0, 0.0D0)
          DO IP=1, NSYS
               AWORK = AMAT(:,:,IP)
               BWORK = BRMAT(:,:,IP)
               XWORK = CMPLX(0.0D0, 0.0D0)
               ! Set up block-diagonal LU preconditioner using per-body self-interaction blocks
               CALL SETUP_BLOCK_PRECONDITIONER(NELEM, AWORK, NBODY, NELEM_MULTI)
               ! Build H-matrix for fast matvec in GMRES
               CALL HMATRIX_BUILD(NELEM, AWORK, XYZ_GLOBAL_MULTI_COMB_P)
               CALL ZITER_SOLVE_MULTI(ISOLV, NELEM, AWORK, BWORK, XWORK, 6*NBODY, INFO)
               BRTMAT(:,:,IP) = XWORK
          ENDDO
          DEALLOCATE(AWORK, BWORK, XWORK)
      ENDIF

      DO 1400 MD=1, 6*NBODY
      DO 1400 IP=1, NSYS
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,IS) !$OMP REDUCTION(+:MXPOT)
      DO 1450 IEL=1, NELEM
         MXPOT(IEL,MD,IP)=CMPLX(0.0D0, 0.0D0)
      DO 1460 IS=1, NSYS
         IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
          MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+BRTMAT(IEL,MD,IS)*RXY_MULTI(IP,IS)
         ELSE
          MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+BRTMAT(IEL,MD,IS)*RXY_MULTI(IP,IS)
         ENDIF
1460  CONTINUE
         MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)/NSYS
1450  CONTINUE
!$OMP END DO NOWAIT
!$OMP END PARALLEL
1400  CONTINUE
      
      DEALLOCATE(ATMAT,BRTMAT)
      
      RETURN
      END SUBROUTINE RADIATION_SOLVER_MULTI
      
! -----------------------------------------------------------------------------
!     Solve the diffraction potentials for not removing irregular frequencies
! -----------------------------------------------------------------------------

      SUBROUTINE DIFFRACTION_SOLVER_MULTI(AMAT,BDMAT,IPIV,MXPOT,NELEM,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM,NSYS
      INTEGER,INTENT(IN):: IPIV(NELEM,NSYS)
      COMPLEX*16,INTENT(IN):: AMAT(NELEM,NELEM,NSYS),BDMAT(NELEM,NSYS)
      COMPLEX*16,INTENT(OUT):: MXPOT(NELEM,6*NBODY+1,NSYS)
      
      INTEGER IEL,IS,IP,MD,INFO
      COMPLEX*16,ALLOCATABLE:: ATMAT(:,:,:),BDTMAT(:,:)
      
      ALLOCATE(ATMAT(NELEM,NELEM,NSYS),BDTMAT(NELEM,NSYS))

      MD=6*NBODY+1

      IF (ISOLV .EQ. 1) THEN
          ATMAT=AMAT
          BDTMAT=BDMAT
!$omp parallel do private(NTHREAD)
          DO IP=1, NSYS
               CALL ZGETRS( 'No transpose', NELEM, 1, ATMAT(:,:,IP), NELEM, IPIV(:,IP), BDTMAT(:,IP), NELEM, INFO )
          ENDDO
!$omp end parallel do
      ELSE
          ! GMRES iterative solver with block-diagonal preconditioner
          DO IP=1, NSYS
               ATMAT(:,:,1) = AMAT(:,:,IP)
               BDTMAT(:,IP) = BDMAT(:,IP)
               ! Block preconditioner already set up during radiation solve
               ! If not yet set up (e.g., diffraction-only), set it up now
               IF (PRECOND_NBLOCKS == 0) THEN
                   CALL SETUP_BLOCK_PRECONDITIONER(NELEM, ATMAT(:,:,1), NBODY, NELEM_MULTI)
               END IF
               CALL ZGMRES_SOLVE(NELEM, ATMAT(:,:,1), BDMAT(:,IP), BDTMAT(:,IP), &
                                 ITER_TOL, ITER_MAXITER, GMRES_RESTART, INFO)
          ENDDO
      ENDIF

       DO 400 IP=1, NSYS
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,IS) !$OMP REDUCTION(+:MXPOT)
       DO 450 IEL=1, NELEM
         MXPOT(IEL,MD,IP)=CMPLX(0.0D0, 0.0D0)
       DO 460 IS=1, NSYS
         IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
          MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+BDTMAT(IEL,IS)*RXY_MULTI(IP,IS)
         ELSE
          MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+BDTMAT(IEL,IS)*RXY_MULTI(IP,IS)
         ENDIF
460    CONTINUE
         MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)/NSYS
450    CONTINUE
!$OMP END DO NOWAIT
!$OMP END PARALLEL
400    CONTINUE

      DEALLOCATE(ATMAT,BDTMAT)
      
      RETURN
      END SUBROUTINE DIFFRACTION_SOLVER_MULTI
!-------------------------------------------------------------------------------
END MODULE AssbMatxMulti
!*******************************************************************************