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
MODULE AssbMatx_irrMulti

   USE HAMS_mod
   USE Body_mod
   USE Const_mod
   USE WaveDyn_mod
   USE PanelMesh_mod
   USE Inerfs_mod
   
   USE BodyIntgr_irrMulti
   USE PatcVelct
   USE IterativeSolvers, ONLY: ZITER_SOLVE_MULTI, SETUP_BLOCK_PRECONDITIONER
   IMPLICIT NONE
   
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: ASSB_LEFT_IRR_MULTI
   PUBLIC :: ASSB_RBC_IRR_MULTI
   
   PUBLIC :: ASSB_DBC_IRR_MULTI
   PUBLIC :: RADIATION_SOLVER_IRR_MULTI
   PUBLIC :: DIFFRACTION_SOLVER_IRR_MULTI
   
CONTAINS
! ------------------------------------------------------------------- 
!    Calculate the element contribution, assembly the left-hand
!    side coefficient matrix [A]
!    This subroutine is for suppression of irregular frequencies
! ------------------------------------------------------------------- 

      SUBROUTINE ASSB_LEFT_IRR_MULTI(AMAT,CMAT,IPIV,NELEM_TOTAL,iNELEM_TOTAL,TNELEM,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM_TOTAL,iNELEM_TOTAL,TNELEM,NSYS
      INTEGER,INTENT(OUT):: IPIV(NELEM_TOTAL,NSYS)
      COMPLEX*16,INTENT(OUT):: AMAT(TNELEM,TNELEM,NSYS),CMAT(NELEM_TOTAL,NELEM_TOTAL,NSYS)

      INTEGER IEL,JEL,KEL,IS,IP,IRR,FLAG,INFO
      REAL*8 DIST
      COMPLEX*16 TINDP(4),SUM

      ! DS_MULTI_COMB / DXYZ_MULTI_COMB / iDS_MULTI_COMB / iDXYZ_MULTI_COMB are populated
      ! once at startup by BuildCombinedArrays (Phase 2.3) — geometry-only, frequency-independent.

      AMAT=CMPLX(0.0D0,0.0D0)

      IRR=1

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,JEL,IP,IS,FLAG,DIST,TINDP) SCHEDULE(DYNAMIC,16)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.

      DO  1000 IEL=1,  NELEM_TOTAL
        
       DO 100  IP=1,  NSYS

        AMAT(IEL,IEL,IP)=2.D0*PI
        
100    CONTINUE
        
       DO 200 JEL=1,  NELEM_TOTAL

        DIST=SQRT((XYZ_GLOBAL_MULTI_COMB_P(IEL,1)-XYZ_GLOBAL_MULTI_COMB_P(JEL,1))**2+(XYZ_GLOBAL_MULTI_COMB_P(IEL,2)-XYZ_GLOBAL_MULTI_COMB_P(JEL,2))**2+(XYZ_GLOBAL_MULTI_COMB_P(IEL,3)-XYZ_GLOBAL_MULTI_COMB_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF

        TINDP=CMPLX(0.0D0, 0.0D0)

        DO  200   IS=1,  NSYS
        
         CALL BODINT_LEFT_IRR_MULTI(IS,IEL,JEL,TINDP,IRR,FLAG)

         DO IP=1, NSYS
          ! Both ISX/ISY branches were identical — collapsed.
          AMAT(IEL,JEL,IP)=AMAT(IEL,JEL,IP)+RXY_MULTI(IS,IP)*TINDP(IS)
         ENDDO
        
200    CONTINUE
     
1000   CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL
       
      IRR=3

!$OMP PARALLEL NUM_THREADS(NTHREAD)  
!$OMP DO PRIVATE(IEL,JEL,IP,IS,FLAG,DIST,TINDP) SCHEDULE(DYNAMIC,16)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.
      
      DO  3000 IEL=NELEM_TOTAL+1,  NELEM_TOTAL+iNELEM_TOTAL
       
       DO 500  IP=1,  NSYS

        AMAT(IEL,IEL,IP)=0.D0*PI
        
500    CONTINUE

       DO 600 JEL=1,  NELEM_TOTAL
            
        DIST=SQRT((iXYZ_GLOBAL_MULTI_COMB_P(IEL-NELEM_TOTAL,1)-XYZ_GLOBAL_MULTI_COMB_P(JEL,1))**2+(iXYZ_GLOBAL_MULTI_COMB_P(IEL-NELEM_TOTAL,2)-XYZ_GLOBAL_MULTI_COMB_P(JEL,2))**2+(iXYZ_GLOBAL_MULTI_COMB_P(IEL-NELEM_TOTAL,3)-XYZ_GLOBAL_MULTI_COMB_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF
        
        TINDP=CMPLX(0.0D0, 0.0D0)
        
        DO  600   IS=1,  NSYS            
        
         CALL BODINT_LEFT_IRR_MULTI(IS,IEL-NELEM_TOTAL,JEL,TINDP,IRR,FLAG)

         DO IP=1, NSYS
          ! Both ISX/ISY branches were identical — collapsed.
          AMAT(IEL,JEL,IP)=AMAT(IEL,JEL,IP)+RXY_MULTI(IS,IP)*TINDP(IS)
         ENDDO
        
600    CONTINUE

3000  CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL
!
!--------------------------------------------------------------
!    Matrix set-up for the least-square problem

      CMAT=CMPLX(0.0D0,0.0D0)

      DO  6000  IP=1, NSYS
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,JEL,KEL,SUM)
! SUM is per-iteration accumulator (PRIVATE); CMAT row-disjoint on IEL — no REDUCTION.
      DO  5000  IEL=1,  NELEM_TOTAL
      DO  5000  JEL=1,  NELEM_TOTAL
        
       SUM=CMPLX(0.0D0,0.0D0)
       DO  KEL=1,  TNELEM
          SUM=SUM+AMAT(KEL,IEL,IP)*AMAT(KEL,JEL,IP)
       ENDDO
       CMAT(IEL,JEL,IP)=SUM
        
5000   CONTINUE
!$OMP END DO NOWAIT
!$OMP END PARALLEL
6000   CONTINUE

! NSYS is 1 or 2 — let MKL's internal threading do the work in ZGETRF.
! Only factorize when direct-LU will be used; ISOLV>1 uses the unfactored normal-equations matrix.
       IF (ISOLV .EQ. 1) THEN
          DO IP=1, NSYS
             CALL ZGETRF( NELEM_TOTAL, NELEM_TOTAL, CMAT(:,:,IP), NELEM_TOTAL, IPIV(:,IP), INFO )
          ENDDO
       ENDIF

      RETURN
      END SUBROUTINE ASSB_LEFT_IRR_MULTI
      
! ------------------------------------------------------------------- 
!    Calculate the element contribution, assembly the right-hand
!    side vector[B]
!    This subroutine is for suppression of irregular frequencies
! ------------------------------------------------------------------- 

      SUBROUTINE ASSB_RBC_IRR_MULTI(BRMAT,DRMAT,AMAT,NELEM_TOTAL,iNELEM_TOTAL,TNELEM,NSYS,NELEM_START,NELEM_END)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM_TOTAL,iNELEM_TOTAL,TNELEM,NSYS,NELEM_START,NELEM_END
      COMPLEX*16,INTENT(IN):: AMAT(TNELEM,TNELEM,NSYS)
      COMPLEX*16,INTENT(OUT):: BRMAT(TNELEM,6,NSYS),DRMAT(NELEM_TOTAL,6,NSYS)
      
      INTEGER IEL,JEL,KEL,IS,IP,IRR,MD,FLAG
      REAL*8 DIST
      COMPLEX*16 TINRD(4,6,4),BTMP(6,4)
!
      IRR=1
  
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,JEL,MD,IP,IS,FLAG,DIST,TINRD,BTMP) SCHEDULE(DYNAMIC,16)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.
      
      DO  1000 IEL=1,  NELEM_TOTAL
            
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
        
          CALL RBC_IRR_MULTI(IS,IEL,JEL,TINRD,IRR,FLAG)
        
         DO MD=1,  6
         DO IP=1, NSYS
          BTMP(MD,IP)=BTMP(MD,IP)+TINRD(IS,MD,IP)
         ENDDO
         ENDDO
        
200     CONTINUE
       
        DO  300  MD=1,  6
        DO  300  IP=1, NSYS

        DO  300  IS=1, NSYS
          ! Both ISX/ISY branches were identical — collapsed.
          BRMAT(IEL,MD,IP)=BRMAT(IEL,MD,IP)+RXY_MULTI(IP,IS)*BTMP(MD,IS)
300     CONTINUE
     
1000    CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL
       
      IRR=3
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,JEL,MD,IP,IS,FLAG,DIST,TINRD,BTMP) SCHEDULE(DYNAMIC,16)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.
      
      DO  3000 IEL=NELEM_TOTAL+1,  NELEM_TOTAL+iNELEM_TOTAL
            
       BTMP=CMPLX(0.0D0,0.0D0)

        DO 600 JEL=NELEM_START,  NELEM_END

           DIST=SQRT((iXYZ_GLOBAL_MULTI_COMB_P(IEL-NELEM_TOTAL,1)-XYZ_GLOBAL_MULTI_COMB_P(JEL,1))**2 +&
                (iXYZ_GLOBAL_MULTI_COMB_P(IEL-NELEM_TOTAL,2)-XYZ_GLOBAL_MULTI_COMB_P(JEL,2))**2 +&
                (iXYZ_GLOBAL_MULTI_COMB_P(IEL-NELEM_TOTAL,3)-XYZ_GLOBAL_MULTI_COMB_P(JEL,3))**2)
         IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
          FLAG=1
         ELSE
          FLAG=0
         ENDIF

         TINRD=CMPLX(0.0D0, 0.0D0)
        
         DO  600   IS=1,  NSYS
        
          CALL RBC_IRR_MULTI(IS,IEL-NELEM_TOTAL,JEL,TINRD,IRR,FLAG)
        
         DO MD=1,  6
         DO IP=1, NSYS
          BTMP(MD,IP)=BTMP(MD,IP)+TINRD(IS,MD,IP)
         ENDDO
         ENDDO
        
600     CONTINUE

       DO  700  MD=1,  6
       DO  700  IP=1, NSYS

       DO  700  IS=1, NSYS
          ! Both ISX/ISY branches were identical — collapsed.
          BRMAT(IEL,MD,IP)=BRMAT(IEL,MD,IP)+RXY_MULTI(IP,IS)*BTMP(MD,IS)
700    CONTINUE
     
3000  CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL
!
!--------------------------------------------------------------
!    Matrix set-up for the least-square problem
!
      !WRITE(6, *) ' Right-hand side matrix setting-up for the least-square problem...'

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IP,IEL,KEL,MD)
! Disjoint writes on (IP,IEL) — no REDUCTION needed.

      DO  5000  IP=1, NSYS
      DO  5000  IEL=1,  NELEM_TOTAL
        
       DO  MD=1,6
       DO  KEL=1,  TNELEM
          DRMAT(IEL,MD,IP)=DRMAT(IEL,MD,IP)+BRMAT(KEL,MD,IP)*AMAT(KEL,IEL,IP)
       ENDDO
       ENDDO

5000  CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL

      RETURN
      END SUBROUTINE ASSB_RBC_IRR_MULTI

! ------------------------------------------------------------------- 
!    Calculate the element contribution, assembly the right-hand
!    side vector[B]
!    This subroutine is for suppression of irregular frequencies
! ------------------------------------------------------------------- 

      SUBROUTINE ASSB_DBC_IRR_MULTI(BDMAT,DDMAT,AMAT,NELEM_TOTAL,iNELEM_TOTAL,TNELEM,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM_TOTAL,iNELEM_TOTAL,TNELEM,NSYS
      COMPLEX*16,INTENT(IN):: AMAT(TNELEM,TNELEM,NSYS)
      COMPLEX*16,INTENT(OUT):: BDMAT(TNELEM,NSYS),DDMAT(NELEM_TOTAL,NSYS)
      
      INTEGER IEL,JEL,KEL,IS,IP,IRR,MD,FLAG
      REAL*8  XP,YP,ZP,DIST
      COMPLEX*16 TINRD(4,4),BTMP(4)

      MD=7
      BDMAT=CMPLX(0.0D0,0.0D0)
!
      IRR=1

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(XP,YP,ZP,IEL,JEL,IP,IS,FLAG,DIST,TINRD,BTMP) SCHEDULE(DYNAMIC,16)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.
      DO  1000 IEL=1,  NELEM_TOTAL

       BTMP=CMPLX(0.0D0,0.0D0)
       
       IF (ISOL.NE.1.AND.ISOL.NE.2) THEN
        PRINT*,"  Error: The input for ISOL should be either 1 or 2."
        STOP
       ENDIF

       ! Both ISOL=1 (total) and ISOL=2 (scattered) solve the same validated
       ! direct BIE for the scattered potential; the total/scattered choice is
       ! applied later by adding back the incident potential F0 (ISOL=1 only).
       DO 200 JEL=1,  NELEM_TOTAL

          DIST=SQRT((XYZ_GLOBAL_MULTI_COMB_P(IEL,1)-XYZ_GLOBAL_MULTI_COMB_P(JEL,1))**2+(XYZ_GLOBAL_MULTI_COMB_P(IEL,2)-XYZ_GLOBAL_MULTI_COMB_P(JEL,2))**2+(XYZ_GLOBAL_MULTI_COMB_P(IEL,3)-XYZ_GLOBAL_MULTI_COMB_P(JEL,3))**2)
          IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
           FLAG=1
          ELSE
           FLAG=0
          ENDIF

          TINRD=CMPLX(0.0D0, 0.0D0)

          DO  200   IS=1,  NSYS

            CALL DBC_IRR_MULTI(IS,IEL,JEL,TINRD,IRR,FLAG)

          DO IP=1, NSYS
            BTMP(IP)=BTMP(IP)+TINRD(IS,IP)
          ENDDO

200     CONTINUE
 
        DO  300  IP=1, NSYS

        DO  300  IS=1, NSYS
          ! Both ISX/ISY branches were identical — collapsed.
          BDMAT(IEL,IP)=BDMAT(IEL,IP)+RXY_MULTI(IP,IS)*BTMP(IS)
300     CONTINUE
     
1000   CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL
       
      IRR=3

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(XP,YP,ZP,IEL,JEL,IP,IS,FLAG,DIST,TINRD,BTMP) SCHEDULE(DYNAMIC,16)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.
      
      DO  3000 IEL=NELEM_TOTAL+1,  NELEM_TOTAL+iNELEM_TOTAL

       BTMP=CMPLX(0.0D0,0.0D0)
       
       IF (ISOL.NE.1.AND.ISOL.NE.2) THEN
        PRINT*,"  Error: The input for ISOL should be either 1 or 2."
        STOP
       ENDIF

       ! Both ISOL=1 (total) and ISOL=2 (scattered) solve the same validated
       ! direct BIE for the scattered potential; the total/scattered choice is
       ! applied later by adding back the incident potential F0 (ISOL=1 only).
       DO 600 JEL=1,  NELEM_TOTAL

        DIST=SQRT((iXYZ_GLOBAL_MULTI_COMB_P(IEL-NELEM_TOTAL,1)-XYZ_GLOBAL_MULTI_COMB_P(JEL,1))**2+(iXYZ_GLOBAL_MULTI_COMB_P(IEL-NELEM_TOTAL,2)-XYZ_GLOBAL_MULTI_COMB_P(JEL,2))**2+(iXYZ_GLOBAL_MULTI_COMB_P(IEL-NELEM_TOTAL,3)-XYZ_GLOBAL_MULTI_COMB_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF

          TINRD=CMPLX(0.0D0, 0.0D0)

          DO  600   IS=1,  NSYS

            CALL DBC_IRR_MULTI(IS,IEL-NELEM_TOTAL,JEL,TINRD,IRR,FLAG)

          DO IP=1, NSYS
            BTMP(IP)=BTMP(IP)+TINRD(IS,IP)
          ENDDO

600     CONTINUE
       
       DO  700  IP=1, NSYS

       DO  700  IS=1, NSYS
          ! Both ISX/ISY branches were identical — collapsed.
          BDMAT(IEL,IP)=BDMAT(IEL,IP)+RXY_MULTI(IP,IS)*BTMP(IS)
700    CONTINUE

3000   CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL
!
!--------------------------------------------------------------
!    Matrix set-up for the least-square problem

      DDMAT=CMPLX(0.0D0,0.0D0)

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IP,IEL,KEL)
! Disjoint writes on (IP,IEL) — no REDUCTION needed.

      DO  5000  IP=1, NSYS
      DO  5000  IEL=1,  NELEM_TOTAL
 
       DO  KEL=1,  TNELEM
          DDMAT(IEL,IP)=DDMAT(IEL,IP)+BDMAT(KEL,IP)*AMAT(KEL,IEL,IP)
       ENDDO

5000  CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL

      RETURN
      END SUBROUTINE ASSB_DBC_IRR_MULTI
      
! -----------------------------------------------------------------------------
!     Solve the radiation potentials after removing irregular frequencies
! -----------------------------------------------------------------------------

      SUBROUTINE RADIATION_SOLVER_IRR_MULTI(CMAT,DRMAT,IPIV,MXPOT,NELEM_TOTAL,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM_TOTAL,NSYS
      INTEGER,INTENT(IN):: IPIV(NELEM_TOTAL,NSYS)
      COMPLEX*16,INTENT(IN):: CMAT(NELEM_TOTAL,NELEM_TOTAL,NSYS),DRMAT(NELEM_TOTAL,6*NBODY,NSYS)
      COMPLEX*16,INTENT(OUT):: MXPOT(NELEM_TOTAL,6*NBODY+1,NSYS)
      
      INTEGER IEL,IS,IP,MD,INFO
      COMPLEX*16,ALLOCATABLE:: DRTMAT(:,:,:)
      COMPLEX*16,ALLOCATABLE:: CWORK(:,:)   ! GMRES preconditioner setup writes into its A argument.
      ! ZGETRS does not modify its A argument — pass CMAT (the LU factorization) directly,
      ! eliminating the previous full-size CTMAT copy (Phase 2.4).

      ALLOCATE(DRTMAT(NELEM_TOTAL,6*NBODY,NSYS))
      DRTMAT=DRMAT

      IF (ISOLV .EQ. 1) THEN
! Direct LU (ISOLV=1). NSYS is 1 or 2 — let MKL's internal threading do the work in ZGETRS.
         DO IP=1, NSYS
              CALL ZGETRS( 'No transpose', NELEM_TOTAL, 6*NBODY, CMAT(:,:,IP), NELEM_TOTAL, IPIV(:,IP), DRTMAT(:,:,IP), NELEM_TOTAL, INFO )
         ENDDO
      ELSE
! GMRES (ISOLV=2) with block-diagonal LU preconditioner on the normal-equations matrix.
         ALLOCATE(CWORK(NELEM_TOTAL,NELEM_TOTAL))
         DO IP=1, NSYS
              CWORK = CMAT(:,:,IP)
              CALL SETUP_BLOCK_PRECONDITIONER(NELEM_TOTAL, CWORK, NBODY, NELEM_MULTI)
              CALL ZITER_SOLVE_MULTI(ISOLV, NELEM_TOTAL, CWORK, DRMAT(:,:,IP), DRTMAT(:,:,IP), 6*NBODY, INFO)
         ENDDO
         DEALLOCATE(CWORK)
      ENDIF
      
      DO 1400 MD=1, 6*NBODY
      DO 1400 IP=1, NSYS
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,IS)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.
      DO 1450 IEL=1, NELEM_TOTAL
         MXPOT(IEL,MD,IP)=CMPLX(0.0D0, 0.0D0)
      DO 1460 IS=1, NSYS
         ! Both ISX/ISY branches were identical — collapsed.
         MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+DRTMAT(IEL,MD,IS)*RXY_MULTI(IP,IS)
1460  CONTINUE
         MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)/NSYS
1450  CONTINUE
!$OMP END DO NOWAIT
!$OMP END PARALLEL
1400  CONTINUE
      
      DEALLOCATE(DRTMAT)

      RETURN
      END SUBROUTINE RADIATION_SOLVER_IRR_MULTI
      
! -----------------------------------------------------------------------------
!     Solve the diffraction potentials after removing irregular frequencies
! -----------------------------------------------------------------------------

      SUBROUTINE DIFFRACTION_SOLVER_IRR_MULTI(CMAT,DDMAT,IPIV,MXPOT,NELEM_TOTAL,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM_TOTAL,NSYS
      INTEGER,INTENT(IN):: IPIV(NELEM_TOTAL,NSYS)
      COMPLEX*16,INTENT(IN):: CMAT(NELEM_TOTAL,NELEM_TOTAL,NSYS),DDMAT(NELEM_TOTAL,NSYS)
      COMPLEX*16,INTENT(OUT):: MXPOT(NELEM_TOTAL,6*NBODY+1,NSYS)
      
      INTEGER IEL,IS,IP,MD,INFO
      COMPLEX*16,ALLOCATABLE:: DDTMAT(:,:)
      COMPLEX*16,ALLOCATABLE:: CWORK(:,:), BSING(:,:), XSING(:,:)
      ! ZGETRS does not modify its A argument — pass CMAT (the LU factorization) directly,
      ! eliminating the previous full-size CTMAT copy (Phase 2.4).

      ALLOCATE(DDTMAT(NELEM_TOTAL,NSYS))
      DDTMAT=DDMAT

      MD=6*NBODY+1

      IF (ISOLV .EQ. 1) THEN
! Direct LU (ISOLV=1). NSYS is 1 or 2 — let MKL's internal threading do the work in ZGETRS.
         DO IP=1, NSYS
              CALL ZGETRS( 'No transpose', NELEM_TOTAL, 1, CMAT(:,:,IP), NELEM_TOTAL, IPIV(:,IP), DDTMAT(:,IP), NELEM_TOTAL, INFO )
         ENDDO
      ELSE
! GMRES (ISOLV=2) with block-diagonal LU preconditioner on the normal-equations matrix.
         ALLOCATE(CWORK(NELEM_TOTAL,NELEM_TOTAL), BSING(NELEM_TOTAL,1), XSING(NELEM_TOTAL,1))
         DO IP=1, NSYS
              CWORK = CMAT(:,:,IP)
              BSING(:,1) = DDMAT(:,IP)
              CALL SETUP_BLOCK_PRECONDITIONER(NELEM_TOTAL, CWORK, NBODY, NELEM_MULTI)
              CALL ZITER_SOLVE_MULTI(ISOLV, NELEM_TOTAL, CWORK, BSING, XSING, 1, INFO)
              DDTMAT(:,IP) = XSING(:,1)
         ENDDO
         DEALLOCATE(CWORK, BSING, XSING)
      ENDIF
      
       DO 400 IP=1, NSYS
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,IS)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.
       DO 450 IEL=1, NELEM_TOTAL
         MXPOT(IEL,MD,IP)=CMPLX(0.0D0, 0.0D0)
       DO 460 IS=1, NSYS
         ! Both ISX/ISY branches were identical — collapsed.
         MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+DDTMAT(IEL,IS)*RXY_MULTI(IP,IS)
460    CONTINUE
         MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)/NSYS
450    CONTINUE
!$OMP END DO NOWAIT
!$OMP END PARALLEL
400    CONTINUE

      DEALLOCATE(DDTMAT)

      RETURN
      END SUBROUTINE DIFFRACTION_SOLVER_IRR_MULTI

!-------------------------------------------------------------------------------
END MODULE AssbMatx_irrMulti
!*******************************************************************************
