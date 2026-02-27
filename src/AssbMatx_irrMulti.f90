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
      INTEGER:: BODY_N,NELEM_GLOBAL,iNELEM_GLOBAL
      INTEGER,INTENT(IN):: NELEM_TOTAL,iNELEM_TOTAL,TNELEM,NSYS
      INTEGER,INTENT(OUT):: IPIV(NELEM_TOTAL,NSYS)
      COMPLEX*16,INTENT(OUT):: AMAT(TNELEM,TNELEM,NSYS),CMAT(NELEM_TOTAL,NELEM_TOTAL,NSYS)
      
      INTEGER IEL,JEL,KEL,IS,IP,IRR,FLAG,INFO
      REAL*8 DIST
      COMPLEX*16 TINDP(4),SUM
      
      NELEM_GLOBAL=0
      iNELEM_GLOBAL=0
      
      ! Creating combination of the variables required for future computations
      DO BODY_N=1,NBODY
        IF (BODY_N.EQ.1) THEN
         DS_MULTI_COMB(1:NELEM_MULTI(BODY_N))=DS_MULTI(BODY_N,:)
         DXYZ_MULTI_COMB(1:NELEM_MULTI(BODY_N),:)=DXYZ_MULTI_P(BODY_N,:,:)
         IF (iNELEM_MULTI(BODY_N).GT.0) THEN
          iDS_MULTI_COMB(1:iNELEM_MULTI(BODY_N))=iDS_MULTI(BODY_N,:)
          iDXYZ_MULTI_COMB(1:iNELEM_MULTI(BODY_N),:)=iDXYZ_MULTI_P(BODY_N,:,:)
         ENDIF
        ELSEIF (BODY_N.EQ.NBODY) THEN
         DS_MULTI_COMB(NELEM_GLOBAL+1:NELEM_TOTAL)=DS_MULTI(BODY_N,:)
         DXYZ_MULTI_COMB(NELEM_GLOBAL+1:NELEM_TOTAL,:)=DXYZ_MULTI_P(BODY_N,:,:)
         IF (iNELEM_MULTI(BODY_N).GT.0) THEN
          iDS_MULTI_COMB(iNELEM_GLOBAL+1:iNELEM_TOTAL)=iDS_MULTI(BODY_N,:)
          iDXYZ_MULTI_COMB(iNELEM_GLOBAL+1:iNELEM_TOTAL,:)=iDXYZ_MULTI_P(BODY_N,:,:)
         ENDIF
        ELSE
         DS_MULTI_COMB(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N))=DS_MULTI(BODY_N,:)
         DXYZ_MULTI_COMB(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N),:)=DXYZ_MULTI_P(BODY_N,:,:)
         IF (iNELEM_MULTI(BODY_N).GT.0) THEN
          iDS_MULTI_COMB(iNELEM_GLOBAL+1:iNELEM_GLOBAL+iNELEM_MULTI(BODY_N))=iDS_MULTI(BODY_N,:)
          iDXYZ_MULTI_COMB(iNELEM_GLOBAL+1:iNELEM_GLOBAL+iNELEM_MULTI(BODY_N),:)=iDXYZ_MULTI_P(BODY_N,:,:)   
         ENDIF
        ENDIF
        NELEM_GLOBAL=NELEM_GLOBAL+NELEM_MULTI(BODY_N)
        iNELEM_GLOBAL=iNELEM_GLOBAL+iNELEM_MULTI(BODY_N)
      ENDDO

      AMAT=CMPLX(0.0D0,0.0D0)

      IRR=1

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,JEL,IP,IS,FLAG,DIST,TINDP) !$OMP REDUCTION(+:AMAT)

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
       
      IRR=3

!$OMP PARALLEL NUM_THREADS(NTHREAD)  
!$OMP DO PRIVATE(IEL,JEL,IP,IS,FLAG,DIST,TINDP) !$OMP REDUCTION(+:AMAT)
      
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
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           AMAT(IEL,JEL,IP)=AMAT(IEL,JEL,IP)+RXY_MULTI(IS,IP)*TINDP(IS)
          ELSE
           AMAT(IEL,JEL,IP)=AMAT(IEL,JEL,IP)+RXY_MULTI(IS,IP)*TINDP(IS)
          ENDIF
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
!$OMP DO PRIVATE(IEL,JEL,KEL,SUM) !$OMP REDUCTION(+:SUM)
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

!$omp parallel do private(NTHREAD)
       DO IP=1, NSYS
          CALL ZGETRF( NELEM_TOTAL, NELEM_TOTAL, CMAT(:,:,IP), NELEM_TOTAL, IPIV(:,IP), INFO )
       ENDDO
!$omp end parallel do

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
!$OMP DO PRIVATE(IEL,JEL,MD,IP,IS,FLAG,DIST,TINRD,BTMP) !$OMP REDUCTION(+:BRMAT)
      
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
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           BRMAT(IEL,MD,IP)=BRMAT(IEL,MD,IP)+RXY_MULTI(IP,IS)*BTMP(MD,IS)
          ELSE
           BRMAT(IEL,MD,IP)=BRMAT(IEL,MD,IP)+RXY_MULTI(IP,IS)*BTMP(MD,IS)
          ENDIF
300     CONTINUE
     
1000    CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL
       
      IRR=3
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,JEL,MD,IP,IS,FLAG,DIST,TINRD,BTMP) !$OMP REDUCTION(+:BRMAT)
      
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
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           BRMAT(IEL,MD,IP)=BRMAT(IEL,MD,IP)+RXY_MULTI(IP,IS)*BTMP(MD,IS)
          ELSE
           BRMAT(IEL,MD,IP)=BRMAT(IEL,MD,IP)+RXY_MULTI(IP,IS)*BTMP(MD,IS)
          ENDIF
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
!$OMP DO PRIVATE(IP,IEL,KEL,MD) !$OMP REDUCTION(+:DRMAT)  

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
!$OMP DO PRIVATE(XP,YP,ZP,IEL,JEL,IP,IS,FLAG,DIST,TINRD,BTMP) !$OMP REDUCTION(+:BDMAT)     
      DO  1000 IEL=1,  NELEM_TOTAL

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
300     CONTINUE
     
1000   CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL
       
      IRR=3

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(XP,YP,ZP,IEL,JEL,IP,IS,FLAG,DIST,TINRD,BTMP) !$OMP REDUCTION(+:BDMAT)
      
      DO  3000 IEL=NELEM_TOTAL+1,  NELEM_TOTAL+iNELEM_TOTAL

       BTMP=CMPLX(0.0D0,0.0D0)
       
       IF (ISOL.EQ.2) THEN
            
         DO 500  IP=1,  NSYS
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           XP=RX_MULTI(IP,1)*iXYZ_GLOBAL_MULTI_COMB_P(IEL-NELEM_TOTAL,1)
           YP=RX_MULTI(IP,2)*iXYZ_GLOBAL_MULTI_COMB_P(IEL-NELEM_TOTAL,2)
           ZP=         iXYZ_GLOBAL_MULTI_COMB_P(IEL-NELEM_TOTAL,3)
          ELSE
           XP=RY_MULTI(IP,1)*iXYZ_GLOBAL_MULTI_COMB_P(IEL-NELEM_TOTAL,1)
           YP=RY_MULTI(IP,2)*iXYZ_GLOBAL_MULTI_COMB_P(IEL-NELEM_TOTAL,2)
           ZP=         iXYZ_GLOBAL_MULTI_COMB_P(IEL-NELEM_TOTAL,3)
          ENDIF
          BTMP(IP)=4.D0*PI*VINP(XP,YP,ZP,XW(1),XW(2),BETA)
500      CONTINUE
         
       ELSEIF (ISOL.EQ.1) THEN

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
        
       ELSE
           
        PRINT*,"  Error: The input for ISOL should be either 1 or 2."
        STOP
        
       ENDIF
       
       DO  700  IP=1, NSYS

       DO  700  IS=1, NSYS
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           BDMAT(IEL,IP)=BDMAT(IEL,IP)+RXY_MULTI(IP,IS)*BTMP(IS)
          ELSE
           BDMAT(IEL,IP)=BDMAT(IEL,IP)+RXY_MULTI(IP,IS)*BTMP(IS)
          ENDIF
700    CONTINUE

3000   CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL
!
!--------------------------------------------------------------
!    Matrix set-up for the least-square problem

      DDMAT=CMPLX(0.0D0,0.0D0)

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IP,IEL,KEL) !$OMP REDUCTION(+:DDMAT)  

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
      COMPLEX*16,ALLOCATABLE:: CTMAT(:,:,:),DRTMAT(:,:,:)
      
      ALLOCATE(CTMAT(NELEM_TOTAL,NELEM_TOTAL,NSYS),DRTMAT(NELEM_TOTAL,6*NBODY,NSYS))

      CTMAT=CMAT
      DRTMAT=DRMAT
      
!$omp parallel do private(NTHREAD)
      DO IP=1, NSYS
           CALL ZGETRS( 'No transpose', NELEM_TOTAL, 6*NBODY, CTMAT(:,:,IP), NELEM_TOTAL, IPIV(:,IP), DRTMAT(:,:,IP), NELEM_TOTAL, INFO )
      ENDDO
!$omp end parallel do
      
      DO 1400 MD=1, 6*NBODY
      DO 1400 IP=1, NSYS
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,IS) !$OMP REDUCTION(+:MXPOT)
      DO 1450 IEL=1, NELEM_TOTAL
         MXPOT(IEL,MD,IP)=CMPLX(0.0D0, 0.0D0)
      DO 1460 IS=1, NSYS
         IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
          MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+DRTMAT(IEL,MD,IS)*RXY_MULTI(IP,IS)
         ELSE
          MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+DRTMAT(IEL,MD,IS)*RXY_MULTI(IP,IS)
         ENDIF
1460  CONTINUE
         MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)/NSYS
1450  CONTINUE
!$OMP END DO NOWAIT
!$OMP END PARALLEL
1400  CONTINUE
      
      DEALLOCATE(CTMAT,DRTMAT)
      
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
      COMPLEX*16,ALLOCATABLE:: CTMAT(:,:,:),DDTMAT(:,:)
      
      ALLOCATE(CTMAT(NELEM_TOTAL,NELEM_TOTAL,NSYS),DDTMAT(NELEM_TOTAL,NSYS))

      CTMAT=CMAT
      DDTMAT=DDMAT

      MD=6*NBODY+1
      
!$omp parallel do private(NTHREAD)
      DO IP=1, NSYS
           CALL ZGETRS( 'No transpose', NELEM_TOTAL, 1, CTMAT(:,:,IP), NELEM_TOTAL, IPIV(:,IP), DDTMAT(:,IP), NELEM_TOTAL, INFO )
      ENDDO
!$omp end parallel do
      
       DO 400 IP=1, NSYS
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,IS) !$OMP REDUCTION(+:MXPOT)
       DO 450 IEL=1, NELEM_TOTAL
         MXPOT(IEL,MD,IP)=CMPLX(0.0D0, 0.0D0)
       DO 460 IS=1, NSYS
         IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
          MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+DDTMAT(IEL,IS)*RXY_MULTI(IP,IS)
         ELSE
          MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+DDTMAT(IEL,IS)*RXY_MULTI(IP,IS)
         ENDIF
460    CONTINUE
         MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)/NSYS
450    CONTINUE
!$OMP END DO NOWAIT
!$OMP END PARALLEL
400    CONTINUE

      DEALLOCATE(CTMAT,DDTMAT)
      
      RETURN
      END SUBROUTINE DIFFRACTION_SOLVER_IRR_MULTI

!-------------------------------------------------------------------------------
END MODULE AssbMatx_irrMulti
!*******************************************************************************
