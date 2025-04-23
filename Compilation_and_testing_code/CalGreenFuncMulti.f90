!  ------------------------------------------------------------------------------------------------------
!                                                               
!    Program HAMS for the diffraction and radiation of waves 
!    by 3D structures.
! 
!  License:
! 
!    This routine is part of HAMS.
!
!    HAMS is a free software framework: you can redistribute it and/or modify it 
!    under the terms of the Apache License, Version 2.0 (the "License"); you may 
!    not use these subroutines except in compliance with the License. The software
!    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND.
!
!    You should have received a copy of the Apache License, Version 2.0, along with 
!    HAMS. If not, see <http://www.apache.org/licenses/LICENSE-2.0>.
!
!  Code Original Author:
!
!    Yingyi Liu
!
!  ------------------------------------------------------------------------------------------------------
MODULE CalGreenFuncMulti

   USE HAMS_mod
   USE Body_mod
   USE Const_mod
   USE WaveDyn_mod
   USE PanelMesh_mod
   USE Inerfs_mod
   USE SingularIntgrMulti
   USE INFG3D_Open
   USE FinGreen3D_Open
   USE Potentials_mod
   USE omp_lib
   IMPLICIT NONE

  ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: CALGREEN_MULTI
   PUBLIC :: CALGREEN_MULTI_IRR

   
CONTAINS
!   ----------------------------------------------------------------------------
!      Calculate WAVE TERM's value of Green function for arbitrary two points                              ! This is without irregular frequencies
!   ----------------------------------------------------------------------------

      SUBROUTINE CALGREEN_MULTI
      IMPLICIT NONE 
      
      INTEGER  IEL,JEL,IS,IP,IRR,FLAG,BODY_N,NELEM_GLOBAL,NTND_GLOBAL
      REAL*8   XQ,YQ,ZQ,XP,YP,ZP,SIJ,DIJ(3),DIST
      COMPLEX*16  GRN(4),DPOX,DPOY,DPOZ

      IRR=1
      NTND_GLOBAL=0                       ! Variable is used to track the number of nodes in the combined coordinate matrix that have been filled in with the global coordinates
      NELEM_GLOBAL=0                      ! Variable is used to track the number of elements in the combined coordinate matrix that have been filled in with the global coordinates
      
      !Filling the Global coordinate array (this is derived from the coordinates of the panel centers), PNSZ, NCON and NCN. The explicit length for the quantities on the right side (example XYZ_GLOBAL_MULTI(BODY_N,1:NTND_MULTI(BODY_N),:)) is done since the number of panels in different bodies can be different.
      DO BODY_N=1,NBODY
        IF (BODY_N.EQ.1) THEN
         XYZ_GLOBAL_MULTI_COMB(1:NTND_MULTI(BODY_N),:)=XYZ_GLOBAL_MULTI(BODY_N,1:NTND_MULTI(BODY_N),:)
         XYZ_GLOBAL_MULTI_COMB_P(1:NELEM_MULTI(BODY_N),:)=XYZ_MULTI_P(BODY_N,1:NELEM_MULTI(BODY_N),:)
         PNSZ_MULTI_COMB(1:NELEM_MULTI(BODY_N))=PNSZ_MULTI(BODY_N,1:NELEM_MULTI(BODY_N))
         NCN_MULTI_COMB(1:NELEM_MULTI(BODY_N))=NCN_MULTI(BODY_N,1:NELEM_MULTI(BODY_N))
         NCON_MULTI_COMB(1:NELEM_MULTI(BODY_N),:)=NCON_MULTI(BODY_N,1:NELEM_MULTI(BODY_N),:)
        ELSEIF (BODY_N.EQ.NBODY) THEN
         XYZ_GLOBAL_MULTI_COMB(NTND_GLOBAL+1:TNTND,:)=XYZ_GLOBAL_MULTI(BODY_N,1:NTND_MULTI(BODY_N),:)
         XYZ_GLOBAL_MULTI_COMB_P(NELEM_GLOBAL+1:TNELEM,:)=XYZ_MULTI_P(BODY_N,1:NELEM_MULTI(BODY_N),:)
         PNSZ_MULTI_COMB(NELEM_GLOBAL+1:TNELEM)=PNSZ_MULTI(BODY_N,1:NELEM_MULTI(BODY_N))
         NCN_MULTI_COMB(NELEM_GLOBAL+1:TNELEM)=NCN_MULTI(BODY_N,1:NELEM_MULTI(BODY_N))
         NCON_MULTI_COMB(NELEM_GLOBAL+1:TNELEM,:)=NCON_MULTI(BODY_N,1:NELEM_MULTI(BODY_N),:)
        ELSE
         XYZ_GLOBAL_MULTI_COMB(NTND_GLOBAL+1:NTND_GLOBAL+NTND_MULTI(BODY_N),:)=XYZ_GLOBAL_MULTI(BODY_N,1:NTND_MULTI(BODY_N),:)
         XYZ_GLOBAL_MULTI_COMB_P(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N),:)=XYZ_MULTI_P(BODY_N,1:NELEM_MULTI(BODY_N),:)
         PNSZ_MULTI_COMB(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N))=PNSZ_MULTI(BODY_N,1:NELEM_MULTI(BODY_N))
         NCN_MULTI_COMB(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N))=NCN_MULTI(BODY_N,1:NELEM_MULTI(BODY_N))
         NCON_MULTI_COMB(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N),:)=NCON_MULTI(BODY_N,1:NELEM_MULTI(BODY_N),:)
        ENDIF
        NTND_GLOBAL=NTND_GLOBAL+NTND_MULTI(BODY_N)
        NELEM_GLOBAL=NELEM_GLOBAL+NELEM_MULTI(BODY_N)
      ENDDO
      
      ! The panel node numbers need to be modified for the combination array since the node numbers in the first and second body for example can be the same if the bodies are exactly the same
      NELEM_GLOBAL=0
      NTND_GLOBAL=0
      DO BODY_N=1,NBODY
        IF (BODY_N.EQ.1) THEN
         NCON_MULTI_COMB(1:NELEM_MULTI(BODY_N),:)=NCON_MULTI_COMB(1:NELEM_MULTI(BODY_N),:)
        ELSEIF (BODY_N.EQ.NBODY) THEN
         NCON_MULTI_COMB(NELEM_GLOBAL+1:TNELEM,:)=NTND_GLOBAL+NCON_MULTI_COMB(NELEM_GLOBAL+1:TNELEM,:)
        ELSE
         NCON_MULTI_COMB(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N),:)=NTND_GLOBAL+NCON_MULTI_COMB(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N),:)
        ENDIF
        NTND_GLOBAL=NTND_GLOBAL+NTND_MULTI(BODY_N)
        NELEM_GLOBAL=NELEM_GLOBAL+NELEM_MULTI(BODY_N)
      ENDDO
      

!$OMP PARALLEL NUM_THREADS(NTHREAD)                                                                        ! This is the syntax for parallelization
!$OMP DO PRIVATE(IEL,JEL,IS,XQ,XP,YQ,YP,ZQ,ZP,SIJ,DIJ,GRN,FLAG)

      DO IEL=1, TNELEM

        DO JEL=1, TNELEM
            
         XQ=XYZ_GLOBAL_MULTI_COMB_P(JEL,1)     ! JEL: source point,  IEL: field point
         YQ=XYZ_GLOBAL_MULTI_COMB_P(JEL,2)
         ZQ=XYZ_GLOBAL_MULTI_COMB_P(JEL,3)

         DIST=SQRT((XYZ_GLOBAL_MULTI_COMB_P(IEL,1)-XYZ_GLOBAL_MULTI_COMB_P(JEL,1))**2+(XYZ_GLOBAL_MULTI_COMB_P(IEL,2)-XYZ_GLOBAL_MULTI_COMB_P(JEL,2))**2+(XYZ_GLOBAL_MULTI_COMB_P(IEL,3)-XYZ_GLOBAL_MULTI_COMB_P(JEL,3))**2)
         IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
          FLAG=1             ! What are these flags for?
         ELSE
          FLAG=0
         ENDIF
         
        DO IS=1, NSYS

         IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
          XP=SY_MULTI(IS,1)*XYZ_GLOBAL_MULTI_COMB_P(IEL,1)
          YP=SX_MULTI(IS,1)*XYZ_GLOBAL_MULTI_COMB_P(IEL,2)
          ZP=         XYZ_GLOBAL_MULTI_COMB_P(IEL,3)
         ELSE                                       
          XP=SX_MULTI(IS,1)*XYZ_GLOBAL_MULTI_COMB_P(IEL,1)
          YP=SY_MULTI(IS,1)*XYZ_GLOBAL_MULTI_COMB_P(IEL,2)
          ZP=         XYZ_GLOBAL_MULTI_COMB_P(IEL,3)
         ENDIF

         IF (NCN_MULTI_COMB(JEL).EQ.3) THEN

          CALL SGLINTBD_MULTI_TRI(IS,IEL,JEL,SIJ,DIJ,1)

         ELSEIF (NCN_MULTI_COMB(JEL).EQ.4) THEN

          CALL SGLINTBD_MULTI_QUAD(IS,IEL,JEL,SIJ,DIJ,1)

         ENDIF

         IF (H.LT.0.D0) THEN
           CALL INFGREEN3D(XQ,XP,YQ,YP,ZQ,ZP,V,GRN,FLAG)
         ELSE
          CALL FINGREEN3D(XQ,XP,YQ,YP,ZQ,ZP,V,WVN,NK,H,GRN,FLAG)
         ENDIF

         RKBN_MULTI_COMB(IEL,JEL,IS,1)=SIJ
         RKBN_MULTI_COMB(IEL,JEL,IS,2)=DIJ(1)
         RKBN_MULTI_COMB(IEL,JEL,IS,3)=DIJ(2)
         RKBN_MULTI_COMB(IEL,JEL,IS,4)=DIJ(3)
         CGRN_MULTI_COMB(IEL,JEL,IS,:)=GRN(:)
      
        ENDDO
        ENDDO
      ENDDO
        
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      RETURN
      END SUBROUTINE CALGREEN_MULTI 
      
!   ----------------------------------------------------------------------------
!      Calculate WAVE TERM's value of Green function for arbitray two points                                      ! This is with Irregular frequencies? IRR is 1 here as well?
!   ----------------------------------------------------------------------------

      SUBROUTINE CALGREEN_MULTI_IRR
      IMPLICIT NONE
      
      INTEGER  IEL,JEL,IS,IP,IRR,FLAG,BODY_N,NELEM_GLOBAL,NTND_GLOBAL,iNELEM_GLOBAL,iNTND_GLOBAL
      REAL*8   XQ,YQ,ZQ,XP,YP,ZP,SIJ,DIJ(3),DIST
      COMPLEX*16  GRN(4),DPOX,DPOY,DPOZ

      NTND_GLOBAL=0                       ! Variable is used to track the number of nodes in the combined coordinate matrix that have been filled in with the global coordinates
      NELEM_GLOBAL=0                      ! Variable is used to track the number of elements in the combined coordinate matrix that have been filled in with the global coordinates
      iNTND_GLOBAL=0                       ! Variable is used to track the number of nodes (water plane mesh) in the combined coordinate matrix that have been filled in with the global coordinates
      iNELEM_GLOBAL=0                      ! Variable is used to track the number of elements (water plane mesh) in the combined coordinate matrix that have been filled in with the global coordinates
      
      !Filling the Global coordinate array (this is derived from the coordinates of the panel centers), PNSZ, NCON and NCN. The explicit length for the quantities on the right side (example XYZ_GLOBAL_MULTI(BODY_N,1:NTND_MULTI(BODY_N),:)) is done since the number of panels in different bodies can be different.
      DO BODY_N=1,NBODY
        IF (BODY_N.EQ.1) THEN
         XYZ_GLOBAL_MULTI_COMB(1:NTND_MULTI(BODY_N),:)=XYZ_GLOBAL_MULTI(BODY_N,1:NTND_MULTI(BODY_N),:)
         XYZ_GLOBAL_MULTI_COMB_P(1:NELEM_MULTI(BODY_N),:)=XYZ_MULTI_P(BODY_N,1:NELEM_MULTI(BODY_N),:)
         PNSZ_MULTI_COMB(1:NELEM_MULTI(BODY_N))=PNSZ_MULTI(BODY_N,1:NELEM_MULTI(BODY_N))
         NCN_MULTI_COMB(1:NELEM_MULTI(BODY_N))=NCN_MULTI(BODY_N,1:NELEM_MULTI(BODY_N))
         NCON_MULTI_COMB(1:NELEM_MULTI(BODY_N),:)=NCON_MULTI(BODY_N,1:NELEM_MULTI(BODY_N),:)
         IF (iNELEM_MULTI(BODY_N).GT.0) THEN
          iXYZ_GLOBAL_MULTI_COMB(1:iNTND_MULTI(BODY_N),:)=iXYZ_GLOBAL_MULTI(BODY_N,1:iNTND_MULTI(BODY_N),:)
          iXYZ_GLOBAL_MULTI_COMB_P(1:iNELEM_MULTI(BODY_N),:)=iXYZ_MULTI_P(BODY_N,1:iNELEM_MULTI(BODY_N),:)
          iPNSZ_MULTI_COMB(1:iNELEM_MULTI(BODY_N))=iPNSZ_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N))
          iNCN_MULTI_COMB(1:iNELEM_MULTI(BODY_N))=iNCN_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N))
          iNCON_MULTI_COMB(1:iNELEM_MULTI(BODY_N),:)=iNCON_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N),:)
         ENDIF
        ELSEIF (BODY_N.EQ.NBODY) THEN
         XYZ_GLOBAL_MULTI_COMB(NTND_GLOBAL+1:NTND_TOTAL,:)=XYZ_GLOBAL_MULTI(BODY_N,1:NTND_MULTI(BODY_N),:)
         XYZ_GLOBAL_MULTI_COMB_P(NELEM_GLOBAL+1:NELEM_TOTAL,:)=XYZ_MULTI_P(BODY_N,1:NELEM_MULTI(BODY_N),:)
         PNSZ_MULTI_COMB(NELEM_GLOBAL+1:NELEM_TOTAL)=PNSZ_MULTI(BODY_N,1:NELEM_MULTI(BODY_N))
         NCN_MULTI_COMB(NELEM_GLOBAL+1:NELEM_TOTAL)=NCN_MULTI(BODY_N,1:NELEM_MULTI(BODY_N))
         NCON_MULTI_COMB(NELEM_GLOBAL+1:NELEM_TOTAL,:)=NCON_MULTI(BODY_N,1:NELEM_MULTI(BODY_N),:)
         IF (iNELEM_MULTI(BODY_N).GT.0) THEN
          iXYZ_GLOBAL_MULTI_COMB(iNTND_GLOBAL+1:iNTND_TOTAL,:)=iXYZ_GLOBAL_MULTI(BODY_N,1:iNTND_MULTI(BODY_N),:)
          iXYZ_GLOBAL_MULTI_COMB_P(iNELEM_GLOBAL+1:iNELEM_TOTAL,:)=iXYZ_MULTI_P(BODY_N,1:iNELEM_MULTI(BODY_N),:)
          iPNSZ_MULTI_COMB(iNELEM_GLOBAL+1:iNELEM_TOTAL)=iPNSZ_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N))
          iNCN_MULTI_COMB(iNELEM_GLOBAL+1:iNELEM_TOTAL)=iNCN_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N))
          iNCON_MULTI_COMB(iNELEM_GLOBAL+1:iNELEM_TOTAL,:)=iNCON_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N),:)
         ENDIF
        ELSE
         XYZ_GLOBAL_MULTI_COMB(NTND_GLOBAL+1:NTND_GLOBAL+NTND_MULTI(BODY_N),:)=XYZ_GLOBAL_MULTI(BODY_N,1:NTND_MULTI(BODY_N),:)
         XYZ_GLOBAL_MULTI_COMB_P(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N),:)=XYZ_MULTI_P(BODY_N,1:NELEM_MULTI(BODY_N),:)
         PNSZ_MULTI_COMB(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N))=PNSZ_MULTI(BODY_N,1:NELEM_MULTI(BODY_N))
         NCN_MULTI_COMB(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N))=NCN_MULTI(BODY_N,1:NELEM_MULTI(BODY_N))
         NCON_MULTI_COMB(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N),:)=NCON_MULTI(BODY_N,1:NELEM_MULTI(BODY_N),:)
         IF (iNELEM_MULTI(BODY_N).GT.0) THEN
          iXYZ_GLOBAL_MULTI_COMB(iNTND_GLOBAL+1:iNTND_GLOBAL+iNTND_MULTI(BODY_N),:)=iXYZ_GLOBAL_MULTI(BODY_N,1:iNTND_MULTI(BODY_N),:)
          iXYZ_GLOBAL_MULTI_COMB_P(iNELEM_GLOBAL+1:iNELEM_GLOBAL+iNELEM_MULTI(BODY_N),:)=iXYZ_MULTI_P(BODY_N,1:iNELEM_MULTI(BODY_N),:)
          iPNSZ_MULTI_COMB(iNELEM_GLOBAL+1:iNELEM_GLOBAL+iNELEM_MULTI(BODY_N))=iPNSZ_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N))
          iNCN_MULTI_COMB(iNELEM_GLOBAL+1:iNELEM_GLOBAL+iNELEM_MULTI(BODY_N))=iNCN_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N))
          iNCON_MULTI_COMB(iNELEM_GLOBAL+1:iNELEM_GLOBAL+iNELEM_MULTI(BODY_N),:)=iNCON_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N),:)   
         ENDIF
        ENDIF
        NTND_GLOBAL=NTND_GLOBAL+NTND_MULTI(BODY_N)
        NELEM_GLOBAL=NELEM_GLOBAL+NELEM_MULTI(BODY_N)
        iNTND_GLOBAL=iNTND_GLOBAL+iNTND_MULTI(BODY_N)
        iNELEM_GLOBAL=iNELEM_GLOBAL+iNELEM_MULTI(BODY_N)
      ENDDO
      
      ! The panel node numbers need to be modified for the combination array since the node numbers in the first and second body for example can be the same if the bodies are exactly the same
      NELEM_GLOBAL=0
      NTND_GLOBAL=0
      iNELEM_GLOBAL=0
      iNTND_GLOBAL=0
      DO BODY_N=1,NBODY
        IF (BODY_N.EQ.1) THEN
         NCON_MULTI_COMB(1:NELEM_MULTI(BODY_N),:)=NCON_MULTI_COMB(1:NELEM_MULTI(BODY_N),:)
         IF (iNELEM_MULTI(BODY_N).GT.0) THEN
          iNCON_MULTI_COMB(1:iNELEM_MULTI(BODY_N),:)=iNCON_MULTI_COMB(1:iNELEM_MULTI(BODY_N),:)
         ENDIF
        ELSEIF (BODY_N.EQ.NBODY) THEN
         NCON_MULTI_COMB(NELEM_GLOBAL+1:NELEM_TOTAL,:)=NTND_GLOBAL+NCON_MULTI_COMB(NELEM_GLOBAL+1:NELEM_TOTAL,:)
         IF (iNELEM_MULTI(BODY_N).GT.0) THEN
          iNCON_MULTI_COMB(iNELEM_GLOBAL+1:iNELEM_TOTAL,:)=iNTND_GLOBAL+iNCON_MULTI_COMB(iNELEM_GLOBAL+1:iNELEM_TOTAL,:)
         ENDIF
        ELSE
         NCON_MULTI_COMB(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N),:)=NTND_GLOBAL+NCON_MULTI_COMB(NELEM_GLOBAL+1:NELEM_GLOBAL+NELEM_MULTI(BODY_N),:)
         IF (iNELEM_MULTI(BODY_N).GT.0) THEN
          iNCON_MULTI_COMB(iNELEM_GLOBAL+1:iNELEM_GLOBAL+iNELEM_MULTI(BODY_N),:)=iNTND_GLOBAL+iNCON_MULTI_COMB(iNELEM_GLOBAL+1:iNELEM_GLOBAL+iNELEM_MULTI(BODY_N),:)
         ENDIF
        ENDIF
        NTND_GLOBAL=NTND_GLOBAL+NTND_MULTI(BODY_N)
        NELEM_GLOBAL=NELEM_GLOBAL+NELEM_MULTI(BODY_N)
        iNTND_GLOBAL=iNTND_GLOBAL+iNTND_MULTI(BODY_N)
        iNELEM_GLOBAL=iNELEM_GLOBAL+iNELEM_MULTI(BODY_N)
      ENDDO
      
      IRR=1
!$OMP PARALLEL NUM_THREADS(NTHREAD)                                                                        ! This is the syntax for parallelization
!$OMP DO PRIVATE(IEL,JEL,IS,XQ,XP,YQ,YP,ZQ,ZP,SIJ,DIJ,GRN,FLAG)

      DO IEL=1, NELEM_TOTAL

        DO JEL=1, NELEM_TOTAL
            
         XQ=XYZ_GLOBAL_MULTI_COMB_P(JEL,1)     ! JEL: source point,  IEL: field point
         YQ=XYZ_GLOBAL_MULTI_COMB_P(JEL,2)
         ZQ=XYZ_GLOBAL_MULTI_COMB_P(JEL,3)

         DIST=SQRT((XYZ_GLOBAL_MULTI_COMB_P(IEL,1)-XYZ_GLOBAL_MULTI_COMB_P(JEL,1))**2+(XYZ_GLOBAL_MULTI_COMB_P(IEL,2)-XYZ_GLOBAL_MULTI_COMB_P(JEL,2))**2+(XYZ_GLOBAL_MULTI_COMB_P(IEL,3)-XYZ_GLOBAL_MULTI_COMB_P(JEL,3))**2)
         IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
          FLAG=1             ! What are these flags for?
         ELSE
          FLAG=0
         ENDIF
         
        DO IS=1, NSYS

         IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
          XP=SY_MULTI(IS,1)*XYZ_GLOBAL_MULTI_COMB_P(IEL,1)
          YP=SX_MULTI(IS,1)*XYZ_GLOBAL_MULTI_COMB_P(IEL,2)
          ZP=         XYZ_GLOBAL_MULTI_COMB_P(IEL,3)
         ELSE                                       
          XP=SX_MULTI(IS,1)*XYZ_GLOBAL_MULTI_COMB_P(IEL,1)
          YP=SY_MULTI(IS,1)*XYZ_GLOBAL_MULTI_COMB_P(IEL,2)
          ZP=         XYZ_GLOBAL_MULTI_COMB_P(IEL,3)
         ENDIF

         IF (NCN_MULTI_COMB(JEL).EQ.3) THEN

          CALL SGLINTBD_MULTI_TRI(IS,IEL,JEL,SIJ,DIJ,IRR)

         ELSEIF (NCN_MULTI_COMB(JEL).EQ.4) THEN

          CALL SGLINTBD_MULTI_QUAD(IS,IEL,JEL,SIJ,DIJ,IRR)

         ENDIF

         IF (H.LT.0.D0) THEN
           CALL INFGREEN3D(XQ,XP,YQ,YP,ZQ,ZP,V,GRN,FLAG)
         ELSE
          CALL FINGREEN3D(XQ,XP,YQ,YP,ZQ,ZP,V,WVN,NK,H,GRN,FLAG)
         ENDIF

         RKBN_MULTI_COMB(IEL,JEL,IS,1)=SIJ
         RKBN_MULTI_COMB(IEL,JEL,IS,2)=DIJ(1)
         RKBN_MULTI_COMB(IEL,JEL,IS,3)=DIJ(2)
         RKBN_MULTI_COMB(IEL,JEL,IS,4)=DIJ(3)
         CGRN_MULTI_COMB(IEL,JEL,IS,:)=GRN(:)
      
        ENDDO
        ENDDO
      ENDDO
        
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IRR=3

!$OMP PARALLEL DO PRIVATE(IEL,JEL,IS,XQ,XP,YQ,YP,ZQ,ZP,SIJ,DIJ,GRN,FLAG)

      DO IEL=1, iNELEM_TOTAL
          
        DO IS=1, NSYS

        IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
         XP=SY_MULTI(IS,1)*iXYZ_GLOBAL_MULTI_COMB_P(IEL,1)
         YP=SX_MULTI(IS,1)*iXYZ_GLOBAL_MULTI_COMB_P(IEL,2)
         ZP=          iXYZ_GLOBAL_MULTI_COMB_P(IEL,3)
        ELSE
         XP=SX_MULTI(IS,1)*iXYZ_GLOBAL_MULTI_COMB_P(IEL,1)
         YP=SY_MULTI(IS,1)*iXYZ_GLOBAL_MULTI_COMB_P(IEL,2)
         ZP=         iXYZ_GLOBAL_MULTI_COMB_P(IEL,3)
        ENDIF

        DO JEL=1, NELEM_TOTAL

         XQ=XYZ_GLOBAL_MULTI_COMB_P(JEL,1)
         YQ=XYZ_GLOBAL_MULTI_COMB_P(JEL,2)
         ZQ=XYZ_GLOBAL_MULTI_COMB_P(JEL,3)
         
         DIST=SQRT((XP-XQ)**2+(YP-YQ)**2+(ZP-ZQ)**2)
         IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
          FLAG=1
         ELSE
          FLAG=0
         ENDIF
         
         IF (NCN_MULTI_COMB(JEL).EQ.3) THEN
             
          CALL SGLINTBD_MULTI_TRI(IS,IEL,JEL,SIJ,DIJ,IRR)
          
         ELSEIF (NCN_MULTI_COMB(JEL).EQ.4) THEN
             
          CALL SGLINTBD_MULTI_QUAD(IS,IEL,JEL,SIJ,DIJ,IRR)

         ENDIF 

         IF (H.LT.0.D0) THEN
           CALL INFGREEN3D(XQ,XP,YQ,YP,ZQ,ZP,V,GRN,FLAG)
         ELSE
          CALL FINGREEN3D(XQ,XP,YQ,YP,ZQ,ZP,V,WVN,NK,H,GRN,FLAG)
         ENDIF

         PKBN_MULTI_COMB(IEL,JEL,IS,1)=SIJ
         PKBN_MULTI_COMB(IEL,JEL,IS,2)=DIJ(1)
         PKBN_MULTI_COMB(IEL,JEL,IS,3)=DIJ(2)
         PKBN_MULTI_COMB(IEL,JEL,IS,4)=DIJ(3)
         DGRN_MULTI_COMB(IEL,JEL,IS,:)=GRN(:)
         
        ENDDO
        ENDDO
      ENDDO

!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE CALGREEN_MULTI_IRR
!-------------------------------------------------------------------------------
END MODULE CalGreenFuncMulti
!*******************************************************************************