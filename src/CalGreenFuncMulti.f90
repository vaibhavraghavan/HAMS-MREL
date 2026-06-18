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
   PUBLIC :: BuildCombinedArrays


CONTAINS
!   ----------------------------------------------------------------------------
!      Pack the per-body padded geometry into the contiguous *_MULTI_COMB
!      arrays. Called ONCE after the mesh-read / normal-calc setup phase, so
!      the COMB construction (previously repeated every frequency in
!      CALGREEN_MULTI / ASSB_LEFT_MULTI) no longer happens in the hot loop.
!   ----------------------------------------------------------------------------

      SUBROUTINE BuildCombinedArrays
      IMPLICIT NONE

      INTEGER BODY_N, IS_N, IE_N, JS_E, JE_E
      INTEGER NTND_GLOBAL, NELEM_GLOBAL, iNTND_GLOBAL, iNELEM_GLOBAL

      ! Hull mesh — always present.
      NTND_GLOBAL  = 0
      NELEM_GLOBAL = 0
      DO BODY_N = 1, NBODY
         IS_N = NTND_GLOBAL  + 1
         IE_N = NTND_GLOBAL  + NTND_MULTI(BODY_N)
         JS_E = NELEM_GLOBAL + 1
         JE_E = NELEM_GLOBAL + NELEM_MULTI(BODY_N)

         XYZ_GLOBAL_MULTI_COMB  (IS_N:IE_N,:) = XYZ_GLOBAL_MULTI(BODY_N, 1:NTND_MULTI(BODY_N), :)
         XYZ_GLOBAL_MULTI_COMB_P(JS_E:JE_E,:) = XYZ_MULTI_P     (BODY_N, 1:NELEM_MULTI(BODY_N), :)
         PNSZ_MULTI_COMB        (JS_E:JE_E)   = PNSZ_MULTI      (BODY_N, 1:NELEM_MULTI(BODY_N))
         NCN_MULTI_COMB         (JS_E:JE_E)   = NCN_MULTI       (BODY_N, 1:NELEM_MULTI(BODY_N))
         ! NCON values are local node IDs per body; offset by NTND_GLOBAL so they index XYZ_GLOBAL_MULTI_COMB.
         NCON_MULTI_COMB        (JS_E:JE_E,:) = NCON_MULTI      (BODY_N, 1:NELEM_MULTI(BODY_N), :) + NTND_GLOBAL
         DS_MULTI_COMB          (JS_E:JE_E)   = DS_MULTI        (BODY_N, 1:NELEM_MULTI(BODY_N))
         DXYZ_MULTI_COMB        (JS_E:JE_E,:) = DXYZ_MULTI_P    (BODY_N, 1:NELEM_MULTI(BODY_N), :)

         NTND_GLOBAL  = NTND_GLOBAL  + NTND_MULTI(BODY_N)
         NELEM_GLOBAL = NELEM_GLOBAL + NELEM_MULTI(BODY_N)
      ENDDO

      ! Waterplane mesh — only allocated and only present when IRSP > 0.
      IF (IRSP .NE. 0) THEN
         iNTND_GLOBAL  = 0
         iNELEM_GLOBAL = 0
         DO BODY_N = 1, NBODY
            IF (iNELEM_MULTI(BODY_N) .GT. 0) THEN
               IS_N = iNTND_GLOBAL  + 1
               IE_N = iNTND_GLOBAL  + iNTND_MULTI(BODY_N)
               JS_E = iNELEM_GLOBAL + 1
               JE_E = iNELEM_GLOBAL + iNELEM_MULTI(BODY_N)

               iXYZ_GLOBAL_MULTI_COMB  (IS_N:IE_N,:) = iXYZ_GLOBAL_MULTI(BODY_N, 1:iNTND_MULTI(BODY_N), :)
               iXYZ_GLOBAL_MULTI_COMB_P(JS_E:JE_E,:) = iXYZ_MULTI_P     (BODY_N, 1:iNELEM_MULTI(BODY_N), :)
               iPNSZ_MULTI_COMB        (JS_E:JE_E)   = iPNSZ_MULTI      (BODY_N, 1:iNELEM_MULTI(BODY_N))
               iNCN_MULTI_COMB         (JS_E:JE_E)   = iNCN_MULTI       (BODY_N, 1:iNELEM_MULTI(BODY_N))
               iNCON_MULTI_COMB        (JS_E:JE_E,:) = iNCON_MULTI      (BODY_N, 1:iNELEM_MULTI(BODY_N), :) + iNTND_GLOBAL
               iDS_MULTI_COMB          (JS_E:JE_E)   = iDS_MULTI        (BODY_N, 1:iNELEM_MULTI(BODY_N))
               iDXYZ_MULTI_COMB        (JS_E:JE_E,:) = iDXYZ_MULTI_P    (BODY_N, 1:iNELEM_MULTI(BODY_N), :)
            ENDIF
            iNTND_GLOBAL  = iNTND_GLOBAL  + iNTND_MULTI(BODY_N)
            iNELEM_GLOBAL = iNELEM_GLOBAL + iNELEM_MULTI(BODY_N)
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE BuildCombinedArrays

!   ----------------------------------------------------------------------------
!      Calculate WAVE TERM's value of Green function for arbitrary two points                              ! This is without irregular frequencies
!   ----------------------------------------------------------------------------

      SUBROUTINE CALGREEN_MULTI
      IMPLICIT NONE

      INTEGER  IEL,JEL,IS,IP,IRR,FLAG
      REAL*8   XQ,YQ,ZQ,XP,YP,ZP,SIJ,DIJ(3),DIST
      COMPLEX*16  GRN(4),DPOX,DPOY,DPOZ

      ! The per-body→COMB packing previously inlined here is now done once in BuildCombinedArrays
      ! at startup (Phase 2.3) — those copies are frequency-independent.
      IRR=1

!$OMP PARALLEL NUM_THREADS(NTHREAD)                                                                        ! This is the syntax for parallelization
!$OMP DO PRIVATE(IEL,JEL,IS,XQ,XP,YQ,YP,ZQ,ZP,SIJ,DIJ,GRN,FLAG) SCHEDULE(DYNAMIC,16)

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

      INTEGER  IEL,JEL,IS,IP,IRR,FLAG
      REAL*8   XQ,YQ,ZQ,XP,YP,ZP,SIJ,DIJ(3),DIST
      COMPLEX*16  GRN(4),DPOX,DPOY,DPOZ

      ! The per-body→COMB packing previously inlined here is now done once in BuildCombinedArrays
      ! at startup (Phase 2.3) — those copies are frequency-independent.
      IRR=1
!$OMP PARALLEL NUM_THREADS(NTHREAD)                                                                        ! This is the syntax for parallelization
!$OMP DO PRIVATE(IEL,JEL,IS,XQ,XP,YQ,YP,ZQ,ZP,SIJ,DIJ,GRN,FLAG) SCHEDULE(DYNAMIC,16)

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

!$OMP PARALLEL DO PRIVATE(IEL,JEL,IS,XQ,XP,YQ,YP,ZQ,ZP,SIJ,DIJ,GRN,FLAG) SCHEDULE(DYNAMIC,16)

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