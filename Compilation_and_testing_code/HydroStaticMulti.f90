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
MODULE HydroStaticMulti

   USE HAMS_mod
   USE Body_mod
   USE Const_mod 
   USE WaveDyn_mod

   IMPLICIT NONE

   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: InitialisationMulti
   PUBLIC :: ReadHydroStaticMulti

CONTAINS
! ---------------------------------------------------------------------------------------------
!      Initialisation of the matrices, variables, and parameters
! ---------------------------------------------------------------------------------------------
! 
      SUBROUTINE InitialisationMulti
      IMPLICIT   NONE
      
      XW=0.D0
      AMP=1.D0
        
      XB_MULTI=0.D0
      MATX_MULTI=0.D0
      CRS_MULTI=0.D0
      BLNR_MULTI=0.D0
      BQDR_MULTI=0.D0
      KSTF_MULTI=0.D0
      
      EXFC_MULTI=CMPLX(0.D0,0.D0)
      AMAS_MULTI=0.D0
      BDMP_MULTI=0.D0
      DSPL_MULTI=0.D0
      
      RETURN
      END SUBROUTINE InitialisationMulti
      
      
! ---------------------------------------------------------------------------------------------
!      Read the data of gravity center, mass matrix and restoring matrix.
! ---------------------------------------------------------------------------------------------
! 
      SUBROUTINE ReadHydroStaticMulti(BODY_N)
      IMPLICIT   NONE  

      INTEGER I,J
      INTEGER,INTENT(IN):: BODY_N
      
      READ(4,*)
      READ(4,*)  (XG_MULTI(BODY_N,I), I=1,3)
      READ(4,*)
      DO I=1,6
       READ(4,120) (MATX_MULTI(BODY_N,I,J), J=1, 6)
      ENDDO
      READ(4,*)
      DO I=1,6
       READ(4,120) (BLNR_MULTI(BODY_N,I,J), J=1, 6)
      ENDDO
      READ(4,*)
      DO I=1,6
       READ(4,120) (BQDR_MULTI(BODY_N,I,J), J=1, 6)
      ENDDO
      READ(4,*)
      DO I=1,6
       READ(4,120) (CRS_MULTI(BODY_N,I,J), J=1, 6)
      ENDDO
      READ(4,*)
      DO I=1,6
       READ(4,120) (KSTF_MULTI(BODY_N,I,J), J=1, 6)
      ENDDO
      
      DO I=1,6
       DO J=1,6
       WRITE(65,130) I,J,CRS_MULTI(BODY_N,I,J)/(RHO*G)
       ENDDO
      ENDDO

120   FORMAT(6(2x,E12.5))
130   FORMAT(2I6,2X,ES14.6)

      RETURN        
      END SUBROUTINE ReadHydroStaticMulti
!-------------------------------------------------------------------------------
END MODULE HydroStaticMulti
!*******************************************************************************
