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
MODULE BodyIntgr_irrMulti

   USE HAMS_mod
   USE Body_mod
   USE Const_mod
   USE WaveDyn_mod
   USE PanelMesh_mod
   USE Inerfs_mod
   USE PatcVelct
   USE Potentials_mod
   USE omp_lib
   IMPLICIT NONE
   
      ! ..... Public Subroutines ...................................................................................................
   
   PUBLIC :: BODINT_LEFT_IRR_MULTI
   PUBLIC :: RBC_IRR_MULTI
   PUBLIC :: DBC_IRR_MULTI
   
CONTAINS


!   --------------------------------------------------------
!     Integration inside an element, which is not singular 
!   --------------------------------------------------------

      SUBROUTINE BODINT_LEFT_IRR_MULTI(IS,IEL,JEL,TINDP,IRR,FLAG)
      IMPLICIT NONE
      
      INTEGER,INTENT(IN):: IS,IEL,JEL,IRR,FLAG
      COMPLEX*16,INTENT(OUT):: TINDP(4)
      
      INTEGER  IP  
      COMPLEX*16  GRN(4)
      REAL*8::   RKN(4),ENV(6),TNV(6),EAR,XEC(3)
      
      IF (IRR.EQ.1) THEN
          
        RKN(:)=RKBN_MULTI_COMB(IEL,JEL,IS,:)
        GRN(:)=CGRN_MULTI_COMB(IEL,JEL,IS,:)
        
        EAR=DS_MULTI_COMB(JEL)
        XEC(:)=XYZ_GLOBAL_MULTI_COMB_P(JEL,:)
        ENV(:)=DXYZ_MULTI_COMB(JEL,:)
        TNV(:)=DXYZ_MULTI_COMB(IEL,:)
        
      ELSEIF (IRR.EQ.3) THEN
          
        RKN(:)=PKBN_MULTI_COMB(IEL,JEL,IS,:)
        GRN(:)=DGRN_MULTI_COMB(IEL,JEL,IS,:)
        
        EAR=DS_MULTI_COMB(JEL)
        XEC(:)=XYZ_GLOBAL_MULTI_COMB_P(JEL,:)
        ENV(:)=DXYZ_MULTI_COMB(JEL,:)
        TNV(:)=iDXYZ_MULTI_COMB(IEL,:)
        
      ENDIF

      IF (FLAG.EQ.1) THEN

        TINDP(IS)=(RKN(2)+GRN(2)*EAR)*ENV(1)   &
                 +(RKN(3)+GRN(3)*EAR)*ENV(2)   &
                 +(RKN(4)+GRN(4)*EAR)*ENV(3)

      ELSE

        TINDP(IS)=(GRN(2)*ENV(1)+GRN(3)*ENV(2)+GRN(4)*ENV(3))*EAR
        
      ENDIF
  
      RETURN
      END SUBROUTINE BODINT_LEFT_IRR_MULTI

!   --------------------------------------------------------
!     Assemble the radiated right-side boundary conditions
!   --------------------------------------------------------

      SUBROUTINE RBC_IRR_MULTI(IS,IEL,JEL,TINRD,IRR,FLAG)
      IMPLICIT NONE
      
      INTEGER,INTENT(IN):: IS,IEL,JEL,IRR,FLAG
      COMPLEX*16,INTENT(OUT):: TINRD(4,6,4)
      
      INTEGER  IP
      REAL*8  XQ,YQ,ZQ,EAR
      REAL*8  RKN(4),ENV(6),TNV(6),XEC(3)
      COMPLEX*16  GRN(4),DUM

      IF (IRR.EQ.1) THEN
          
        RKN(:)=RKBN_MULTI_COMB(IEL,JEL,IS,:)
        GRN(:)=CGRN_MULTI_COMB(IEL,JEL,IS,:)
        
        EAR=DS_MULTI_COMB(JEL)
        XEC(:)=XYZ_GLOBAL_MULTI_COMB_P(JEL,:)
        ENV(:)=DXYZ_MULTI_COMB(JEL,:)
        TNV(:)=DXYZ_MULTI_COMB(IEL,:)
        
      ELSEIF (IRR.EQ.3) THEN
          
        RKN(:)=PKBN_MULTI_COMB(IEL,JEL,IS,:)
        GRN(:)=DGRN_MULTI_COMB(IEL,JEL,IS,:)
        
        EAR=DS_MULTI_COMB(JEL)
        XEC(:)=XYZ_GLOBAL_MULTI_COMB_P(JEL,:)
        ENV(:)=DXYZ_MULTI_COMB(JEL,:)
        TNV(:)=iDXYZ_MULTI_COMB(IEL,:)
        
      ENDIF
    
!  Above----
!  Field point changes, source element keeps the same as in the first quadrant

      IF (FLAG.EQ.1) THEN
        DUM=RKN(1)+GRN(1)*EAR
      ELSE
        DUM=GRN(1)*EAR
      ENDIF

      DO 200 IP=1, NSYS
       IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
        TINRD(IS,1,IP)=TINRD(IS,1,IP)+DUM*SY_MULTI(IS,IP)*ENV(1)
        TINRD(IS,2,IP)=TINRD(IS,2,IP)+DUM*SX_MULTI(IS,IP)*ENV(2)
        TINRD(IS,3,IP)=TINRD(IS,3,IP)+DUM*ENV(3)
        TINRD(IS,4,IP)=TINRD(IS,4,IP)+DUM*SX_MULTI(IS,IP)*ENV(4)
        TINRD(IS,5,IP)=TINRD(IS,5,IP)+DUM*SY_MULTI(IS,IP)*ENV(5)
        TINRD(IS,6,IP)=TINRD(IS,6,IP)+DUM*SY_MULTI(IS,IP)*SX_MULTI(IS,IP)*ENV(6)
       ELSE
        TINRD(IS,1,IP)=TINRD(IS,1,IP)+DUM*SX_MULTI(IS,IP)*ENV(1)
        TINRD(IS,2,IP)=TINRD(IS,2,IP)+DUM*SY_MULTI(IS,IP)*ENV(2)
        TINRD(IS,3,IP)=TINRD(IS,3,IP)+DUM*ENV(3)
        TINRD(IS,4,IP)=TINRD(IS,4,IP)+DUM*SY_MULTI(IS,IP)*ENV(4)
        TINRD(IS,5,IP)=TINRD(IS,5,IP)+DUM*SX_MULTI(IS,IP)*ENV(5)
        TINRD(IS,6,IP)=TINRD(IS,6,IP)+DUM*SX_MULTI(IS,IP)*SY_MULTI(IS,IP)*ENV(6)
       ENDIF

200   CONTINUE
   
      RETURN
      END SUBROUTINE RBC_IRR_MULTI
      

!   --------------------------------------------------------
!     Assemble the diffracted right-side boundary conditions
!   --------------------------------------------------------

      SUBROUTINE DBC_IRR_MULTI(IS,IEL,JEL,TINRD,IRR,FLAG)
      IMPLICIT NONE 
      
      INTEGER,INTENT(IN):: IS,IEL,JEL,IRR,FLAG
      COMPLEX*16,INTENT(OUT):: TINRD(4,4)
      
      INTEGER  IP
      REAL*8  XQ,YQ,ZQ,EAR
      REAL*8  RKN(4),ENV(6),TNV(6),XEC(3)
      
      COMPLEX*16  GRN(4),DUM(2)
      COMPLEX*16  DPOX,DPOY,DPOZ

      IF (IRR.EQ.1) THEN      

        RKN(:)=RKBN_MULTI_COMB(IEL,JEL,IS,:)
        GRN(:)=CGRN_MULTI_COMB(IEL,JEL,IS,:)
        
        EAR=DS_MULTI_COMB(JEL)
        XEC(:)=XYZ_GLOBAL_MULTI_COMB_P(JEL,:)
        ENV(:)=DXYZ_MULTI_COMB(JEL,:)
        TNV(:)=DXYZ_MULTI_COMB(IEL,:)
        
      ELSEIF (IRR.EQ.3) THEN  
          
        RKN(:)=PKBN_MULTI_COMB(IEL,JEL,IS,:)
        GRN(:)=DGRN_MULTI_COMB(IEL,JEL,IS,:)
        
        EAR=DS_MULTI_COMB(JEL)
        XEC(:)=XYZ_GLOBAL_MULTI_COMB_P(JEL,:)
        ENV(:)=DXYZ_MULTI_COMB(JEL,:)
        TNV(:)=iDXYZ_MULTI_COMB(IEL,:)
        
      ENDIF    

!  Above----
!  Field point changes, source element keeps the same as in the first quadrant

      IF (FLAG.EQ.1) THEN
        DUM(1)=RKN(1)+GRN(1)*EAR
      ELSE
        DUM(1)=GRN(1)*EAR
      ENDIF

      DO 200 IP=1, NSYS

       IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
        XQ=SY_MULTI(IS,IP)*XEC(1)
        YQ=SX_MULTI(IS,IP)*XEC(2)
        ZQ=          XEC(3)
        CALL DINP(XQ,YQ,ZQ,XW(1),XW(2),BETA,DPOX,DPOY,DPOZ)
        DUM(2)=SY_MULTI(IS,IP)*DPOX*ENV(1)+SX_MULTI(IS,IP)*DPOY*ENV(2)+DPOZ*ENV(3)          
        TINRD(IS,IP)=TINRD(IS,IP)-DUM(1)*DUM(2)
       ELSE
        XQ=SX_MULTI(IS,IP)*XEC(1)
        YQ=SY_MULTI(IS,IP)*XEC(2)
        ZQ=          XEC(3)
        CALL DINP(XQ,YQ,ZQ,XW(1),XW(2),BETA,DPOX,DPOY,DPOZ)
        DUM(2)=SX_MULTI(IS,IP)*DPOX*ENV(1)+SY_MULTI(IS,IP)*DPOY*ENV(2)+DPOZ*ENV(3)
        TINRD(IS,IP)=TINRD(IS,IP)-DUM(1)*DUM(2)
       ENDIF

200   CONTINUE
      
      RETURN
      END SUBROUTINE DBC_IRR_MULTI
!-------------------------------------------------------------------------------
END MODULE BodyIntgr_irrMulti
!*******************************************************************************