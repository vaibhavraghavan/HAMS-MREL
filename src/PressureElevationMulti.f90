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

MODULE PressureElevationMulti
    
   USE HAMS_mod
   USE Body_mod
   USE Const_mod
   USE WaveDyn_mod
   USE PanelMesh_mod
   USE Inerfs_mod
   USE Potentials_mod

   USE INFG3D_Open
   USE FinGreen3D_Open
   USE ImplementSubs
   USE SingularIntgr
   USE SingularIntgrMulti
   USE PatcVelct
   USE FieldOutput_mod

   IMPLICIT NONE

      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: CalPotentialMulti
   PUBLIC :: CalPressureMulti
   PUBLIC :: CalElevationMulti
   PUBLIC :: CalPotentialIncMulti
   
CONTAINS
! ----------------------------------------------------------------------------------------------------------------
!       Compute the velocity potential due to the separate radiation or diffraction modes for each multi-body mode
! ----------------------------------------------------------------------------------------------------------------
      SUBROUTINE CalPotentialMulti(XET,RDFLG,MD,POT,NBODY,PHIDOF,FACTOR_NORMAL)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: MD,NBODY,PHIDOF,FACTOR_NORMAL(:)
      REAL*8,INTENT(IN):: XET(3)
      COMPLEX*16,INTENT(OUT):: POT
      CHARACTER(*),INTENT(IN)::RDFLG

      INTEGER  JEL,IS,IRR,FLAG
      REAL*8   XQ(3),XP(3),XT(3),DIST,EAR,RKN(4),ENV(3),ENT(3),SLD
      COMPLEX*16  F0,DPOX,DPOY,DPOZ,DINCP,TERM1,TERM2,GRN(4),DUM(2)
      COMPLEX*16, ALLOCATABLE:: XPOT(:)
      
      ALLOCATE(XPOT(NELEM_PE))
        
      IRR=1
      XPOT=DCMPLX(0.0D0, 0.0D0)

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(JEL,XQ,ENV,EAR,IS,XP,DIST,FLAG,RKN,GRN,DUM,XT,ENT,DPOX,DPOY,DPOZ,DINCP,TERM1,TERM2) !$OMP REDUCTION(+:XPOT)
      
      DO JEL=1, NELEM_PE

       XQ(1)=XYZ_GLOBAL_MULTI_COMB_P(JEL,1)     ! XQ: source point,  XP: field point
       XQ(2)=XYZ_GLOBAL_MULTI_COMB_P(JEL,2)
       XQ(3)=XYZ_GLOBAL_MULTI_COMB_P(JEL,3)
       ENV(1)=DXYZ_MULTI_COMB(JEL,1)
       ENV(2)=DXYZ_MULTI_COMB(JEL,2)
       ENV(3)=DXYZ_MULTI_COMB(JEL,3)
       EAR=DS_MULTI_COMB(JEL)

       DO IS=1, NSYS

        IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
         XP(1)=SY_MULTI(IS,1)*XET(1)
         XP(2)=SX_MULTI(IS,1)*XET(2)
         XP(3)=         XET(3)
        ELSE
         XP(1)=SX_MULTI(IS,1)*XET(1)
         XP(2)=SY_MULTI(IS,1)*XET(2)
         XP(3)=         XET(3)
        ENDIF

        DIST=SQRT((XP(1)-XQ(1))**2+(XP(2)-XQ(2))**2+(XP(3)-XQ(3))**2)
        IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF

         IF (NCN_MULTI_COMB(JEL).EQ.3) THEN
          CALL SGLINTBD_MULTI_TRI2(XP,JEL,RKN,IRR)                                 ! Again obtaining the singularity term
         ELSEIF (NCN_MULTI_COMB(JEL).EQ.4) THEN
          CALL SGLINTBD_MULTI_QUAD2(XP,JEL,RKN,IRR)
         ENDIF

         IF (H.LT.0.D0) THEN
          CALL INFGREEN3D(XQ(1),XP(1),XQ(2),XP(2),XQ(3),XP(3),V,GRN,FLAG)    ! Obtaining the local term and oscillatory term of the Green's function
         ELSE
          CALL FINGREEN3D(XQ(1),XP(1),XQ(2),XP(2),XQ(3),XP(3),V,WVN,NK,H,GRN,FLAG)
         ENDIF

         IF (FLAG.EQ.1) THEN
          DUM(1)= RKN(1)+GRN(1)*EAR
          DUM(2)=(RKN(2)+GRN(2)*EAR)*ENV(1)   &
                +(RKN(3)+GRN(3)*EAR)*ENV(2)   &
              +(RKN(4)+GRN(4)*EAR)*ENV(3)
         ELSE
          DUM(1)= GRN(1)*EAR
          DUM(2)=(GRN(2)*ENV(1)+GRN(3)*ENV(2)+GRN(4)*ENV(3))*EAR
         ENDIF
         ! For multi-bodies, diffraction would be due to all elements, while the radiation would be due to only the elements of the body in the specific degree of freedom? 
         ! Yes, because diffraction field should be just all elements (it is just one mode)
         
         IF (MD.EQ.(6*NBODY+1)) THEN  
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN   
           XT(1)=SY_MULTI(1,IS)*XQ(1)
           XT(2)=SX_MULTI(1,IS)*XQ(2)
           XT(3)=         XQ(3)
           ENT(1)=SY_MULTI(1,IS)*ENV(1)
           ENT(2)=SX_MULTI(1,IS)*ENV(2)
           ENT(3)=         ENV(3)
          ELSE
           XT(1)=SX_MULTI(1,IS)*XQ(1)
           XT(2)=SY_MULTI(1,IS)*XQ(2)
           XT(3)=         XQ(3)
           ENT(1)=SX_MULTI(1,IS)*ENV(1)
           ENT(2)=SY_MULTI(1,IS)*ENV(2)
           ENT(3)=         ENV(3)
          ENDIF
          CALL DINP(XT(1),XT(2),XT(3),XW(1),XW(2),BETA,DPOX,DPOY,DPOZ)   ! Calculate the diffraction potential
          DINCP= DPOX*ENT(1)+DPOY*ENT(2)+DPOZ*ENT(3)
          TERM1=-DUM(1)*DINCP
          TERM2= DUM(2)*MXPOT_MULTI_COMB(JEL,MD,IS)
         ELSEIF (MD.EQ.1.or.MD.EQ.5) THEN
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           TERM1=DUM(1)*SY_MULTI(1,IS)*DXYZ_MULTI_COMB(JEL,MD)*FACTOR_NORMAL(JEL) !DUM(1) stays only when the normal vector is non-zero, which is for the considered DOF motion. For the rest of the elements, the MXPOT should only remain. For example for DOF1 of body 1 which is the first column of the MXPOT_MULTI, the normal vector is non zero and thus DUM(1) and DUM(2) are non-zero, while for Body2 and Body3 only the DUM(2) is non-zero. 
          ELSE
           TERM1=DUM(1)*SX_MULTI(1,IS)*DXYZ_MULTI_COMB(JEL,MD)*FACTOR_NORMAL(JEL)
          ENDIF
          TERM2=DUM(2)*MXPOT_MULTI_COMB(JEL,6*(PHIDOF-1)+MD,IS)
         ELSEIF (MD.EQ.2.or.MD.EQ.4) THEN
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           TERM1=DUM(1)*SX_MULTI(1,IS)*DXYZ_MULTI_COMB(JEL,MD)*FACTOR_NORMAL(JEL)
          ELSE
           TERM1=DUM(1)*SY_MULTI(1,IS)*DXYZ_MULTI_COMB(JEL,MD)*FACTOR_NORMAL(JEL)
          ENDIF
          TERM2=DUM(2)*MXPOT_MULTI_COMB(JEL,6*(PHIDOF-1)+MD,IS)
         ELSEIF (MD.EQ.3) THEN
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           TERM1=DUM(1)         *DXYZ_MULTI_COMB(JEL,MD)*FACTOR_NORMAL(JEL)
          ELSE
           TERM1=DUM(1)         *DXYZ_MULTI_COMB(JEL,MD)*FACTOR_NORMAL(JEL)
          ENDIF
          TERM2=DUM(2)*MXPOT_MULTI_COMB(JEL,6*(PHIDOF-1)+MD,IS)
         ELSEIF (MD.EQ.6) THEN
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           TERM1=DUM(1)*SX_MULTI(1,IS)*SY_MULTI(1,IS)*DXYZ_MULTI_COMB(JEL,MD)*FACTOR_NORMAL(JEL)
          ELSE
           TERM1=DUM(1)*SX_MULTI(1,IS)*SY_MULTI(1,IS)*DXYZ_MULTI_COMB(JEL,MD)*FACTOR_NORMAL(JEL)
          ENDIF
          TERM2=DUM(2)*MXPOT_MULTI_COMB(JEL,6*(PHIDOF-1)+MD,IS)
         ENDIF
        

         IF (ISOL.EQ.1) THEN                               
          XPOT(JEL)=XPOT(JEL)+TERM1-TERM2
         ELSEIF (ISOL.EQ.2) THEN
          XPOT(JEL)=XPOT(JEL)-TERM2 
         ENDIF
        
       ENDDO
      ENDDO

!$OMP END DO NOWAIT
!$OMP END PARALLEL
      
      POT=0.D0
      
      DO JEL=1,NELEM_PE
       POT=POT+XPOT(JEL)
      ENDDO


      SLD=4.D0*PI                                                                     
      POT=POT/SLD
      
      IF (MD.EQ.(6*NBODY+1)) THEN
        IF (ABS(TP+1.D0).LT.1.E-6.OR.ABS(TP).LT.1.E-6) THEN                          ! This is for 0 or infinity wave frequency
         F0=0.D0
        ELSE
         F0=VINP(XET(1),XET(2),XET(3),XW(1),XW(2),BETA)
        ENDIF
        POT=POT+F0
      ENDIF

      DEALLOCATE(XPOT)
      
      RETURN
      END SUBROUTINE CalPotentialMulti
      
! ----------------------------------------------------------------------------------------------------------------
!       Compute the velocity potential due to the incident wave
! ----------------------------------------------------------------------------------------------------------------
      SUBROUTINE CalPotentialIncMulti(XET,MD,POT,NBODY)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: MD,NBODY
      REAL*8,INTENT(IN):: XET(3)
      COMPLEX*16,INTENT(OUT):: POT
      COMPLEX*16  F0
      
      POT=0.D0
      
      IF (MD.EQ.(6*NBODY+1)) THEN
        IF (ABS(TP+1.D0).LT.1.E-6.OR.ABS(TP).LT.1.E-6) THEN                          ! This is for 0 or infinity wave frequency
         F0=0.D0
        ELSE
         F0=VINP(XET(1),XET(2),XET(3),XW(1),XW(2),BETA)
        ENDIF
        POT=POT+F0
      ENDIF
      
      RETURN
      END SUBROUTINE CalPotentialIncMulti

! -----------------------------------------------------------------------------------------------------
!       Compute the pressure at a point due to the separate radiation or diffraction modes
! ----------------------------------------------------------------------------------------------------- 
      SUBROUTINE CalPressureMulti(XP,RDFLG,MD,PRS,NBODY,PHIDOF,FACTOR_NORMAL)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: MD,NBODY,PHIDOF,FACTOR_NORMAL(:)
      REAL*8,INTENT(IN):: XP(3)
      COMPLEX*16,INTENT(OUT):: PRS
      CHARACTER(*),INTENT(IN)::RDFLG
      
      COMPLEX*16  XPOT
      
      CALL CalPotentialMulti(XP,RDFLG,MD,XPOT,NBODY,PHIDOF,FACTOR_NORMAL)
      
      IF (ABS(TP+1.D0).LT.1.E-6.OR.ABS(TP).LT.1.E-6) THEN                            ! This is the condition for 0 and infinite frequency
       PRS=RHO*XPOT
      ELSE
       PRS=CI*W1*RHO*XPOT                                                            ! Isn't this for the diffraction force?
      ENDIF

      RETURN
      END SUBROUTINE CalPressureMulti

! -----------------------------------------------------------------------------------------------------
!    Compute the free-surface elevation at a point due to the separate radiation or diffraction modes                              ! Can be directly implemented for computing RAOs
! ----------------------------------------------------------------------------------------------------- 
      SUBROUTINE CalElevationMulti(XP,RDFLG,MD,ELV,NBODY,PHIDOF,FACTOR_NORMAL)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: MD,NBODY,PHIDOF,FACTOR_NORMAL(:)
      REAL*8,INTENT(IN):: XP(3)
      COMPLEX*16,INTENT(OUT):: ELV
      CHARACTER(*),INTENT(IN)::RDFLG
      
      COMPLEX*16  XPOT
      
      CALL CalPotentialMulti(XP,RDFLG,MD,XPOT,NBODY,PHIDOF,FACTOR_NORMAL)
      
      IF (ABS(TP+1.D0).LT.1.E-6.OR.ABS(TP).LT.1.E-6) THEN
       ELV=XPOT
      ELSE
       ELV=CI*W1/G*XPOT
      ENDIF

      RETURN
      END SUBROUTINE CalElevationMulti
      
! -----------------------------------------------------------------------------------------------------
!    Compute the free-surface elevation at a point due to the separate radiation or diffraction modes                              ! This is the implementation in a dimensionless scheme in WAMIT format
! ----------------------------------------------------------------------------------------------------- 
      SUBROUTINE WamitNondimensMulti(VCP,PEFLG,RDFLG,MD,NVCP)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: MD
      INTEGER MEXP
      
      COMPLEX*16,INTENT(IN):: VCP
      COMPLEX*16,INTENT(OUT):: NVCP
      CHARACTER(*),INTENT(IN)::PEFLG,RDFLG
      REAL*8 NFAC

      IF (adjustl(trim(PEFLG)).EQ.'Pressure') THEN
       IF (ABS(TP+1.D0).LT.1.E-6.OR.ABS(TP).LT.1.E-6) THEN
        NFAC=RHO*AMP
       ELSE
        NFAC=RHO*G*AMP
       ENDIF
      ELSEIF (adjustl(trim(PEFLG)).EQ.'Elevation') THEN
       NFAC=AMP/W1**2
      ENDIF

      IF (adjustl(trim(RDFLG)).EQ.'Diffraction') THEN
       NVCP=VCP/AMP
      ELSEIF (adjustl(trim(RDFLG)).EQ.'Radiation') THEN
       IF (MD.LE.3) THEN
        MEXP=0
       ELSEIF (MD.GE.4) THEN
        MEXP=1
       ENDIF
       IF (ABS(TP+1.D0).LT.1.E-6.OR.ABS(TP).LT.1.E-6) THEN
        NVCP=VCP/NFAC*AMP/REFL**(MEXP+1)
       ELSE
        NVCP=VCP/NFAC*(-CI/W1)*AMP/REFL**MEXP
        !PRINT*, NFAC,REFL,MEXP
        !PAUSE
       ENDIF
      ENDIF

      IF (ABS(NVCP).LT.1.E-15) NVCP=0.D0
      
      IF (adjustl(trim(RDFLG)).EQ.'Diffraction') THEN
       NVCP=CMPLX(-IMAG(NVCP),-REAL(NVCP))
      ELSEIF (adjustl(trim(RDFLG)).EQ.'Radiation') THEN
       NVCP=CMPLX(REAL(NVCP),-IMAG(NVCP))
      ENDIF
          
      !NVCP=CONJG (NVCP)     ! Because of the Wamit format

      RETURN
      END SUBROUTINE WamitNondimensMulti
      
! -----------------------------------------------------------
!    Output pressures and elevations into the WAMIT format for radiation
! -----------------------------------------------------------
      SUBROUTINE OutputPressureElevation_RadiationMulti(NFILE,NBODY)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: NFILE,NBODY
      INTEGER IPT,MD,ND,EMD,IHD,FILE_NUMBER,NELEM_START,NELEM_END,NELEM_GLOBAL
      REAL*8 XP(3)
      COMPLEX*16 VCP
      COMPLEX*16,ALLOCATABLE:: VCPX(:,:)
      INTEGER,ALLOCATABLE:: FACTOR_NORMAL(:)
      
      ALLOCATE(FACTOR_NORMAL(NELEM_PE))
      NELEM_GLOBAL=0
      
      DO FILE_NUMBER = 1,NBODY
        ALLOCATE(VCPX(NFP,6))
        FACTOR_NORMAL = 0
        IF (FILE_NUMBER.EQ.1) THEN
         NELEM_START = 1
         NELEM_END = NELEM_MULTI(FILE_NUMBER)
         FACTOR_NORMAL(NELEM_START:NELEM_END) = 1
        ELSEIF (FILE_NUMBER.EQ.NBODY) THEN
         NELEM_START = NELEM_PE-NELEM_MULTI(FILE_NUMBER)+1
         NELEM_END = NELEM_PE
         FACTOR_NORMAL(NELEM_START:NELEM_END) = 1
        ELSE
         NELEM_START = NELEM_GLOBAL+1
         NELEM_END = NELEM_GLOBAL+NELEM_MULTI(FILE_NUMBER)
         FACTOR_NORMAL(NELEM_START:NELEM_END) = 1
        ENDIF
        NELEM_GLOBAL=NELEM_GLOBAL+NELEM_MULTI(FILE_NUMBER)
    
       DO IPT=1,NFP
        XP=XFP(IPT,:)
          DO MD=1,6
           IF (ABS(XP(3)).GT.1.E-6) THEN
            CALL CalPressureMulti(XP,'Radiation',MD,VCP,NBODY,FILE_NUMBER,FACTOR_NORMAL)
            CALL WamitNondimensMulti(VCP,'Pressure','Radiation',MD,VCPX(IPT,MD)) ! Check this after fixing the potential function
           ELSE
            CALL CalElevationMulti(XP,'Radiation',MD,VCP,NBODY,FILE_NUMBER,FACTOR_NORMAL)
            CALL WamitNondimensMulti(VCP,'Elevation','Radiation',MD,VCPX(IPT,MD))
           ENDIF
      
!   ===================================================
!    Write WAMIT-style output files (Both pressure and displacement fields)
!   
      
          ENDDO
          WRITE(NFILE+FILE_NUMBER,1000) OUFR,IPT,(VCPX(IPT,ND),ND=1,6)
      
       ENDDO
       DEALLOCATE(VCPX)
      ENDDO
      
1000  FORMAT(ES14.6,I10,12ES14.6)
      RETURN
      END SUBROUTINE OutputPressureElevation_RadiationMulti
      
! -----------------------------------------------------------
!    Output pressures and elevations into the WAMIT format for diffraction
! -----------------------------------------------------------
      SUBROUTINE OutputPressureElevation_DiffractionMulti(NFILE,NBODY)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: NFILE,NBODY
      INTEGER IPT,MD,EMD,IHD
      REAL*8 XP(3),REL,IMG,MOD,PHS
      COMPLEX*16 VCP,NVCP
      INTEGER,ALLOCATABLE:: FACTOR_NORMAL(:)
      
      ALLOCATE(FACTOR_NORMAL(NELEM_PE))
      FACTOR_NORMAL = 0

      DO IPT=1,NFP
       XP=XFP(IPT,:)
       IF (ABS(XP(3)).GT.1.E-6) THEN
        CALL CalPressureMulti(XP,'Diffraction',6*NBODY+1,VCP,NBODY,0,FACTOR_NORMAL)
        CALL WamitNondimensMulti(VCP,'Pressure','Diffraction',0,NVCP)
       ELSE
        CALL CalElevationMulti(XP,'Diffraction',6*NBODY+1,VCP,NBODY,0,FACTOR_NORMAL)
        CALL WamitNondimensMulti(VCP,'Elevation','Diffraction',0,NVCP)
       ENDIF
       
       !WRITE(NFILE,1020) OUFR,BETA*180.0D0/PI,IPT,NVCP
       REL=REAL(NVCP)
       IMG=IMAG(NVCP)
       MOD=ABS(NVCP)
       PHS=ATAN2D(IMG,REL)
       WRITE(NFILE,1020) OUFR,BETA*180.0D0/PI,IPT,MOD,PHS,REL,IMG
       
      ENDDO
      
1020  FORMAT(2ES14.6,I10,4ES14.6)
      RETURN
      END SUBROUTINE OutputPressureElevation_DiffractionMulti
      
      ! -----------------------------------------------------------------------------------------------------
      !       Compute the pressure at a point due to the incidence wave
      ! ----------------------------------------------------------------------------------------------------- 
      SUBROUTINE CalPressureIncMulti(XP,MD,PRS,NBODY)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: MD,NBODY
      REAL*8,INTENT(IN):: XP(3)
      COMPLEX*16,INTENT(OUT):: PRS
      
      COMPLEX*16  XPOT
      
      CALL CalPotentialIncMulti(XP,MD,XPOT,NBODY)
      
      IF (ABS(TP+1.D0).LT.1.E-6.OR.ABS(TP).LT.1.E-6) THEN                            ! This is the condition for 0 and infinite frequency
       PRS=RHO*XPOT
      ELSE
       PRS=CI*W1*RHO*XPOT                                                            ! Isn't this for the diffraction force?
      ENDIF

      RETURN
      END SUBROUTINE CalPressureIncMulti

! -----------------------------------------------------------------------------------------------------
!    Compute the free-surface elevation at a point due to the incidence wave                              ! Can be directly implemented for computing RAOs
! ----------------------------------------------------------------------------------------------------- 
      SUBROUTINE CalElevationIncMulti(XP,MD,ELV,NBODY)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: MD,NBODY
      REAL*8,INTENT(IN):: XP(3)
      COMPLEX*16,INTENT(OUT):: ELV
      
      COMPLEX*16  XPOT
      
      CALL CalPotentialIncMulti(XP,MD,XPOT,NBODY)
      
      IF (ABS(TP+1.D0).LT.1.E-6.OR.ABS(TP).LT.1.E-6) THEN
       ELV=XPOT
      ELSE
       ELV=CI*W1/G*XPOT
      ENDIF

      RETURN
      END SUBROUTINE CalElevationIncMulti
      
      ! -----------------------------------------------------------
      !    Output pressures and elevations into the WAMIT format for Incident wave
      ! -----------------------------------------------------------
      SUBROUTINE OutputPressureElevation_IncidenceMulti(NFILE,NBODY)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: NFILE,NBODY
      INTEGER IPT,MD,EMD,IHD
      REAL*8 XP(3),REL,IMG,MOD,PHS
      COMPLEX*16 VCP,NVCP

      DO IPT=1,NFP
       XP=XFP(IPT,:)
       IF (ABS(XP(3)).GT.1.E-6) THEN
        CALL CalPressureIncMulti(XP,6*NBODY+1,VCP,NBODY)
        CALL WamitNondimensMulti(VCP,'Pressure','Diffraction',0,NVCP)
       ELSE
        CALL CalElevationIncMulti(XP,6*NBODY+1,VCP,NBODY)
        CALL WamitNondimensMulti(VCP,'Elevation','Diffraction',0,NVCP)
       ENDIF
       
       !WRITE(NFILE,1020) OUFR,BETA*180.0D0/PI,IPT,NVCP
       REL=REAL(NVCP)
       IMG=IMAG(NVCP)
       MOD=ABS(NVCP)
       PHS=ATAN2D(IMG,REL)
       WRITE(NFILE,1020) OUFR,BETA*180.0D0/PI,IPT,MOD,PHS,REL,IMG
       
      ENDDO
      
1020  FORMAT(2ES14.6,I10,4ES14.6)
      RETURN
      END SUBROUTINE OutputPressureElevation_IncidenceMulti
!-------------------------------------------------------------------------------
END MODULE PressureElevationMulti
!*******************************************************************************
