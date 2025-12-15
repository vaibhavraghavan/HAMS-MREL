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
MODULE PotentWavForceMulti

   USE HAMS_mod
   USE Body_mod
   USE Const_mod
   USE PanelMesh_mod
   USE Inerfs_mod

   USE PatcVelct
   USE Potentials_mod
   IMPLICIT NONE

      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: EFORCE_MULTI
   PUBLIC :: RFORCE_MULTI
   PUBLIC :: ReorderAmssDamp
   
CONTAINS
! ------------------------------------------------------------------- 
!          Compute the wave exciting force on a 3D body
! ------------------------------------------------------------------- 
      SUBROUTINE EFORCE_MULTI(WK,W1,TP,BETA,AMP,EXFC,NELEM_START,NELEM_END,NBODY,DOF1)
      IMPLICIT NONE
!
      REAL*8,INTENT(IN):: WK,W1,TP,BETA,AMP
      COMPLEX*16,INTENT(OUT):: EXFC(6)
      INTEGER,INTENT(IN):: NELEM_START,NELEM_END,NBODY,DOF1
      
      INTEGER IEL,IP,MD
      REAL*8:: XP,YP,ZP,AMFJ(6)
      REAL*8:: MOD,PHS(6),REL,IMG,NREL,NIMG
      COMPLEX*16 PHI,FORCE(6,4)

      MD=6*NBODY+1
      FORCE=CMPLX(0.0D0, 0.0D0)
      
!!$CALL OMP_SET_NUM_THREADS(NTHREAD)
!!$OMP PARALLEL DO PRIVATE(IEL,IP,XP,YP,ZP,PHI) !$OMP REDUCTION(+:FORCE)

      DO 100  IEL=NELEM_START,  NELEM_END               ! Here the variation should also be for all elements, but the normal vector is only on the surface of the body considered.
      DO 100  IP=1,  NSYS
          
       IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
        XP=RX_MULTI(IP,1)*XYZ_GLOBAL_MULTI_COMB_P(IEL,1)
        YP=RX_MULTI(IP,2)*XYZ_GLOBAL_MULTI_COMB_P(IEL,2)
        ZP=         XYZ_GLOBAL_MULTI_COMB_P(IEL,3)
        IF (ISOL.EQ.1) THEN
         PHI=(MXPOT_MULTI_COMB(IEL,MD,IP)+VINP(XP,YP,ZP,XW(1),XW(2),BETA))*DS_MULTI_COMB(IEL)
        ELSEIF (ISOL.EQ.2) THEN
         PHI=MXPOT_MULTI_COMB(IEL,MD,IP)*DS_MULTI_COMB(IEL)
        ENDIF
        FORCE(1,IP)=FORCE(1,IP)+PHI*RX_MULTI(IP,1)*DXYZ_MULTI_COMB(IEL,1)
        FORCE(2,IP)=FORCE(2,IP)+PHI*RX_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,2)
        FORCE(3,IP)=FORCE(3,IP)+PHI*         DXYZ_MULTI_COMB(IEL,3)
        FORCE(4,IP)=FORCE(4,IP)+PHI*RX_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,4)
        FORCE(5,IP)=FORCE(5,IP)+PHI*RX_MULTI(IP,1)*DXYZ_MULTI_COMB(IEL,5)
        FORCE(6,IP)=FORCE(6,IP)+PHI*RX_MULTI(IP,1)*RX_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,6)
       ELSE
        XP=RY_MULTI(IP,1)*XYZ_GLOBAL_MULTI_COMB_P(IEL,1)
        YP=RY_MULTI(IP,2)*XYZ_GLOBAL_MULTI_COMB_P(IEL,2)
        ZP=         XYZ_GLOBAL_MULTI_COMB_P(IEL,3)
        IF (ISOL.EQ.1) THEN
         PHI=(MXPOT_MULTI_COMB(IEL,MD,IP)+VINP(XP,YP,ZP,XW(1),XW(2),BETA))*DS_MULTI_COMB(IEL)
        ELSEIF (ISOL.EQ.2) THEN
         PHI=MXPOT_MULTI_COMB(IEL,MD,IP)*DS_MULTI_COMB(IEL)
        ENDIF
        FORCE(1,IP)=FORCE(1,IP)+PHI*RY_MULTI(IP,1)*DXYZ_MULTI_COMB(IEL,1)
        FORCE(2,IP)=FORCE(2,IP)+PHI*RY_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,2)
        FORCE(3,IP)=FORCE(3,IP)+PHI*          DXYZ_MULTI_COMB(IEL,3)
        FORCE(4,IP)=FORCE(4,IP)+PHI*RY_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,4)
        FORCE(5,IP)=FORCE(5,IP)+PHI*RY_MULTI(IP,1)*DXYZ_MULTI_COMB(IEL,5)
        FORCE(6,IP)=FORCE(6,IP)+PHI*RY_MULTI(IP,1)*RY_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,6)
       ENDIF

100   CONTINUE
!!$OMP END PARALLEL DO
! -------------------------------------------------
!
      DO 200 MD=1, 6
        EXFC(MD)=CMPLX(0.D0,0.D0)
      DO  200  IP=1,  NSYS
	    EXFC(MD)=EXFC(MD)+FORCE(MD,IP)
200   CONTINUE

      EXFC(:)=CI*W1*RHO*EXFC(:)
      AMFJ(:)=ABS(EXFC(:))

      DO MD=1,6
          
       REL=REAL(EXFC(MD))/(RHO*G*AMP)
       IMG=IMAG(EXFC(MD))/(RHO*G*AMP)
       MOD=ABS(EXFC(MD))/(RHO*G*AMP) !SQRT(REL**2+MDMG**2)
       NREL=-IMG
       NIMG=-REL
       PHS(MD)=ATAN2D(NIMG,NREL)*180.0D0/PI
       WRITE(20+MD,1010) WK,W1,REAL(EXFC(MD)),IMAG(EXFC(MD)) ! These are different from the values of WAMIT? These change the RAOs. Why is that happening?
       
       IF (ABS(TP+1.D0).GT.1.E-6.AND.ABS(TP).GT.1.E-6) THEN
        WRITE(62,1030)  OUFR, BETA*180.0D0/PI, (DOF1 - 1)*6 + MD, MOD, PHS(MD), NREL, NIMG
       ENDIF
           
      ENDDO
!
!   =================================================== 
!
1010  FORMAT(F7.3,1x,F7.3,1x,6E14.6)
1030  FORMAT(2ES14.6,I6,4ES14.6)
      
      RETURN 
      END SUBROUTINE EFORCE_MULTI 

! ------------------------------------------------------------------- 
!          Compute the wave radiation force on a 3D body
! ------------------------------------------------------------------- 
      SUBROUTINE RFORCE_MULTI(WK,W1,TP,AMAS,BDMP,NELEM_START,NELEM_END,PHIDOF,DOF1,DOF2)
      IMPLICIT   NONE
!
      REAL*8,INTENT(IN):: WK,W1,TP
      REAL*8,INTENT(OUT):: AMAS(6,6),BDMP(6,6)
      INTEGER,INTENT(IN):: NELEM_START,NELEM_END,PHIDOF,DOF1,DOF2
      
      INTEGER IEL,JEL,MD,MD1,MD2,IP
      REAL*8:: AMFJ(6),NAMAS(6,6),NBDMP(6,6)
      COMPLEX*16 RPHI,IPHI
      
      ! Initialize arrays to zero to prevent NaN from uninitialized memory
      AMAS = 0.0D0
      BDMP = 0.0D0
      
      !!$CALL OMP_SET_NUM_THREADS(NTHREAD)
      !!$OMP PARALLEL DO PRIVATE(IEL,IP,MD,RPHI,IPHI) !$OMP REDUCTION(+:AMAS,BDMP)

      DO 100  MD=1,6
        DO 100  IEL=NELEM_START,  NELEM_END
        DO 100  IP=1,  NSYS

         IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
          RPHI=CI*RHO*MXPOT_MULTI_COMB(IEL,6*(PHIDOF-1)+MD,IP)*DS_MULTI_COMB(IEL)
          IPHI=CI*RHO*W1*MXPOT_MULTI_COMB(IEL,6*(PHIDOF-1)+MD,IP)*DS_MULTI_COMB(IEL)

          AMAS(1,MD)=AMAS(1,MD)+IMAG(RPHI*RX_MULTI(IP,1)*DXYZ_MULTI_COMB(IEL,1))
          AMAS(2,MD)=AMAS(2,MD)+IMAG(RPHI*RX_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,2))
          AMAS(3,MD)=AMAS(3,MD)+IMAG(RPHI*           DXYZ_MULTI_COMB(IEL,3))
          AMAS(4,MD)=AMAS(4,MD)+IMAG(RPHI*RX_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,4))
          AMAS(5,MD)=AMAS(5,MD)+IMAG(RPHI*RX_MULTI(IP,1)*DXYZ_MULTI_COMB(IEL,5))
          AMAS(6,MD)=AMAS(6,MD)+IMAG(RPHI*RX_MULTI(IP,1)*RX_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,6))
         
          BDMP(1,MD)=BDMP(1,MD)-DBLE(IPHI*RX_MULTI(IP,1)*DXYZ_MULTI_COMB(IEL,1))
          BDMP(2,MD)=BDMP(2,MD)-DBLE(IPHI*RX_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,2))
          BDMP(3,MD)=BDMP(3,MD)-DBLE(IPHI*           DXYZ_MULTI_COMB(IEL,3))
          BDMP(4,MD)=BDMP(4,MD)-DBLE(IPHI*RX_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,4))
          BDMP(5,MD)=BDMP(5,MD)-DBLE(IPHI*RX_MULTI(IP,1)*DXYZ_MULTI_COMB(IEL,5))
          BDMP(6,MD)=BDMP(6,MD)-DBLE(IPHI*RX_MULTI(IP,1)*RX_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,6))
         ELSE
          
          RPHI=CI*RHO*MXPOT_MULTI_COMB(IEL,6*(PHIDOF-1)+MD,IP)*DS_MULTI_COMB(IEL)
          IPHI=CI*RHO*W1*MXPOT_MULTI_COMB(IEL,6*(PHIDOF-1)+MD,IP)*DS_MULTI_COMB(IEL)

          AMAS(1,MD)=AMAS(1,MD)+IMAG(RPHI*RY_MULTI(IP,1)*DXYZ_MULTI_COMB(IEL,1))
          AMAS(2,MD)=AMAS(2,MD)+IMAG(RPHI*RY_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,2))
          AMAS(3,MD)=AMAS(3,MD)+IMAG(RPHI*          DXYZ_MULTI_COMB(IEL,3))  
          AMAS(4,MD)=AMAS(4,MD)+IMAG(RPHI*RY_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,4))
          AMAS(5,MD)=AMAS(5,MD)+IMAG(RPHI*RY_MULTI(IP,1)*DXYZ_MULTI_COMB(IEL,5))
          AMAS(6,MD)=AMAS(6,MD)+IMAG(RPHI*RY_MULTI(IP,1)*RY_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,6)) 
         
          BDMP(1,MD)=BDMP(1,MD)-DBLE(IPHI*RY_MULTI(IP,1)*DXYZ_MULTI_COMB(IEL,1))
          BDMP(2,MD)=BDMP(2,MD)-DBLE(IPHI*RY_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,2))
          BDMP(3,MD)=BDMP(3,MD)-DBLE(IPHI*          DXYZ_MULTI_COMB(IEL,3))  
          BDMP(4,MD)=BDMP(4,MD)-DBLE(IPHI*RY_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,4))
          BDMP(5,MD)=BDMP(5,MD)-DBLE(IPHI*RY_MULTI(IP,1)*DXYZ_MULTI_COMB(IEL,5))
          BDMP(6,MD)=BDMP(6,MD)-DBLE(IPHI*RY_MULTI(IP,1)*RY_MULTI(IP,2)*DXYZ_MULTI_COMB(IEL,6)) 
         ENDIF

      100   CONTINUE
      
      !!$OMP END PARALLEL DO
      
       DO MD1=1,6
        WRITE(30+MD1,1000) WK,W1,(AMAS(MD1,MD2),MD2=1,6)
        WRITE(40+MD1,1000) WK,W1,(BDMP(MD1,MD2),MD2=1,6)
       ENDDO  
!
!   =================================================== 
!    Write WAMIT-style output files
!
       DO MD1=1,6
          DO MD2=1,6
            IF (ABS(TP+1.D0).LT.1.E-6.OR.ABS(TP).LT.1.E-6) THEN
             WRITE(61,1020) OUFR, (DOF1 - 1)*6 + MD1, (DOF2 - 1)*6 + MD2, AMAS(MD1,MD2)/RHO
            ELSE
             WRITE(61,1020) OUFR, (DOF1 - 1)*6 + MD1, (DOF2 - 1)*6 + MD2, AMAS(MD1,MD2)/RHO, BDMP(MD1,MD2)/(RHO*W1)
            ENDIF
          ENDDO
       ENDDO
!   =================================================== 
!    Write HydroStar-style output files
!

1000    FORMAT(F7.3,1x,F7.3,1x,6E14.5)
1020    FORMAT(ES14.6,2I6,2ES14.6)

      RETURN
      END SUBROUTINE RFORCE_MULTI
      
!  ------------------------------------------------------------------------------------------------------
!    Subroutine to reorder WAMIT-style hydrodynamic coefficient files
!    Reorders from body-pair block format to row-major format
!    (all columns for row 1, then all columns for row 2, etc.)
!    Cross-platform: works on both Windows and Linux
!  ------------------------------------------------------------------------------------------------------

      SUBROUTINE ReorderAmssDamp(filename, NFREQ, NDOF)
          IMPLICIT NONE
      
          CHARACTER(LEN=*), INTENT(IN) :: filename
          INTEGER, INTENT(IN) :: NFREQ    ! Number of frequencies
          INTEGER, INTENT(IN) :: NDOF     ! Total DOFs (e.g., 36 for 6 bodies x 6 DOF)
      
          ! Local variables
          INTEGER :: I, IOS, IFREQ, IROW, ICOL, LINE_COUNT, NFREQ_FOUND
          REAL*8, ALLOCATABLE :: FREQ(:,:,:)      ! (NFREQ, NDOF, NDOF)
          REAL*8, ALLOCATABLE :: AMAS(:,:,:)      ! (NFREQ, NDOF, NDOF)
          REAL*8, ALLOCATABLE :: BDMP(:,:,:)      ! (NFREQ, NDOF, NDOF)
          REAL*8, ALLOCATABLE :: FREQ_LIST(:)     ! List of frequencies
      
          REAL*8 :: FREQ_VAL, AMAS_VAL, BDMP_VAL, LAST_FREQ
          INTEGER :: ROW_VAL, COL_VAL
          LOGICAL :: file_exists, file_open
      
          ! Check if file exists
          INQUIRE(FILE=TRIM(filename), EXIST=file_exists)
          IF (.NOT. file_exists) RETURN
          
          ! Make sure unit 61 is closed (the original write unit for AmssDamp.1)
          INQUIRE(UNIT=61, OPENED=file_open)
          IF (file_open) CLOSE(61)
      
          ! Allocate arrays
          ALLOCATE(FREQ(NFREQ, NDOF, NDOF), STAT=IOS)
          IF (IOS /= 0) RETURN
          ALLOCATE(AMAS(NFREQ, NDOF, NDOF), STAT=IOS)
          IF (IOS /= 0) THEN
             DEALLOCATE(FREQ)
             RETURN
          ENDIF
          ALLOCATE(BDMP(NFREQ, NDOF, NDOF), STAT=IOS)
          IF (IOS /= 0) THEN
             DEALLOCATE(FREQ, AMAS)
             RETURN
          ENDIF
          ALLOCATE(FREQ_LIST(NFREQ), STAT=IOS)
          IF (IOS /= 0) THEN
             DEALLOCATE(FREQ, AMAS, BDMP)
             RETURN
          ENDIF
      
          ! Initialize
          FREQ = 0.0D0
          AMAS = 0.0D0
          BDMP = 0.0D0
          FREQ_LIST = 0.0D0
      
          ! Read the original file
          OPEN(UNIT=901, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', IOSTAT=IOS)
          IF (IOS /= 0) THEN
             DEALLOCATE(FREQ, AMAS, BDMP, FREQ_LIST)
             RETURN
          ENDIF
          
          ! Read all data, detecting frequency changes
          IFREQ = 0
          LAST_FREQ = -999.0D0
          
          DO
             READ(901, *, IOSTAT=IOS) FREQ_VAL, ROW_VAL, COL_VAL, AMAS_VAL, BDMP_VAL
             IF (IOS /= 0) EXIT
             
             ! Detect new frequency
             IF (ABS(FREQ_VAL - LAST_FREQ) > 1.0D-10) THEN
                IFREQ = IFREQ + 1
                IF (IFREQ > NFREQ) IFREQ = NFREQ
                FREQ_LIST(IFREQ) = FREQ_VAL
                LAST_FREQ = FREQ_VAL
             ENDIF
             
             ! Store data if indices are valid
             IF (ROW_VAL >= 1 .AND. ROW_VAL <= NDOF .AND. &
                 COL_VAL >= 1 .AND. COL_VAL <= NDOF .AND. &
                 IFREQ >= 1 .AND. IFREQ <= NFREQ) THEN
                FREQ(IFREQ, ROW_VAL, COL_VAL) = FREQ_VAL
                AMAS(IFREQ, ROW_VAL, COL_VAL) = AMAS_VAL
                BDMP(IFREQ, ROW_VAL, COL_VAL) = BDMP_VAL
             ENDIF
          ENDDO
      
          CLOSE(901)
          NFREQ_FOUND = IFREQ
      
          ! Write reordered data back to original file
          OPEN(UNIT=902, FILE=TRIM(filename), STATUS='REPLACE', ACTION='WRITE', IOSTAT=IOS)
          IF (IOS /= 0) THEN
             DEALLOCATE(FREQ, AMAS, BDMP, FREQ_LIST)
             RETURN
          ENDIF
      
          DO I = 1, NFREQ_FOUND
             DO IROW = 1, NDOF          ! Row as outer loop
                DO ICOL = 1, NDOF       ! Column as inner loop
                   WRITE(902, '(ES14.6, 2I6, 2ES14.6)') &
                      FREQ_LIST(I), IROW, ICOL, &
                      AMAS(I, IROW, ICOL), BDMP(I, IROW, ICOL)
                ENDDO
             ENDDO
          ENDDO
      
          CLOSE(902)
      
          WRITE(*,*) ' AmssDamp.1 has been rearranged to row-major format.'
      
          ! Deallocate
          DEALLOCATE(FREQ, AMAS, BDMP, FREQ_LIST)
      
      END SUBROUTINE ReorderAmssDamp

END MODULE PotentWavForceMulti
!*******************************************************************************
