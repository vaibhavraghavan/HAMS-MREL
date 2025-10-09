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
!---------------------------------------------------------------------------------------------
!  BLOCK-DIAGONAL ASSEMBLY: (NBODY,6,6)  --->  (NBODY*6, NBODY*6)
!---------------------------------------------------------------------------------------------
SUBROUTINE LOCALTOGLOBAL(LOCALMATRIX, NBODY, GLOBALMATRIX)
  IMPLICIT NONE
  INTEGER, INTENT(IN)            :: NBODY
  REAL*8,  INTENT(IN)            :: LOCALMATRIX(NBODY,6,6)
  REAL*8,  INTENT(OUT)           :: GLOBALMATRIX(6*NBODY,6*NBODY)
  INTEGER                        :: B, R0

  ! INITIALIZE THE GLOBAL MATRIX TO ZERO
  GLOBALMATRIX = 0.0D0

  ! ASSEMBLE BLOCK-DIAGONAL GLOBAL MATRIX
  DO B = 1, NBODY
     R0 = (B-1)*6
     GLOBALMATRIX(R0+1:R0+6, R0+1:R0+6) = LOCALMATRIX(B,:,:)
  END DO

END SUBROUTINE LOCALTOGLOBAL


!---------------------------------------------------------------------------------------------
!       Calculate the motion response in frequency domain by panel model.
!---------------------------------------------------------------------------------------------
SUBROUTINE SolveMotionMulti(W1,TP,OUFR,BETA,AMP,AMAS,BDMP, &
     BLNR,BQDR,EXFC,DSPL,NBODY)
  USE HAMS_mod
  USE Const_mod
  USE Body_mod
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: NBODY
  REAL*8,INTENT(IN)  :: W1,TP,OUFR,BETA,AMP
  REAL*8,INTENT(IN)  :: AMAS(NBODY*6,NBODY*6), BDMP(NBODY*6,NBODY*6), BLNR(NBODY,6,6), BQDR(NBODY,6,6)
  COMPLEX*16,INTENT(IN)  :: EXFC(NBODY*6)
  COMPLEX*16,INTENT(OUT) :: DSPL(NBODY*6)

  REAL*8 DLANGE

  INTEGER :: INFO, MD, MEXP, I, J, NB6
  INTEGER :: IPV(NBODY*6)
  INTEGER :: LD, BODYID
  REAL*8  :: BLNR_FINAL(NBODY*6,NBODY*6), BQDR_FINAL(NBODY*6,NBODY*6)
  REAL*8  :: NORM, WORK(NBODY*6), RERR
  REAL*8  :: MAGN, PHS, NREL, NIMG, NFAC
  REAL*8  :: VDMP(NBODY*6,NBODY*6)                 ! real equivalent damping from quadratic term
  COMPLEX*16 :: LEFT(NBODY*6,NBODY*6), RIGHT(NBODY*6), DSPL1(NBODY*6), DX(NBODY*6)

  NB6 = NBODY*6

  ! ========================================================
  CALL LOCALTOGLOBAL(BLNR, NBODY, BLNR_FINAL)
  CALL LOCALTOGLOBAL(BQDR, NBODY, BQDR_FINAL)

  ! Norm of quadratic damping (decide whether to iterate)
  NORM = DLANGE('M', NB6, NB6, BQDR_FINAL, NB6, WORK)

  ! ---------- Linear solve (no quadratic) ----------
  LEFT  = -W1**2*(MATX + AMAS) + CI*W1*(BDMP + BLNR_FINAL) + CRS + KSTF
  RIGHT = EXFC
  CALL ZGESV(NB6, 1, LEFT, NB6, IPV, RIGHT, NB6, INFO)
  DSPL  = RIGHT

  ! ---------- Iterate if quadratic term is non-negligible ----------
  IF (NORM .GE. 1.0D-6) THEN
     RERR = 100.0D0
     DO WHILE (RERR .GT. 1.0D-6)
        VDMP = 0.0D0
        DO I = 1, NB6
           DO J = 1, NB6
              VDMP(I,J) = BQDR_FINAL(I,J) * W1 * ABS(DSPL(J))
           END DO
        END DO

        LEFT  = -W1**2*(MATX + AMAS) + CI*W1*(BDMP + BLNR_FINAL + VDMP) + CRS + KSTF
        RIGHT = EXFC
        CALL ZGESV(NB6, 1, LEFT, NB6, IPV, RIGHT, NB6, INFO)
        DSPL1 = RIGHT

        DX = DSPL1 - DSPL
        RERR = 0.0D0
        DO J = 1, NB6
           IF (ABS(DSPL1(J)) .GT. 0.0D0) RERR = RERR + ABS(DX(J))/ABS(DSPL1(J))
        END DO
        DSPL = DSPL1
     END DO
  END IF

  !   =================================================== 
  !    Write WAMIT-style output files (ALL bodies & DOFs)
  !
  DO MD = 1, NB6
     ! Local DOF within the body (1..6) and body index (1..NBODY)
     ! LD: 1=surge, 2=sway, 3=heave, 4=roll, 5=pitch, 6=yaw
     LD     = MOD(MD-1, 6) + 1
     BODYID = (MD-1)/6 + 1

     IF (LD .LE. 3) THEN
        MEXP = 2
     ELSE
        MEXP = 3
     END IF

     NFAC = (RHO*G*AMP) * REFL**MEXP

     NREL = REAL(DSPL(MD)) / NFAC
     NIMG = AIMAG(DSPL(MD)) / NFAC
     MAGN  = SQRT(NREL*NREL + NIMG*NIMG)
     PHS  = ATAN2D(NIMG, NREL)

     IF (ABS(TP+1.0D0) .GT. 1.0D-6 .AND. ABS(TP) .GT. 1.0D-6) THEN
        ! Writes global DOF index MD (1..NB6). If you prefer BODYID/LD, adjust here.
        WRITE(63,1030) OUFR, BETA*180.0D0/PI, MD, MAGN, PHS, NREL, NIMG
     END IF
  END DO
  !   =================================================== 
1030 FORMAT(2ES14.6,I6,4ES14.6)

  RETURN
END SUBROUTINE SolveMotionMulti

