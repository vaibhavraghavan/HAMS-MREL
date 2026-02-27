!  ------------------------------------------------------------------------------------------------------
!                                                               
!    Program HAMS-MREL for the diffraction and radiation of waves 
!    by 3D structures.
! 
!  License:
! 
!    This routine is part of HAMS and HAMS-MREL.
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
!    Yingyi Liu
!
!  ------------------------------------------------------------------------------------------------------
MODULE ReadPanelMesh

   USE Body_mod
   USE PanelMesh_mod
   USE Inerfs_mod
   USE NormalProcess
   
   IMPLICIT NONE

      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: ReadBodyMesh
   PUBLIC :: ReadWTPLMesh
   
   PUBLIC :: CalNormals

CONTAINS
!=======================================================================
!SUBROUTINE ReadPanelMesh

!   !   ReadPanelMesh is used to read nodal and elemental data on the 
!   !   body surface of an arbitrary body

      SUBROUTINE ReadBodyMesh
      IMPLICIT   NONE

      INTEGER M,N,IND,IEL,J
      REAL*8::  DX,DY,DZ

        DO 10 IND=1,NTND                                                              ! The 10 is a label. Using this, statements to be continued if the loop does not progress can be mentioned with the CONTINUE command. Reference: http://www.personal.psu.edu/jhm/f90/lectures/15.html
          READ(2,*) M, XYZ(IND,1), XYZ(IND,2), XYZ(IND,3)
10      CONTINUE

         DO J=1,3
          READ(2,*)
         ENDDO
         
        DO 30  IEL=1, NELEM
	    READ(2, *) M, NCN(IEL), (NCON(IEL,J), J=1, NCN(IEL))
30      CONTINUE

      RETURN
      END SUBROUTINE ReadBodyMesh

!=======================================================================
!SUBROUTINE ReadWTPLMesh


!   !   ReadWTPLMesh is used to read nodal and elemental data on the                     
!   !   inner water plane

      SUBROUTINE ReadWTPLMesh
      IMPLICIT   NONE  

      INTEGER M,N,IND,IEL,J
      REAL*8::  DX,DY,DZ
!
! -------------------------------------------------------------------------
! 
      DO 10 IND=1,iNTND
        READ(5,*) M, iXYZ(IND,1), iXYZ(IND,2), iXYZ(IND,3)
        IF (ABS(iXYZ(IND,3)).GT.1.E-10) THEN
         Print *,' Error: Z Coordinate is not zero at Node No.',IND
         STOP
        ENDIF
10    CONTINUE

         DO J=1,3
          READ(5,*)
         ENDDO
         
      DO 30  IEL=1, iNELEM
	    READ(5, *) M, iNCN(IEL),(iNCON(IEL,J), J=1, iNCN(IEL))    
30    CONTINUE

       RETURN
       END SUBROUTINE ReadWTPLMesh

        
!=======================================================================
!SUBROUTINE CalNormals


!   !   CalNormals is used to calculate normals on the immersed body surface 
!   !   and on the inner water plane


      SUBROUTINE CalNormals(IFLAG,INPUTDIR)
      IMPLICIT   NONE

      INTEGER IEL,IND
      INTEGER,INTENT(IN):: IFLAG
      CHARACTER(LEN=*),INTENT(IN):: INPUTDIR

      ! Local variables for generalized modes
      INTEGER:: GEN_DOF, GEN_PID, GEN_IEL, GEN_IOST
      REAL*8:: GEN_NORMAL
      CHARACTER(LEN=512):: GEN_FILENAME
      CHARACTER(LEN=10):: DOF_STR
      LOGICAL:: GEN_FILE_EXISTS

! -------------------------------------------------------------------------
! 
!      Calculate some panel properties
!
      CALL CalPanelCentre( XYZ, NTND, NELEM, NCN, NCON, XYZ_P)

      CALL CalPanelArea( XYZ, NTND, NELEM, NCN, NCON, DS)

      CALL CalPanelChartLength( XYZ, NTND, NELEM, NCN, NCON, PNSZ)

      CALL CalTransNormals( XYZ, NTND, NELEM, NCN, NCON, DXYZ_P)

      CALL CalRotNormals(XR, XYZ_P, DXYZ_P, NELEM)

      IF (IFLAG.NE.0) THEN                                                            ! IFLAG is used as an indicator for irregular frequency removal. If 1, these will be removed

         CALL CalPanelCentre(iXYZ,INTND,INELEM,INCN,INCON,iXYZ_P)

         CALL CalPanelArea(iXYZ,INTND,INELEM,INCN,INCON,IDS)

         CALL CalPanelChartLength(iXYZ,INTND,INELEM,INCN,INCON,IPNSZ)

         CALL CalTransNormals(iXYZ,INTND,INELEM,INCN,INCON,IDXYZ_P)

         CALL CalRotNormals(XR,iXYZ_P,IDXYZ_P,INELEM)
      ENDIF
!
! -------------------------------------------------------------------------
!     Read generalized mode normals from gen_mod files if present
!     File naming convention: gen_mod_1_{dof}.txt  (body number is always 1)
!     File format per line:   panel_id    normal_value
! -------------------------------------------------------------------------

      DO GEN_DOF = 1, 6
        WRITE(DOF_STR,'(I0)') GEN_DOF
        GEN_FILENAME = TRIM(INPUTDIR)//'/gen_mod_1_'// &
                        TRIM(ADJUSTL(DOF_STR))//'.txt'

        INQUIRE(FILE=TRIM(GEN_FILENAME), EXIST=GEN_FILE_EXISTS)

        IF (GEN_FILE_EXISTS) THEN
          PRINT *, ' Found generalized mode file: ', TRIM(GEN_FILENAME)
          PRINT *, '   Overriding DOF ', GEN_DOF, ' normals for body 1'

          OPEN(99, FILE=TRIM(GEN_FILENAME), STATUS='OLD', ACTION='READ')

          DO GEN_IEL = 1, NELEM
            READ(99, *, IOSTAT=GEN_IOST) GEN_PID, GEN_NORMAL
            IF (GEN_IOST .NE. 0) THEN
              PRINT *, '  Warning: Could not read line ', GEN_IEL, &
                       ' from ', TRIM(GEN_FILENAME)
              EXIT
            ENDIF
            DXYZ_P(GEN_PID, GEN_DOF) = GEN_NORMAL
          ENDDO

          CLOSE(99)
          PRINT *, '   Loaded ', NELEM, ' generalized normals.'
        ENDIF
      ENDDO

      Print *,' Calculating panel normals is finished...'
      Print * 
!

      RETURN
      END SUBROUTINE CalNormals
      
!-------------------------------------------------------------------------------
    END MODULE ReadPanelMesh
!*******************************************************************************
