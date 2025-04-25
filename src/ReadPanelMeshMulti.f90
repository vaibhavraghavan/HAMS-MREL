!  ------------------------------------------------------------------------------------------------------
!                                                               
!    Program HAMS for the diffraction and radiation of waves 
!    by 3D structures.
! 
!  License:
! 
!    This routine is part of HAMS. This deals with reading the panel meshes for the multibodies.
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
MODULE ReadPanelMeshMulti

   USE Body_mod
   USE PanelMesh_mod
   USE Inerfs_mod
   USE NormalProcess
   
   IMPLICIT NONE

      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: ReadBodyMeshMulti
   PUBLIC :: ReadWTPLMeshMulti
   !
   !PUBLIC :: CalNormals                                           ! This subroutine is already available in ReadPanelMesh

CONTAINS
!=======================================================================
!SUBROUTINE ReadPanelMeshMulti

!   !   ReadPanelMeshMulti is used to read nodal and elemental data on the 
!   !   body surface of each arbitrary body

      SUBROUTINE ReadBodyMeshMulti(BODY_N,NBODY,LCS_MULTI)
      IMPLICIT   NONE

      INTEGER M,IND,IEL,J
      INTEGER,INTENT(IN):: BODY_N,NBODY
      REAL*8,INTENT(IN):: LCS_MULTI(NBODY,4)                                   !Array of location of the Local Coordinate System (LCS) along with the rotation of the LCS w.r.t the Global Coordinate System (GCS) 
      
        DO 10 IND=1,NTND_MULTI(BODY_N)                                                              ! The 10 is a label. Using this, statements to be continued if the loop does not progress can be mentioned with the CONTINUE command. Reference: http://www.personal.psu.edu/jhm/f90/lectures/15.html
          READ(2,*) M, XYZ_LOCAL_MULTI(BODY_N,IND,1), XYZ_LOCAL_MULTI(BODY_N,IND,2), XYZ_LOCAL_MULTI(BODY_N,IND,3)
10      CONTINUE
        
        ! Creating the array containing the Global coordinates
        DO IND=1,NTND_MULTI(BODY_N)                                                   
          XYZ_GLOBAL_MULTI(BODY_N,IND,1) = XYZ_GLOBAL_MULTI(BODY_N,IND,1)+LCS_MULTI(BODY_N,1)+COS(LCS_MULTI(BODY_N,4))*XYZ_LOCAL_MULTI(BODY_N,IND,1)-SIN(LCS_MULTI(BODY_N,4))*XYZ_LOCAL_MULTI(BODY_N,IND,2)
          XYZ_GLOBAL_MULTI(BODY_N,IND,2) = XYZ_GLOBAL_MULTI(BODY_N,IND,2)+LCS_MULTI(BODY_N,2)+SIN(LCS_MULTI(BODY_N,4))*XYZ_LOCAL_MULTI(BODY_N,IND,1)+COS(LCS_MULTI(BODY_N,4))*XYZ_LOCAL_MULTI(BODY_N,IND,2)
          XYZ_GLOBAL_MULTI(BODY_N,IND,3) = XYZ_LOCAL_MULTI(BODY_N,IND,3) ! This needs to be modified so any z-coordinate can be given
        ENDDO
        
         DO J=1,3
          READ(2,*)
         ENDDO
         
        DO 30  IEL=1, NELEM_MULTI(BODY_N)
         READ(2, *) M, NCN_MULTI(BODY_N,IEL), (NCON_MULTI(BODY_N,IEL,J), J=1, NCN_MULTI(BODY_N,IEL))
30      CONTINUE

      RETURN
      END SUBROUTINE ReadBodyMeshMulti

!=======================================================================
!SUBROUTINE ReadWTPLMeshMulti


!   !   ReadWTPLMesh is used to read nodal and elemental data on the                     
!   !   inner water plane for multi-bodies

      SUBROUTINE ReadWTPLMeshMulti(BODY_N,NBODY,LCS_MULTI)
      IMPLICIT   NONE  

      INTEGER M,N,IND,IEL,J
      INTEGER,INTENT(IN):: BODY_N,NBODY
      REAL*8,INTENT(IN):: LCS_MULTI(NBODY,4)                                   !Array of location of the Local Coordinate System (LCS) along with the rotation of the LCS w.r.t the Global Coordinate System (GCS) 
!
! -------------------------------------------------------------------------
! 
      DO 10 IND=1,iNTND_MULTI(BODY_N)
        READ(5,*) M, iXYZ_LOCAL_MULTI(BODY_N,IND,1), iXYZ_LOCAL_MULTI(BODY_N,IND,2), iXYZ_LOCAL_MULTI(BODY_N,IND,3)
        IF (ABS(iXYZ_LOCAL_MULTI(BODY_N,IND,3)).GT.1.E-10) THEN
         Print *,' Error: Z Coordinate is not zero at Node No.',IND
         STOP
        ENDIF
10    CONTINUE
      
      ! Creating the array containing the Global coordinates
        DO IND=1,iNTND_MULTI(BODY_N)                                                   
          iXYZ_GLOBAL_MULTI(BODY_N,IND,1) = iXYZ_GLOBAL_MULTI(BODY_N,IND,1)+LCS_MULTI(BODY_N,1)+COS(LCS_MULTI(BODY_N,4))*iXYZ_LOCAL_MULTI(BODY_N,IND,1)-SIN(LCS_MULTI(BODY_N,4))*iXYZ_LOCAL_MULTI(BODY_N,IND,2)
          iXYZ_GLOBAL_MULTI(BODY_N,IND,2) = iXYZ_GLOBAL_MULTI(BODY_N,IND,2)+LCS_MULTI(BODY_N,2)+SIN(LCS_MULTI(BODY_N,4))*iXYZ_LOCAL_MULTI(BODY_N,IND,1)+COS(LCS_MULTI(BODY_N,4))*iXYZ_LOCAL_MULTI(BODY_N,IND,2)
          iXYZ_GLOBAL_MULTI(BODY_N,IND,3) = iXYZ_LOCAL_MULTI(BODY_N,IND,3) ! This needs to be modified so any z-coordinate can be given
        ENDDO

         DO J=1,3
          READ(5,*)
         ENDDO
         
      DO 30  IEL=1, iNELEM_MULTI(BODY_N)
         READ(5, *) M, iNCN_MULTI(BODY_N,IEL), (iNCON_MULTI(BODY_N,IEL,J), J=1, iNCN_MULTI(BODY_N,IEL))
30    CONTINUE

       RETURN
       END SUBROUTINE ReadWTPLMeshMulti

        
!=======================================================================
!SUBROUTINE CalNormalsMulti


!   !   CalNormalsMulti is used to calculate normals on the immersed body surface 
!   !   and on the inner water plane for all bodies


      SUBROUTINE CalNormalsMulti(BODY_N,IFLAG)
      IMPLICIT   NONE

      INTEGER IEL,IND
      INTEGER,INTENT(IN):: BODY_N,IFLAG

! -------------------------------------------------------------------------
! 
!      Calculate some panel properties
      
      CALL CalPanelCentre(XYZ_GLOBAL_MULTI(BODY_N,1:NTND_MULTI(BODY_N),:), NTND_MULTI(BODY_N), NELEM_MULTI(BODY_N), NCN_MULTI(BODY_N,1:NELEM_MULTI(BODY_N)), NCON_MULTI(BODY_N,1:NELEM_MULTI(BODY_N),:), XYZ_MULTI_P(BODY_N,1:NELEM_MULTI(BODY_N),:))
      
      CALL CalPanelArea(XYZ_GLOBAL_MULTI(BODY_N,1:NTND_MULTI(BODY_N),:), NTND_MULTI(BODY_N), NELEM_MULTI(BODY_N), NCN_MULTI(BODY_N,1:NELEM_MULTI(BODY_N)), NCON_MULTI(BODY_N,1:NELEM_MULTI(BODY_N),:), DS_MULTI(BODY_N,1:NELEM_MULTI(BODY_N)))
      
      CALL CalPanelChartLength(XYZ_GLOBAL_MULTI(BODY_N,1:NTND_MULTI(BODY_N),:), NTND_MULTI(BODY_N), NELEM_MULTI(BODY_N), NCN_MULTI(BODY_N,1:NELEM_MULTI(BODY_N)), NCON_MULTI(BODY_N,1:NELEM_MULTI(BODY_N),:), PNSZ_MULTI(BODY_N,1:NELEM_MULTI(BODY_N)))
      
      CALL CalTransNormals(XYZ_GLOBAL_MULTI(BODY_N,1:NTND_MULTI(BODY_N),:), NTND_MULTI(BODY_N), NELEM_MULTI(BODY_N), NCN_MULTI(BODY_N,1:NELEM_MULTI(BODY_N)), NCON_MULTI(BODY_N,1:NELEM_MULTI(BODY_N),:), DXYZ_MULTI_P(BODY_N,1:NELEM_MULTI(BODY_N),:))
      
      CALL CalRotNormals(XR_MULTI(BODY_N,:), XYZ_MULTI_P(BODY_N,1:NELEM_MULTI(BODY_N),:), DXYZ_MULTI_P(BODY_N,1:NELEM_MULTI(BODY_N),:), NELEM_MULTI(BODY_N))                    ! The center of gravity should be used here, since the rotational normals are w.r.t COG

      IF (IFLAG.NE.0) THEN                                                            ! IFLAG is used as an indicator for irregular frequency removal. If 1, these will be removed
       IF (iNELEM_MULTI(BODY_N).GT.0) THEN
      
        CALL CalPanelCentre(iXYZ_GLOBAL_MULTI(BODY_N,1:iNTND_MULTI(BODY_N),:), iNTND_MULTI(BODY_N), iNELEM_MULTI(BODY_N), iNCN_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N)), iNCON_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N),:), iXYZ_MULTI_P(BODY_N,1:iNELEM_MULTI(BODY_N),:))
      
        CALL CalPanelArea(iXYZ_GLOBAL_MULTI(BODY_N,1:iNTND_MULTI(BODY_N),:), iNTND_MULTI(BODY_N), iNELEM_MULTI(BODY_N), iNCN_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N)), iNCON_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N),:), iDS_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N)))
      
        CALL CalPanelChartLength(iXYZ_GLOBAL_MULTI(BODY_N,1:iNTND_MULTI(BODY_N),:), iNTND_MULTI(BODY_N), iNELEM_MULTI(BODY_N), iNCN_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N)), iNCON_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N),:), iPNSZ_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N)))
      
        CALL CalTransNormals(iXYZ_GLOBAL_MULTI(BODY_N,1:iNTND_MULTI(BODY_N),:), iNTND_MULTI(BODY_N), iNELEM_MULTI(BODY_N), iNCN_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N)), iNCON_MULTI(BODY_N,1:iNELEM_MULTI(BODY_N),:), iDXYZ_MULTI_P(BODY_N,1:iNELEM_MULTI(BODY_N),:))
      
        CALL CalRotNormals(XR_MULTI(BODY_N,:), iXYZ_MULTI_P(BODY_N,1:iNELEM_MULTI(BODY_N),:), iDXYZ_MULTI_P(BODY_N,1:iNELEM_MULTI(BODY_N),:), iNELEM_MULTI(BODY_N))                    ! The center of gravity should be used here, since the rotational normals are w.r.t COG
       
       ENDIF
      ENDIF 
!

      RETURN
      END SUBROUTINE CalNormalsMulti
      
!-------------------------------------------------------------------------------
    END MODULE ReadPanelMeshMulti
!*******************************************************************************
