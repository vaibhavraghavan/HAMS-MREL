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

!  -------------------------------------------------------------------------------------------
!      Open files, and read the main input file
!---------------------------------------------------------------------------------------------
! 

module IO

    implicit none

    public :: ReadInputFile
    public :: CreateOutputFiles

    contains

    subroutine ReadInputFile(dir, success)
        
        use HAMS_mod
        use Body_mod
        use WaveDyn_mod
        use FieldOutput_mod
      
        implicit none
        character(len=*), intent(in) :: dir
        logical, intent(out) :: success
        logical :: exists

        integer I,J,err,IFS,NPET,NB
        character(len=100) FILE_RADIATION,FILE_NUMBER,BODY_CHECK

        success = .false.

        ! Open input file
        inquire(file=dir//"/ControlFile.in", exist=exists)
        if (exists) then
            open(1, file=dir//"/ControlFile.in", status="OLD", action="READ")
        else
            print *, "Error opening ControlFile.in for reading"
            return
        end if

        ! Open file to write rotation center
        open(9, file=dir//'/ErrorCheck.txt', status='UNKNOWN', action="WRITE")


        !  H :   Water depth; H<0: For infinite water depth; H>0: For finite water depth; 
        !  AMP:  Wave amplitude
        !  BETA: Wave incident angle with repect to the x-direction
        !  INFT=0,  input wave number; INFT=1, input wave frequency
        !  XC,YC,ZC: the coordinates of body rotation center or bouyancy center

        read(1,*) 
        read(1,*) 
        read(1,'(14x,f30.15)')     H
        read(1,*) 
        read(1,*) 
        read(1,'(27x,i16)')        SYBO
        read(1,'(25x,i16)')        INFT
        read(1,'(25x,i17)')        OUFT
        read(1,'(26x,i16)')        NPET
        
        if (SYBO.EQ.0) then
            IFS = 0
        else if (SYBO.EQ.1) then
            IFS = 2
        else
            print*, 'Warning: SYBO must be 0 or 1.'
            print*
            IFS = 0
        end if
        
        IF (NPET.GE.0) THEN
         NPER=IFS+NPET
         ALLOCATE(WVNB(NPER))
         READ(1,*) (WVNB(I),I=IFS+1,NPER)
        ELSEIF (NPET.LT.0) THEN
         NPET=ABS(NPET)
         NPER=IFS+NPET
         ALLOCATE(WVNB(NPER))
         READ(1,'(27x,f30.15)')     WK1
         READ(1,'(19x,f30.15)')     DWK
         DO I=IFS+1,NPER
          WVNB(I)=WK1+(I-IFS-1)*DWK
         ENDDO
        ENDIF
        
        READ(1,*)
        READ(1,*)
        READ(1,'(A100)') BODY_CHECK
        
        ! Reading of the number of multi-bodies
        IF (INDEX(BODY_CHECK,"multi")>0) THEN
            READ(1,'(24x,i16)')        NBODY
        ELSE 
            NBODY = 1
        ENDIF
        
        IF (NBODY.LT.0) THEN
         PRINT*, 'ERROR:The number of bodies must be greater than or equal to 1'
         PRINT*
         pause
         stop
        ENDIF
        ! Obtaining the origin coordinates of the LCS per mesh as well as the rotation w.r.t GCS
        IF (NBODY.GT.1) THEN
         ALLOCATE(LCS_MULTI(NBODY,4))
         DO NB=1,NBODY
          READ(1,'(26x,4f12.3)')      (LCS_MULTI(NB,I), I=1,4)
         ENDDO
        ENDIF
        IF (NBODY.GT.1) THEN        
        READ(1,*) 
        READ(1,*)
        READ(1,*) 
        ENDIF
        READ(1,'(23x,i16)')        NBETA
        IF (NBETA.GT.0) THEN
         ALLOCATE(WVHD(NBETA))
         READ(1,*) (WVHD(I),I=1,NBETA)
        ELSEIF (NBETA.LT.0) THEN
         NBETA=ABS(NBETA)
         ALLOCATE(WVHD(NBETA))
         READ(1,'(20x,f30.15)')     BETA1
         READ(1,'(17x,f30.15)')     DBETA
         DO I=1,NBETA
          WVHD(I)=BETA1+(I-1)*DBETA
         ENDDO
        ENDIF
        
        READ(1,*) 
        READ(1,*)
        IF (NBODY.EQ.1) THEN
         READ(1,'(30x,3f12.3)')      (XR(I), I=1,3)                          
         WRITE(9,*) "The rotation center is input as (please confirm if it is correct):" ! This is written into the Errorcheck.txt
         WRITE(9,'(3f12.3)') (XR(I), I=1,3)
        ELSEIF (NBODY.GT.1) THEN
         ALLOCATE(XR_MULTI(NBODY,3))
         WRITE(9,*) "The rotation center for the bodies are input as (please confirm if it is correct):" ! This is written into the Errorcheck.txt
         DO NB=1,NBODY
          READ(1,'(28x,3f12.3)')      (XR_MULTI(NB,I), I=1,3)
          WRITE(9,'(3f12.3)') (XR_MULTI(NB,I), I=1,3)
         ENDDO
        ENDIF
        READ(1,'(26x,f30.15)')     REFL
        READ(1,'(26x,i16)')        ISOL
        READ(1,'(23x,i16)')        IRSP
        READ(1,'(23x,i16)')        NTHREAD
        
        if (IRSP .NE. 0) then ! Checks if the waterplane mesh is there for a single body. For multi-bodies, this is added to the HAMS_Prog
            if (NBODY .EQ. 1) then
                open(5, file=dir//'/WaterPlaneMesh.pnl', status='OLD', iostat=err)
                if (err /= 0) then
                    print*, 'Error: The waterplane mesh file does not exist.'
                    print*
                    stop
                endif
            end if
        end if

        READ(1,*) 
        READ(1,*) 
        READ(1,'(27x,i16)')        NFP
        ALLOCATE(XFP(NFP,3))
        DO I=1,NFP
            ! READ(1,'(26x,3(1x,f10.4))')     (XFP(I,J), J=1,3)
           READ(1,*)     (XFP(I,J), J=1,3)
        ENDDO

        success = .true.

    end subroutine ReadInputFile

    subroutine CreateFile(filepath, success)
        implicit none
        character(len=*), intent(in) :: filepath
        logical, intent(out) :: success

        OPEN(856, FILE=filepath, STATUS='UNKNOWN')

        success = .true.

    end subroutine CreateFile

    subroutine CreateDirectory(dir, success)
        implicit none
        character(len=*), intent(in) :: dir
        logical, intent(out) :: success
        logical :: exists
        integer :: cstat
        character(100) :: cmsg

        inquire(directory=dir, exist=exists)
        if (.not. exists) then
            call execute_command_line("mkdir -p " // dir, CMDSTAT=cstat, CMDMSG=cmsg)
            if (cstat /= 0) then
                print *, "Failed to create directory ", dir, ". Error: ", trim(cmsg)
                success = .false.
                return
            end if
        end if

        success = .true.

    end subroutine CreateDirectory

    subroutine CreateHamsFiles(dir)
        implicit none
        character(len=*), intent(in) :: dir
        character(len=1) :: istr
        integer :: i

        do i = 1, 6
            write(istr, '(I1)') i   ! Convert i to str
            open(757, file=dir//'/OEXFOR'//istr//'.txt', status='UNKNOWN')
            open(758, file=dir//'/OAMASS'//istr//'.txt', status='UNKNOWN')
            open(759, file=dir//'/ODAMPING'//istr//'.txt', status='UNKNOWN')
        end do
    end subroutine CreateHamsFiles

    subroutine CreateWamitFiles(dir, numbodies)
        implicit none
        character(len=*), intent(in) :: dir
        integer, intent(in) :: numbodies
        character(len=5) :: istr ! len=5 assumes <= 99999 bodies
        integer :: i

        if (numbodies == 1) then
            open(61, file=dir//'/AmssDamp.1', status='UNKNOWN')
            open(62, file=dir//'/ExcForce.3', status='UNKNOWN')
            open(63, file=dir//'/Motion.4', status='UNKNOWN')
            open(64, file=dir//'/PressureElevation.6p', status='UNKNOWN')
            open(65, file=dir//'/Hydrostat.hst', status='UNKNOWN')
        else if (numbodies > 1) then
            open(61, file=dir//'/Buoy.1', status='UNKNOWN')
            open(62, file=dir//'/Buoy.3', status='UNKNOWN')
            open(63, file=dir//'/Buoy.4', status='UNKNOWN')
            open(64, file=dir//'/Buoy_Diffraction.6p', status='UNKNOWN')
            open(65, file=dir//'/Buoy.hst', status='UNKNOWN')
            open(66, file=dir//'/Buoy_Incidence.6p', status='UNKNOWN')
            do i = 1, numbodies
                write(istr, '(I5)') i
                open(200+i, file=dir//'/Buoy_Radiation_'//trim(adjustl(istr))//'.6p', status='UNKNOWN') 
            end do
        end if
    end subroutine CreateWamitFiles

    subroutine CreateHydrostarFiles(dir)
        implicit none
        character(len=*), intent(in) :: dir
        character(len=1) :: istr, jstr
        integer :: i, j

        do i = 1, 6
            write(istr, '(I1)') i   ! Convert i to str
            do j = 1, 6
                write(jstr, '(I1)') j   ! Convert j to str
                open(857, file=dir//'/AddedMass_'//istr//jstr//'.rao', status='UNKNOWN')
                open(858, file=dir//'/WaveDamping_'//istr//jstr//'.rao', status='UNKNOWN')
            end do
            open(859, file=dir//'/Excitation_'//istr//'.rao', status='UNKNOWN')
            open(859, file=dir//'/Motion_'//istr//'.rao', status='UNKNOWN')
            end do
    end subroutine CreateHydrostarFiles

    ! Create output files and directories
    subroutine CreateOutputFiles(outputdir, numbodies, success)
        implicit none
        character(len=*), intent(in) :: outputdir
        integer, intent(in) :: numbodies
        logical, intent(out) :: success

        ! Create output directories
        call CreateDirectory(outputdir, success)
        if (success) then
            call CreateDirectory(outputdir // "/Hams_format", success)
            call CreateDirectory(outputdir // "/Wamit_format", success)
            call CreateDirectory(outputdir // "/Hydrostar_format", success)
        end if

        ! Exit if directories were not created
        if (.not. success) then
            print *, "Could not create output directories."
            return
        end if

        call CreateHamsFiles(outputdir // "/Hams_format")
        call CreateWamitFiles(outputdir // "/Wamit_format", numbodies)
        if (numbodies == 1) then
            call CreateHydrostarFiles(outputdir // "/Hydrostar_format")
        end if
        
    end subroutine CreateOutputFiles

end module IO