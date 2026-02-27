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
!  Code Original Authors:
!
!    Vaibhav Raghavan and Yingyi Liu
!
!  ------------------------------------------------------------------------------------------------------

!  -------------------------------------------------------------------------------------------
!      Open files, and read the main input file
!---------------------------------------------------------------------------------------------
! 

module IO

    implicit none

    public :: ReadControlFile
    public :: VerifyInputFilesExist
    public :: CreateOutputFiles
    public :: CreateErrorCheckFile

    ! Specify if WK1 and DWK were intialized (for ErrorCheck.txt)
    logical, private :: allocated_WK1_and_DWK

    ! Write one PressureElevation.6p file per body 
    logical, public :: separate_wamit_diffraction_radiation_files

    contains

    subroutine ReadControlFile(dir, success)
        
        use HAMS_mod
        use Body_mod
        use WaveDyn_mod
        use FieldOutput_mod
      
        implicit none
        character(len=*), intent(in) :: dir
        logical, intent(out) :: success
        logical :: exists

        integer I,J,err,IFS,NPET,NB
        character(len=100) FILE_RADIATION,FILE_NUMBER,BODY_CHECK,line

        success = .false.

        ! Open input file
        inquire(file=dir//"/ControlFile.in", exist=exists)
        if (exists) then
            open(1, file=dir//"/ControlFile.in", status="OLD", action="READ")
        else
            print *, "Error opening ControlFile.in for reading"
            return
        end if

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
        
        if (SYBO == 0) then
            IFS = 0
        else if (SYBO == 1) then
            IFS = 2
        else
            print*, 'Warning: SYBO must be 0 or 1.', new_line('a')
            IFS = 0
        end if
        
        if (NPET > 0) then
            allocated_WK1_and_DWK = .false.
            NPER = IFS+NPET
            allocate(WVNB(NPER))
            read(1,*) (WVNB(I), I = IFS+1, NPER)
        else if (NPET < 0) then
            allocated_WK1_and_DWK = .true.
            NPET = ABS(NPET)
            NPER = IFS+NPET
            allocate(WVNB(NPER))
            read(1,'(27x,f30.15)') WK1
            read(1,'(19x,f30.15)') DWK
            do I= IFS+1, NPER
                WVNB(I) = WK1 + (I-IFS-1)*DWK
            end do
        end if
        
        read(1,*)
        read(1,*)
        read(1,'(A100)') BODY_CHECK
        
        ! Reading of the number of multi-bodies
        if (INDEX(BODY_CHECK,"multi") > 0) then
            READ(1,'(24x,i16)') NBODY
        else 
            NBODY = 1
        end if
        
        if (NBODY < 0) then
            print*, "ERROR: The number of bodies must be greater than or equal to 1.", new_line('a')
            print*, "Terminating application."
            stop
        end if
        
        ! Obtaining the origin coordinates of the LCS per mesh as well as the rotation w.r.t GCS
        if (NBODY > 1) then
            allocate(LCS_MULTI(NBODY,4))
            do NB = 1, NBODY
                read(1,'(26x,4f12.3)')      (LCS_MULTI(NB,I), I=1,4)
            end do
        end if
        if (NBODY > 1) then        
            read(1,*) 
            read(1,*)
            read(1,*) 
        end if
        
        read(1,'(23x,i16)') NBETA
        
        if (NBETA > 0) then
            allocate(WVHD(NBETA))
            read(1,*) (WVHD(I),I=1,NBETA)
        else if (NBETA < 0) then
            NBETA = ABS(NBETA)
            allocate(WVHD(NBETA))
            read(1,'(20x,f30.15)') BETA1
            read(1,'(17x,f30.15)') DBETA
            do I = 1, NBETA
                WVHD(I)=BETA1+(I-1)*DBETA
            end do
        end if
        
        read(1,*) 
        read(1,*)

        ! Read rotation centers
        if (NBODY == 1) then
            read(1,'(28x,3f12.3)') (XR(I), I=1,3)                          
        else if (NBODY > 1) then
            allocate(XR_MULTI(NBODY,3))
            do NB = 1, NBODY
                read(1,'(28x,3f12.3)') (XR_MULTI(NB,I), I=1,3)
            end do
        end if

        read(1,'(26x,f30.15)')     REFL
        read(1,'(26x,i16)')        ISOL
        read(1,'(23x,i16)')        IRSP
        read(1,'(23x,i16)')        NTHREAD

        ! Field Points
        read(1,*) 
        read(1,*)
        ! Read number of field points
        read(1,'(27x,i16)')        NFP
        ! Read field points
        allocate(XFP(NFP,3))
        do I = 1,NFP
            ! READ(1,'(26x,3(1x,f10.4))')     (XFP(I,J), J=1,3)
            read(1,*) (XFP(I,J), J=1,3)
        end do

        ! Wamit output config
        read(1,*)
        read(1,*)
        ! Check if user defined "separate_wamit_diffraction_radiation_files'
        read(1,'(A)') line
        if (index(line, "separate") > 0 .or. index(line, "Separate") > 0) then
            ! Parse optional line
            read(line, '(50x,i16)') separate_wamit_diffraction_radiation_files
        else
            ! User did not define option
            separate_wamit_diffraction_radiation_files = .false.
        end if
        if (NBODY == 1) then
            ! Always false for single-body simulationas.
            separate_wamit_diffraction_radiation_files = .false.
        end if

        success = .true.

    end subroutine ReadControlFile

    subroutine VerifyInputFilesExist(dir, numbodies, irrfreq, success)
        implicit none
        character(len=*), intent(in) :: dir
        integer, intent(in) :: numbodies, irrfreq
        logical, intent(out) :: success
        logical :: hullexists, hydroexists, wpmexists
        integer :: i
        character(len=5) :: istr ! len=5 assumes <= 99999 bodies

        success = .true.
        wpmexists = .true.

        if (numbodies == 1) then
            istr = '1'
            inquire(file=dir//"/HullMesh.pnl", exist=hullexists)
            inquire(file=dir//"/Hydrostatic.in", exist=hydroexists)
            if (irrfreq .ne. 0) then
                inquire(file=dir//"/WaterPlaneMesh.pnl", exist=wpmexists)
            end if
        else if (numbodies > 1) then
            do i = 1, numbodies
                write(istr, '(I5)') i   ! Convert i to str
                inquire(file=dir//"/HullMesh_"//trim(adjustl(istr))//".pnl", exist=hullexists)
                inquire(file=dir//"/Hydrostatic_"//trim(adjustl(istr))//".in", exist=hydroexists)
                if (irrfreq .ne. 0) then
                    inquire(file=dir//"/WaterPlaneMesh_"//trim(adjustl(istr))//".pnl", exist=wpmexists)
                end if
                if ((.not. hullexists) .or. (.not. hydroexists) .or. (.not. wpmexists)) then
                    exit
                end if
            end do
        end if

        if (.not. hullexists) then 
            print*, "HullMesh.pnl input file missing for body ", istr, new_line('a')
            success = .false.
        end if
        if (.not. hydroexists) then
            print*, "Hydrostatic.in input file missing for body ", istr, new_line('a')
            success = .false.
        end if
        if (.not. wpmexists) then
            print *, "WaterPlaneMesh.pnl input file missing for body ", istr, new_line('a')
            success = .false.
        end if

    end subroutine VerifyInputFilesExist

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

        ! IDs used in open statements are also used in write statements
        ! Don't modify their values or write statements will not find the files
        do i = 1, 6
            write(istr, '(I1)') i   ! Convert i to str
            open(20+i, file=dir//'/OEXFOR'//istr//'.txt', status='UNKNOWN')
            open(30+i, file=dir//'/OAMASS'//istr//'.txt', status='UNKNOWN')
            open(40+i, file=dir//'/ODAMPING'//istr//'.txt', status='UNKNOWN')
        end do
    end subroutine CreateHamsFiles

    subroutine CreateWamitFiles(dir, numbodies)
        implicit none
        character(len=*), intent(in) :: dir
        integer, intent(in) :: numbodies
        integer :: i
        character(len=5) :: istr ! len=5 assumes <= 99999 bodies

        ! IDs used in open statements are also used in write statements
        ! Don't modify their values or write statements will not find the files
        open(61, file=dir//'/AmssDamp.1', status='UNKNOWN')
        open(62, file=dir//'/ExcForce.3', status='UNKNOWN')
        open(63, file=dir//'/Motion.4', status='UNKNOWN')
        open(65, file=dir//'/Hydrostat.hst', status='UNKNOWN')
        if (numbodies > 1 .AND. separate_wamit_diffraction_radiation_files) then
            open(639, file=dir//'/PressureElevationDiffraction.6p', status='UNKNOWN')
            do i = 1, numbodies
                write(istr, '(I5)') i   ! Convert i to str
                open(640+i, file=dir//'/PressureElevationRadiation_'//trim(adjustl(istr))//'.6p', status='UNKNOWN')
            end do
        else
            open(640, file=dir//'/PressureElevation.6p', status='UNKNOWN')
        end if
        if (numbodies > 1) then
            open(66, file=dir//'/PressureElevationIncidence.6p', status='UNKNOWN')
        end if
    end subroutine CreateWamitFiles

    subroutine CreateHydrostarFiles(dir)
        implicit none
        character(len=*), intent(in) :: dir
        character(len=1) :: istr, jstr
        integer :: i, j

        ! IDs used in open statements are also used in write statements
        ! Don't modify their values or write statements will not find the files
        do i = 1, 6
            write(istr, '(I1)') i   ! Convert i to str
            do j = 1, 6
                write(jstr, '(I1)') j   ! Convert j to str
                open((60+(i*10))+j, file=dir//'/AddedMass_'//istr//jstr//'.rao', status='UNKNOWN')
                open((120+(i*10))+j, file=dir//'/WaveDamping_'//istr//jstr//'.rao', status='UNKNOWN')
            end do
            open(190+i, file=dir//'/Excitation_'//istr//'.rao', status='UNKNOWN')
            open(200+i, file=dir//'/Motion_'//istr//'.rao', status='UNKNOWN')
            end do
    end subroutine CreateHydrostarFiles

    subroutine CreateErrorCheckFile(dir, numbodies)
        use Body_mod, only : XR, XR_MULTI, XG, XG_MULTI     ! WavDynMods.f90
        use HAMS_mod, only : WK1, DWK                   ! WavDynMods.f90
        use WaveDyn_mod, only : H                       ! WavDynMods.f90

        implicit none
        character(len=*), intent(in) :: dir
        integer, intent(in) :: numbodies
        integer :: i, j

        open(9, file=dir//'/ErrorCheck.txt', status="UNKNOWN", action="WRITE")

        write(9,*) "Water depth: ", H, new_line('a')

        if (allocated_WK1_and_DWK .eqv. .true.) then
            write(9,*) "Starting Frequency: ", WK1, new_line('a')
            write(9,*) "Spacing Between Frequencies: ", DWK, new_line('a')
        end if

        if(numbodies == 1) then
            write(9,*) "The rotation center is input as (please confirm if it is correct):"
            write(9,'(3f12.3)') (XR(i), i=1,3)
            write(9,*)
            write(9,*) "The center of gravity is input as (please confirm if it is correct):"
            write(9,'(3f12.3)') (XG(i), i=1,3)
        else
            write(9,*) "The rotation center for the bodies are input as (please confirm if it is correct):"
            do i = 1, numbodies
                write(9,'(3f12.3)') (XR_MULTI(i,j), j=1,3)
            end do
            write(9,*) "The center of gravity per body is input as (please confirm if it is correct):"
            do i = 1, numbodies
                write(9,'(3f12.3)') (XG_MULTI(i,j), j=1,3)
            end do
        end if

        close(9)

    end subroutine CreateErrorCheckFile

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
            if (numbodies == 1) then
                call CreateDirectory(outputdir // "/Hydrostar_format", success)
            end if
        end if

        ! Exit if directories were not created
        if (.not. success) then
            success = .false.
            return
        end if
        
        call CreateHamsFiles(outputdir//"/Hams_format")
        call CreateWamitFiles(outputdir//"/Wamit_format", numbodies)
        if (numbodies == 1) then
            call CreateHydrostarFiles(outputdir // "/Hydrostar_format")
        end if
        
    end subroutine CreateOutputFiles

end module IO