module test_integration_config1
  
    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed, skip_test
    use testing_utilities

    implicit none
    public :: setup_config1, collect_tests_config1
    character(len=*), parameter :: test_data_dir = "test/data/config1/"
    logical :: simulation_execution_failed = .false. ! Flag to control further tests

    contains

    !> Collect all exported tests for config1
    subroutine setup_config1(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("run_simulation_config1", run_simulation_config1) &
        ]
    end subroutine setup_config1

    !> Collect all exported tests for config1
    subroutine collect_tests_config1(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("test_oamass_config1", test_oamass_config1) &
            !new_unittest("test_odamping_config1", test_odamping), &
            !new_unittest("test_oexfor_config1", test_oexfor) &
        ]
    end subroutine collect_tests_config1

    subroutine run_simulation_config1(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: inputdir, outputdir, command, errcode, simlog
        integer :: exitcode
        
        inputdir = get_repo_root_path()//test_data_dir//"input/"
        outputdir = get_running_exe_path()//"output_config1/"
        
        print*, "Running simulation config1 using '"//inputdir//"' as input directory."

        command = hamsexe_path//" "//inputdir//" "//outputdir//" > "//get_running_exe_path()//"config1.log 2>&1"
        
        ! Run the simulation
        print *, "Waiting for simulation config1 to complete", command
        call execute_command_line(command, exitstat=exitcode, wait=.true.)
        if (exitcode /= 0) then
            simulation_execution_failed = .true.
            write(errcode, '(I0)') exitcode
            call test_failed(error, "Error: could not execute HAMS for integration test config1. Failed with exit code: "//errcode)
        end if
    end subroutine run_simulation_config1

    subroutine test_oamass_config1(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: filename_expected, filename_actual
        !real(kind=8), allocatable :: expected(:,:), actual(:,:)
        character(len=1) :: istr
        integer :: i

        ! Skip test is simulation did not succeed
        if (simulation_execution_failed) then
            call skip_test(error, "Skipping OAMASS comparison for config1 because simulation failed")
            return
        end if

        do i = 1, 6
            ! Convert i to str
            write(istr, '(I1)') i
            ! Get full paths to expected and actual oamass datafiles
            filename_expected = get_repo_root_path()//test_data_dir//'expected_output/OAMASS'//istr//'.txt'
            filename_actual = get_running_exe_path()//'output_config1/Hams_format/OAMASS'//istr//'.txt'
            ! Check data in files is equal
            call test_oamass_files_are_equal(filename_expected, filename_actual, 1.0d-6, error)
            if (allocated(error)) return   ! stop further comparisons on first failure
        end do

    end subroutine test_oamass_config1





end module test_integration_config1