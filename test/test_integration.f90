module test_integration_config1
  
    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
    use testing_utilities

    implicit none
    public :: collect_tests_integration_config1
    character(len=*), parameter :: test_data_dir = "test/data/config1/"

    contains

    !> Collect all exported tests for config1
    subroutine collect_tests_integration_config1(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("run_simulation_config1", run_simulation_config1), &
            new_unittest("test_oamass_config1", test_oamass_config1) &
            !new_unittest("test_odamping_config1", test_odamping), &
            !new_unittest("test_oexfor_config1", test_oexfor) &
        ]
    end subroutine collect_tests_integration_config1

    subroutine run_simulation_config1(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: command, errcode
        integer :: exitcode
        
        print*, "Running simulation config1 using test/data/config1/input"
        
        command = hamsexe_path//" "//get_repo_root_path()//test_data_dir//"input "//get_running_exe_path()//"config1_output"
        print *, "COMMAND: ", command
        ! Run the simulation
        call execute_command_line(command, exitstat=exitcode)
        if (exitcode /= 0) then
            write(errcode, '(I0)') exitcode
            call test_failed(error, "Error: could not execute HAMS for integration test config1. Failed with exit code: "//errcode)
        end if

    end subroutine run_simulation_config1

    subroutine test_oamass_config1(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: filename_expected, filename_actual
        real(kind=8), allocatable :: expected(:,:), actual(:,:)
        character(len=1) :: istr
        integer :: i

        do i = 1, 6
            ! Convert i to str
            write(istr, '(I1)') i
            ! Get full paths to expected and actual oamass datafiles
            filename_expected = get_repo_root_path()//test_data_dir//'expected_output/OAMASS'//istr//'.txt'
            filename_actual = get_running_exe_path()//'output/Hams_format/OAMASS'//istr//'.txt'
            ! Load data into matrices
            call read_matrix_with_format(filename_expected, '(F7.3,1X,F7.3,1X,6E14.5)', 8, expected)
            call read_matrix_with_format(filename_actual, '(F7.3,1X,F7.3,1X,6E14.5)', 8, actual)
            print *, "Actual Size: ", size(actual)
            print *, "Excpected Size: ", size(expected)

            
            !print *, trim(filename_expected)

            !get_test_data_path('OAMASS'//istr//'.txt')
            !print*, trim(filename_expected)
            !call read_matrix_with_format(trim(filename_expected), '(F7.3,1X,F7.3,1X,6E14.5)', 8, expected)
            !print*, expected(:,1)
        end do


    end subroutine test_oamass_config1





end module test_integration_config1