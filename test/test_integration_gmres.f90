module test_integration_gmres

    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed, skip_test
    use testing_utilities

    implicit none
    public :: setup_gmres, collect_tests_gmres
    character(len=*), parameter :: test_data_dir = "test/data/config1_gmres/"
    character(len=*), parameter :: baseline_dir  = "test/data/config1/"
    logical :: simulation_execution_failed = .false.

    contains

    subroutine setup_gmres(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("run_simulation_gmres", run_simulation_gmres) &
        ]
    end subroutine setup_gmres

    subroutine collect_tests_gmres(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("test_oamass_gmres_vs_lu", test_oamass_gmres), &
            new_unittest("test_odamping_gmres_vs_lu", test_odamping_gmres) &
        ]
    end subroutine collect_tests_gmres

    subroutine run_simulation_gmres(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: inputdir, outputdir, command, errcode
        integer :: exitcode

        inputdir  = get_repo_root_path()//test_data_dir//"input/"
        outputdir = get_running_exe_path()//"output_gmres/"
        command   = hamsexe_path//" "//inputdir//" "//outputdir//" > "//get_running_exe_path()//"gmres.log 2>&1"

        print *, "Running simulation config1_gmres (ISOLV=2, H-matrix accelerated GMRES)"
        call execute_command_line(command, exitstat=exitcode, wait=.true.)
        if (exitcode /= 0) then
            simulation_execution_failed = .true.
            write(errcode, '(I0)') exitcode
            call test_failed(error, "Error: could not execute HAMS for GMRES integration test. Exit code: "//errcode)
        end if
    end subroutine run_simulation_gmres

    ! Cross-check GMRES output vs the direct-LU baseline shipped with config1.
    ! GMRES uses ITER_TOL = 1e-4 (relative), so we compare with a RELATIVE
    ! Frobenius-norm check at 1e-3 (test_oamass_files_are_close), not the
    ! absolute check used by the LU regression test — OAMASS entries span ~1e5,
    ! so an absolute tolerance would demand digit-exact agreement.
    subroutine test_oamass_gmres(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: filename_expected, filename_actual
        character(len=1) :: istr
        integer :: i

        if (simulation_execution_failed) then
            call skip_test(error, "Skipping OAMASS comparison (GMRES) because simulation failed")
            return
        end if

        do i = 1, 6
            write(istr, '(I1)') i
            filename_expected = get_repo_root_path()//baseline_dir//'expected_output/OAMASS'//istr//'.txt'
            filename_actual   = get_running_exe_path()//'output_gmres/Hams_format/OAMASS'//istr//'.txt'
            call test_oamass_files_are_close(filename_expected, filename_actual, 1.0d-3, error)
            if (allocated(error)) return
        end do
    end subroutine test_oamass_gmres

    ! Damping tolerance is looser than added mass. At the low frequencies of
    ! this test case the damping peak (~0.3) sits six orders of magnitude below
    ! the added-mass scale (~1e5), so a GMRES residual of 1e-4 relative to the
    ! full solution scale legitimately produces a few-percent relative deviation
    ! in damping. Measured deviation on this case is ~1.5e-2; 5e-2 gives
    ! headroom while still catching real solver failures (sign flips, garbage,
    ! order-of-magnitude errors).
    subroutine test_odamping_gmres(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: filename_expected, filename_actual
        character(len=1) :: istr
        integer :: i

        if (simulation_execution_failed) then
            call skip_test(error, "Skipping ODAMPING comparison (GMRES) because simulation failed")
            return
        end if

        do i = 1, 6
            write(istr, '(I1)') i
            filename_expected = get_repo_root_path()//baseline_dir//'expected_output/ODAMPING'//istr//'.txt'
            filename_actual   = get_running_exe_path()//'output_gmres/Hams_format/ODAMPING'//istr//'.txt'
            call test_oamass_files_are_close(filename_expected, filename_actual, 5.0d-2, error)
            if (allocated(error)) return
        end do
    end subroutine test_odamping_gmres

end module test_integration_gmres
