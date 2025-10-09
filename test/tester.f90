program tester
  
    ! Testing framework imports
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use testdrive_version, only : get_testdrive_version

    ! Tests imports (add only collect functions)
    !use test_output_hams_format, only : collect_testsuite_output_hams_format
    use test_integration_config1, only : collect_tests_integration_config1
    use test_suite2, only : collect_suite2
    use testing_utilities

    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'
    integer :: major, minor, patch
    character(len=:), allocatable :: version_string

    call get_testdrive_version(major, minor, patch, version_string)
    print *, "Test-Drive version: ", version_string

    ! Collect path to HAMS executable from command line arg. Needed for integration tests
    call  get_hams_path_from_command_line_args()

    ! Specify the tests that need to be run
    testsuites = [ &
        !new_testsuite("test_output_hams_format", collect_testsuite_output_hams_format), &
        new_testsuite("test_integration_config1", collect_tests_integration_config1), &
        new_testsuite("suite2", collect_suite2) &
    ]

    stat = 0

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if

end program tester