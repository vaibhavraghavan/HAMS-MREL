module testing_utilities

    implicit none

    public :: get_running_exe_path
    public :: get_repo_root_path
    public :: read_matrix_with_format

    ! Path to HAMS suppied by the user
    public :: hamsexe_path
    character(len=:), allocatable :: hamsexe_path
    public :: get_hams_path_from_command_line_args

    contains

    subroutine get_hams_path_from_command_line_args()
        integer :: n
        logical :: exists
        character(len=1024) :: arg

        n = command_argument_count()
        if (n >= 1) then 
            if (n > 1) then
                print *, "Warning: more than one argument supplied to test app."
            end if
            ! Get HAMS exe
            call get_command_argument(1, arg)
            hamsexe_path = trim(arg)
            ! Check file exists in path provided by ser
            inquire(file=hamsexe_path, exist=exists)
            if (.not. exists) then
                ! Check file exists in repo-root//path-by-user
                inquire(file=get_repo_root_path()//hamsexe_path, exist=exists)
                if (.not. exists) then
                    print *, "Error: HAMS executable needed for integration tests not found in: ", hamsexe_path
                    print *, "Expected usage: fpm test -- <path to hams relative to root of repo>"
                    stop 1
                end if
            else
                print*, "Integration tests: using HAMS executable in ", hamsexe_path
            end if
        else
            print*, "Error: path to HAMS executable needed for integration tests not found. Expected usage: fpm test -- <path to hams relative to root of repo>"
            stop 1
        end if
    end subroutine get_hams_path_from_command_line_args


    function get_running_exe_path() result(exepath)
        implicit none
        character(:), allocatable :: exepath
        character(len=1024) :: fullpath
        integer :: i

        ! Get path of running executable
        call get_command_argument(0, fullpath)

        ! Strip executable name (Find the position of the last '/' or '\')
        i = len_trim(fullpath)
        do while (i > 0)
            if (fullpath(i:i) == '/' .or. fullpath(i:i) == '\') exit
            i = i - 1
        end do
        exepath = trim(fullpath(1:i-1))//'/'
    end function get_running_exe_path

    function get_repo_root_path() result(reporootpath)
        implicit none
        character(:), allocatable :: reporootpath

        ! Build full path to top level of the repository,
        ! assuming exe is in repo/build/comp_id/test/
        ! which is fpm's default location for tests
        reporootpath = get_running_exe_path()//'../../../'
    end function get_repo_root_path


    ! Reads numeric data from a text file into a matrix using a user-specified Fortran format string.
    ! Arguments:
    !   fname (in)  - character string, path to input file
    !   fmt   (in)  - character string, Fortran format string for reading
    !   ncol  (in)  - integer, number of columns in the file per row
    !   A     (out) - allocatable real(kind=8) matrix
    subroutine read_matrix_with_format(filename, formatstr, ncol, A)
        implicit none
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: formatstr
        integer, intent(in) :: ncol
        real(kind=8), allocatable, intent(out) :: A(:,:)

        integer, parameter :: maxrows = 10000
        real(kind=8), allocatable :: temp(:,:)
        real(kind=8) :: row(ncol)
        integer :: ios, nrows

        ! Allocate temporary array
        allocate(temp(maxrows, ncol))
        nrows = 0

        ! Open file
        open(10, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening file:", trim(filename), "iostat=", ios
            stop 1
        end if

        ! Read file line by line into temporary array
        do
            read(10, formatstr, iostat=ios) row
            if (ios /= 0) exit  ! Reached end-of-file
            nrows = nrows + 1
            if (nrows > maxrows) then
                print *, "Error: exceeded maximum number of rows =", maxrows
                stop 1
            end if
            temp(nrows,1:ncol) = row
        end do

        close(10)

        ! Allocate exact-size array and copy data
        allocate(A(nrows,ncol))
        A(:,:) = temp(1:nrows,1:ncol)
        deallocate(temp)
    end subroutine read_matrix_with_format


    function frobenius_norm(A) result(norm)
        implicit none
        real(kind=8), intent(in) :: A(:,:)
        real :: sum, norm
        integer :: i, j

        sum = 0.0
        do i = 1, size(A, 1)
            do j = 1, size(A, 2)
                sum = sum + A(i,j)**2
            end do
        end do

        norm = sqrt(sum)
    end function frobenius_norm 


end module testing_utilities