module test_suite2
  use testdrive, only : new_unittest, unittest_type, error_type, check
  implicit none
  private

  public :: collect_suite2

contains

!> Collect all exported unit tests
subroutine collect_suite2(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("valid2", test_valid2), &
    new_unittest("invalid2", test_invalid2, should_fail=.true.) &
    ]

end subroutine collect_suite2

subroutine test_valid2(error)
  type(error_type), allocatable, intent(out) :: error
  call check(error, 1 + 2 == 3)
  if (allocated(error)) return
end subroutine test_valid2

subroutine test_invalid2(error)
  type(error_type), allocatable, intent(out) :: error
  call check(error, 1 + 2, 4)
  if (allocated(error)) return
end subroutine test_invalid2

end module test_suite2