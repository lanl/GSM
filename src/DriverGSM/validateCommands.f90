
  subroutine validateCommands( gsmInput, numCmdErr )

! ====================================================================
!
! Check for command argument errors, such as file not existing, etc.,
! and return the number of encountered unhandlable errors
!
!
! Written by CMJ, XCP-3 (04/2019)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32

    implicit none
    character(LEN=*), intent(in   ) :: gsmInput
    integer(int32),   intent(  out) :: numCmdErr

    logical :: fileExists

! ====================================================================

    ! Reset number of errors encountered:
    numCmdErr = 0_int32

    ! Check validity of the input file (existence):
    inquire(file=trim(gsmInput), exist=fileExists)
    if ( .not. fileExists ) then
       write(*, 1000) trim(gsmInput)
       numCmdErr = numCmdErr + 1
    end if

    return
! ====================================================================
1000 format("The file '", A, "' does not exist.")
! ====================================================================
  end subroutine validateCommands
