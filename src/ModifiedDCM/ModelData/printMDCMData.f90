
  subroutine printMDCMData ( lineVerbose, lineType, lineStr )

! ======================================================================
!
! Prints "lineStr" based on the verbosity of the line, assuming client
! program provides no I/O handler
!
!
! Written by CMJ, XCP-3, 8/2018
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, output_unit, error_unit

    implicit none
    integer(int32),   intent(in) :: lineVerbose   ! How 'important' the string is to print
    integer(int32),   intent(in) :: lineType      ! Comment, Warning, Error, or Fatal Error
    character(len=*), intent(in) :: lineStr       ! String to print

    integer(int32)      :: erunit, j
    character(len=1024) :: bigMessage

! ======================================================================

    ! Declare print "types"
    ! (Type index in the "typeLabel" array)
    integer(int32), parameter :: &
         & erprnt_fatal    = 1_int32,   &
         & erprnt_error    = 2_int32,   &
         & erprnt_warn     = 3_int32,   &
         & erprnt_comment  = 4_int32
    ! (Type "flag")
    character(len=*), parameter, dimension(4) :: typeLabel = [ &
         & "Fatal Error: ", &
         & "      Error: ", &
         & "    Warning: ", &
         & "    Comment: "   ]

! ======================================================================

    ! Default to follow line type passed in [Note: lineType < 0 prints strictly to error_unit]
    j = abs(lineType)
    select case( j )
       case ( erprnt_fatal   );   j = lineType
       case ( erprnt_error   );   j = lineType
       case ( erprnt_warn    );   j = lineType
       case ( erprnt_comment );   j = lineType
       case default;              j = erprnt_fatal   ! Assume fatal error for bad values of (lineType)
    end select


    if ( j == erprnt_fatal .or. mDCMDataVerbose >= lineVerbose ) then
       ! Print all fatal errors, and when the line's verbosity rank is allowed
       bigMessage = typeLabel(j) // lineStr
       if ( j == erprnt_fatal .or. j == erprnt_error .or. lineType < 0 ) then
          ! Print all fatal errors and normal errors to error unit
          erunit = error_unit
       else
          ! Print all other messages (comments and warnings) to terminal
          erunit = output_unit
       end if

       ! Print message
       write(erunit, 1000) trim(bigMessage)

    else
       ! String does not have high enough verbosity; print nothing.
    end if


    ! Reset string to empty
!    lineStr = ""

    return
! ======================================================================
1000 format(A)
! ======================================================================
  end subroutine printMDCMData
