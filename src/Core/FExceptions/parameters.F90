! =============================================================================
!>
!> \file
!> \brief  Contains the FExceptions parameters
!> \author CMJ
!>
! =============================================================================

  !> \brief Logical flag to include or exclude the location of the exception
  logical, private, parameter :: includeLocation = .TRUE.

  !> \brief Character for a new line
  character(*), private, parameter :: newLine = NEW_LINE("A")

  !> \brief Text that is prepended to the Fortran Exception message
  character(*), private, parameter :: fmtFExcept = "Fortran Exception: "
  
  !> \brief Format strings for each of the exception types
  character(*), private, parameter :: &
       & fmtGeneral        = fmtFExcept // "A general exception occurred: ", &
       & fmtUnreachable    = fmtFExcept // "The code should be unreachable: ", &
       & fmtNotImplemented = fmtFExcept // "The method has not been implemented: ", &
       & fmtBadAllocation  = fmtFExcept // "A bad object allocation occurred: ", &
       & fmtUnderFlow      = fmtFExcept // "An IEEE overflow error occurred: ", &
       & fmtOverFlow       = fmtFExcept // "An IEEE underflow error occurred: ", &
       & fmtLogicException = fmtFExcept // "A logic error occured: ", &
       & fmtInvalidArg     = fmtFExcept // "The parameter supplied is not valid for ", &
       & fmtOutOfRange     = fmtFExcept // "Out of bounds for the array ", &
       & fmtRunTime        = fmtFExcept // "A general runtime error occurred: "
       ! & fmt       = fmtFExcept // "" // newLine, &

