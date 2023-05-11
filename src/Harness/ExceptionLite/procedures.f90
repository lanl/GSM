! ==============================================================================
!
!> \file
!> \brief Contains the ExceptionLite procedures and private module procedures
!> \author CMJ
!
! ==============================================================================

! ==============================================================================
!
!> \fn     what
!> \brief Returns the exception message
!
!> Returns the message given to the exception. This will include the location
!> that the exception was created at if provided and applicable for the
!> compilation arguments.
!
! ARGUMENTS:
!> \param[in] ex   The exception object
!
! ==============================================================================
function what(ex)
   class(ExceptionLite), intent(in) :: ex
   return ex%msg
end function what

! ==============================================================================
!
!> \fn
!> \brief
!
!>
!
! ARGUMENTS:
!> \param[inout]
!
! ==============================================================================
!> Changes the internal exception flag to mark the exception as handled
procedure, public :: markHandled

! ==============================================================================
!
!> \fn
!> \brief
!
!>
!
! ARGUMENTS:
!> \param[inout]
!
! ==============================================================================
!> Returns a flag that indicates if the exception was handled or not
procedure, public :: handled

! ==============================================================================
!
!> \fn
!> \brief
!
!>
!
! ARGUMENTS:
!> \param[inout]
!
! ==============================================================================
!> Returns the exception ID
procedure, public :: id

! ==============================================================================
!
!> \fn
!> \brief
!
!>
!
! ARGUMENTS:
!> \param[inout]
!
! ==============================================================================
!> Returns the exception type
procedure, public :: kind

! ==============================================================================
!
!> \fn
!> \brief
!
!>
!
! ARGUMENTS:
!> \param[inout]
!
! ==============================================================================
!> Class destructor to terminate the program if exception is unhandled
final => throwUnhandled

