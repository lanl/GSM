! ==============================================================================
!
!> \file
!> \brief Contains the implementation of the ExceptionLite methods
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
function what(ex) result(message)
   class(ExceptionLite), intent(in) :: ex
   character(:), allocatable :: message = ex%message
end function what

! ==============================================================================
!
!> \fn    markHandled
!> \brief Marks the exception as handled
!>
!> Marks the exception as handled. This is important because, if unhandled, the
!> exception should terminate execution at destruction if unhandled.
!
! ARGUMENTS:
!> \param[inout] The exception object
!
! ==============================================================================
subroutine markHandled(ex)
    class(ExceptionLite), intent(inout) :: ex
    ex%exceptionHandled = .true.
end subroutine markHandled

! ==============================================================================
!
!> \fn    handled
!> \brief Indicates if the exception was handled or not
!>
!> Returns a flag that indicates if the exception was handled. This is useful
!> at destruction and for clients if they were attempting to track exceptions
!> and already handled one.
!
! ARGUMENTS:
!> \param[in] The exception object
!
! ==============================================================================
function handled (ex) result(isHandled)
    class(ExceptionLite), intent(in) :: ex
    logical isHandled = ex%exceptionHandled
end function handled

! ==============================================================================
!
!> \fn    id
!> \brief Returns the identification number of the Exception
!>
!> Returns the assigned identification number of the Exception object.
!
! ARGUMENTS:
!> \param[inout] An ExceptionLite class object
!
! ==============================================================================
!> Returns the exception ID
function id(ex) result(id)
    class(ExceptionLite), intent(in) :: ex
    integer(int16) :: id = ex%identification
procedure, public :: id

! ==============================================================================
!
!> \fn    kind
!> \brief Returns a string stating what kind of exception this is.
!>
!> Indicates what kind of exception the object is, e.g., runtime, logical, etc.
!
! ARGUMENTS:
!> \param[in] The exception object
!
! ==============================================================================
function kind(ex) result(kindString)
    class(ExceptionLite), intent(in) :: ex
    character(:), allocatable :: kindString

    select case (ex%kind)
        case (GENERAL_EXCEPTION)
            kindString = "General"
        case (LOGIC_EXCEPTION)
            kindString = "Logic"
        case (RUNTIME_EXCEPTION)
            kindString = "Runtime"
        case default
            error stop "The exception type is invalid."
   end select
end function kind

