!==============================================================================
!>
!> Contains procedure implementation for the ExceptionLite module.
!>
!==============================================================================

!==============================================================================
!
!> \fn throwUnhandled
!> \brief Halt program execution (for when consumer doesn't handle an exception)
!
!> Throws an exception if it is unhandled. This will update the global exception
!> stack as well to mark an exception as handled.
!
! ARGUMENTS:
!> \param[in] ex   An ExceptionLite class
!
!==============================================================================
   subroutine throwUnhandled(ex)
      class(ExceptionLite), intent(in) :: ex

      error stop ex%what()

   end subroutine

   !> Create an runtime exception (e.g., issues for ranges, over/underflows)
   public :: throwRunTimeException

   !> Create a logical exception (e.g., issues for arguments, ranges, etc.)
   public :: throwLogicException

   !> Creates a general exception
   public :: throwException

   !> Retrieve an exception
   public :: retrieveException

   !> Clear a single exception given a retrieved_exception
   public :: handleException

   !> Returns the ID of the next exception type and predicts the next exception
   !> ID
   private :: retrieveID

   !> Creates an exception of the appropriate type
   private :: throwExceptionInternal

