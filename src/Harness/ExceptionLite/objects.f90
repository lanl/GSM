! ==============================================================================
!
!> \file
!> \brief Contains the ExceptionLite module object definitions
!> \author CMJ
!
! ==============================================================================

! ==============================================================================
!
!> ExceptionLite class
!
!> Provides basic exception behaviors, specifically to create an error
!> message for a consumer to handle. Properties can only be set on object
!> construction, and only the \c handled flag may be adjusted via the \c
!> markHandled method. Otherwise, exception properties may be retrieved as
!> the consumer desires.
!
! ==============================================================================
type, public :: ExceptionLite(id, msg, exKind)
   private

   !> Indicates the ID of the exception
   integer(int16), private :: identification = 0

   !> Error message for the exception
   character(:), allocatable, private :: message

   !> Indicates the exception type
   integer(kind(GENERAL_EXCEPTION)), private :: exKind = GENERAL_EXCEPTION

   !> Flags that the exception was handled properly
   logical, private :: exceptionHandled = .false.

   contains
   private

   !> Returns the exception ID
   procedure, public :: id

   !> Returns the internal exception message
   procedure, public :: what

   !> Returns the exception type
   procedure, public :: kind

   !> Changes the internal exception flag to mark the exception as handled
   procedure, public :: markHandled

   !> Returns a flag that indicates if the exception was handled or not
   procedure, public :: handled

   !> Class destructor to terminate the program if exception is unhandled
   final => throwUnhandled

end type ExceptionLite

!> Internal exception list for exception aggregation
!
!> Exception list for internally managing and tracking thrown exceptions.
type(ExceptionLite), private, dimension(:) :: exceptionList

