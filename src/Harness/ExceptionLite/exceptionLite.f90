! ==============================================================================
!
!> \file
!> \brief Contains the ExceptionLite class implementation
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

   !> Error message for the exception
   character(:), allocatable, private :: msg

   !> Flags that the exception was handled properly
   logical, private :: handled = .false.

   !> Indicates the ID of the exception
   integer(int64), private :: id = 0

   !> Indicates the exception type
   integer(kind(GENERAL_EXCEPTION)), private :: exKind = GENERAL_EXCEPTION

   contains
   private

   !> Returns the internal exception message
   procedure, public :: what

   !> Changes the internal exception flag to mark the exception as handled
   procedure, public :: markHandled

   !> Returns a flag that indicates if the exception was handled or not
   procedure, public :: handled

   !> Returns the exception ID
   procedure, public :: id

   !> Returns the exception type
   procedure, public :: kind

   !> Class destructor to terminate the program if exception is unhandled
   final => throwUnhandled

end type ExceptionLite

!> TODO: Create exception list
